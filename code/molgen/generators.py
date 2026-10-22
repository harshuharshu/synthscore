# -*- coding: utf-8 -*-
"""
generators.py - pluggable generation backends.

Two sources, same interface: each yields candidate SMILES strings.

  1. BRICSGenerator    - BRICS fragment recombination (RDKit). No training; the
                         output is a recombination of the supplied fragment
                         library, so its chemistry is bounded by that library.
                         Note that BRICS recombination over-produces chained
                         polycyclic systems relative to approved-drug space
                         (measured 14.0% vs 6.2%; see AUDIT.md).

  2. CharRNNGenerator  - a SMILES character-level RNN prior (PyTorch), trained
                         here on a SMILES set the caller supplies, then sampled.
                         Small enough to train on a single consumer GPU.

Reproducibility: both backends take a seed. BRICSGenerator uses its own
random.Random instance; CharRNNGenerator seeds torch, CUDA and its own shuffle
RNG (Defect 6.2 in AUDIT.md was that neither was seeded and training shuffled
the global RNG). Determinism still depends on the cuDNN backend; pass
deterministic=True to constrain it at some cost in speed.

Stereochemistry is retained in emitted SMILES (isomericSmiles=True). The original
code emitted flat SMILES while the property panel counted stereocentres
(Defect 6.1).
"""

import random

import torch
import torch.nn as nn
import torch.nn.functional as F

from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS

RDLogger.DisableLog("rdApp.*")


# ===========================================================================
# BRICS backend
# ===========================================================================
class BRICSGenerator:
    def __init__(self, fragments, n_sample=4, seed=42):
        self.fragments = [f for f in fragments if Chem.MolFromSmiles(f) is not None]
        self.n_sample = n_sample
        self.rng = random.Random(seed)
        if len(self.fragments) < n_sample:
            raise ValueError(f"Need >= {n_sample} valid fragments, have {len(self.fragments)}")

    def sample(self, max_yield):
        """Yield up to max_yield canonical SMILES (caller dedupes/filters)."""
        produced = 0
        guard = 0
        while produced < max_yield and guard < max_yield * 50:
            guard += 1
            chosen = self.rng.sample(self.fragments, self.n_sample)
            mols = [Chem.MolFromSmiles(f) for f in chosen]
            mols = [m for m in mols if m is not None]
            if len(mols) < 2:
                continue
            try:
                for prod in BRICS.BRICSBuild(mols):
                    prod.UpdatePropertyCache(strict=False)
                    smi = Chem.MolToSmiles(prod, isomericSmiles=True, canonical=True)
                    if Chem.MolFromSmiles(smi) is not None:
                        produced += 1
                        yield smi
                        if produced >= max_yield:
                            return
            except Exception:
                continue


# ===========================================================================
# Char-RNN backend
# ===========================================================================
PAD, BOS, EOS = "<pad>", "^", "$"


class SmilesVocab:
    def __init__(self, smiles_list):
        chars = set()
        for s in smiles_list:
            chars.update(s)
        self.itos = [PAD, BOS, EOS] + sorted(chars)
        self.stoi = {c: i for i, c in enumerate(self.itos)}

    def __len__(self):
        return len(self.itos)

    def encode(self, smi):
        return [self.stoi[BOS]] + [self.stoi[c] for c in smi] + [self.stoi[EOS]]

    def decode(self, idxs):
        out = []
        for i in idxs:
            c = self.itos[i]
            if c == EOS:
                break
            if c in (BOS, PAD):
                continue
            out.append(c)
        return "".join(out)


class CharRNN(nn.Module):
    def __init__(self, vocab_size, emb=128, hidden=512, layers=3, dropout=0.2):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, emb, padding_idx=0)
        self.lstm = nn.LSTM(emb, hidden, layers, batch_first=True,
                            dropout=dropout if layers > 1 else 0.0)
        self.fc = nn.Linear(hidden, vocab_size)

    def forward(self, x, hidden=None):
        e = self.embed(x)
        out, hidden = self.lstm(e, hidden)
        return self.fc(out), hidden


class CharRNNGenerator:
    """Train a SMILES prior once, then sample. Persist via save()/load()."""

    def __init__(self, vocab=None, model=None, device=None, seed=42,
                 deterministic=False):
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.vocab = vocab
        self.model = model.to(self.device) if model is not None else None
        self.seed = seed
        self.rng = random.Random(seed)
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        if deterministic:
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False

    # --- training ---
    def fit(self, smiles_list, epochs=20, batch_size=128, lr=1e-3,
            max_len=120, hidden=512, layers=3, log=print):
        # reseed here too: weight initialization happens below, and callers may
        # have consumed the RNG between __init__ and fit().
        torch.manual_seed(self.seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(self.seed)
        self.rng = random.Random(self.seed)

        smiles_list = [s for s in smiles_list
                       if s and len(s) <= max_len and Chem.MolFromSmiles(s) is not None]
        if len(smiles_list) < 50:
            raise ValueError(f"Need a reasonable training set; got {len(smiles_list)} valid SMILES.")
        self.vocab = SmilesVocab(smiles_list)
        self.model = CharRNN(len(self.vocab), hidden=hidden, layers=layers).to(self.device)

        seqs = [torch.tensor(self.vocab.encode(s), dtype=torch.long) for s in smiles_list]
        opt = torch.optim.Adam(self.model.parameters(), lr=lr)
        pad_idx = self.vocab.stoi[PAD]

        self.model.train()
        for ep in range(1, epochs + 1):
            self.rng.shuffle(seqs)
            total, n = 0.0, 0
            for i in range(0, len(seqs), batch_size):
                batch = seqs[i:i + batch_size]
                lengths = [len(s) for s in batch]
                maxl = max(lengths)
                inp = torch.full((len(batch), maxl), pad_idx, dtype=torch.long)
                for j, s in enumerate(batch):
                    inp[j, :len(s)] = s
                inp = inp.to(self.device)
                x, y = inp[:, :-1], inp[:, 1:]
                logits, _ = self.model(x)
                loss = F.cross_entropy(
                    logits.reshape(-1, logits.size(-1)), y.reshape(-1),
                    ignore_index=pad_idx)
                opt.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 5.0)
                opt.step()
                total += loss.item() * len(batch)
                n += len(batch)
            log(f"  epoch {ep:3d}/{epochs}  loss={total / n:.4f}")
        return self

    # --- sampling ---
    @torch.no_grad()
    def sample(self, max_yield, temperature=1.0, max_len=120, batch=256):
        """Yield valid canonical SMILES sampled from the prior."""
        self.model.eval()
        produced = 0
        bos, eos, pad = self.vocab.stoi[BOS], self.vocab.stoi[EOS], self.vocab.stoi[PAD]
        guard = 0
        while produced < max_yield and guard < max_yield * 20:
            guard += batch
            x = torch.full((batch, 1), bos, dtype=torch.long, device=self.device)
            hidden = None
            done = torch.zeros(batch, dtype=torch.bool, device=self.device)
            seqs = [[] for _ in range(batch)]
            for _ in range(max_len):
                logits, hidden = self.model(x, hidden)
                logits = logits[:, -1, :] / max(temperature, 1e-3)
                probs = F.softmax(logits, dim=-1)
                nxt = torch.multinomial(probs, 1).squeeze(1)
                for b in range(batch):
                    if not done[b]:
                        if nxt[b].item() == eos:
                            done[b] = True
                        else:
                            seqs[b].append(nxt[b].item())
                x = nxt.unsqueeze(1)
                if done.all():
                    break
            for b in range(batch):
                smi = self.vocab.decode(seqs[b])
                m = Chem.MolFromSmiles(smi)
                if m is not None:
                    can = Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
                    produced += 1
                    yield can
                    if produced >= max_yield:
                        return

    # --- persistence ---
    def save(self, path):
        torch.save({"itos": self.vocab.itos,
                    "state_dict": self.model.state_dict(),
                    "hidden": self.model.lstm.hidden_size,
                    "layers": self.model.lstm.num_layers}, path)

    @classmethod
    def load(cls, path, device=None):
        dev = device or ("cuda" if torch.cuda.is_available() else "cpu")
        ckpt = torch.load(path, map_location=dev)
        vocab = SmilesVocab.__new__(SmilesVocab)
        vocab.itos = ckpt["itos"]
        vocab.stoi = {c: i for i, c in enumerate(vocab.itos)}
        model = CharRNN(len(vocab), hidden=ckpt["hidden"], layers=ckpt["layers"])
        model.load_state_dict(ckpt["state_dict"])
        return cls(vocab=vocab, model=model, device=dev)
