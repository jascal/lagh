"""A minimal transformer for modular addition mod p -- the canonical grokking idiom.

Trained to exact behaviour it computes (a+b) mod p; queried as a black box it is an
exact, real, undocumented integer oracle. An undertrained checkpoint is the negative
control (its map has errors, so no clean quasi-polynomial exists -> lagh must abstain).
"""

from __future__ import annotations

import torch
import torch.nn as nn


class TinyTransformer(nn.Module):
    def __init__(self, p: int, d: int = 64, heads: int = 4):
        super().__init__()
        self.p = p
        self.embed = nn.Embedding(p + 1, d)          # p symbols + a '=' token
        self.pos = nn.Parameter(torch.randn(3, d) * 0.02)
        layer = nn.TransformerEncoderLayer(d, heads, 4 * d, batch_first=True,
                                           dropout=0.0)
        self.enc = nn.TransformerEncoder(layer, num_layers=1)
        self.head = nn.Linear(d, p)

    def forward(self, ab: torch.Tensor) -> torch.Tensor:
        # ab: (n, 2) integer pairs; append the '=' token
        eq = torch.full((ab.shape[0], 1), self.p, dtype=torch.long)
        x = torch.cat([ab, eq], dim=1)
        h = self.embed(x) + self.pos[None, :, :]
        h = self.enc(h)
        return self.head(h[:, -1, :])                # predict from the '=' position


def train(p: int, steps: int, seed: int = 0) -> TinyTransformer:
    torch.manual_seed(seed)
    a, b = torch.meshgrid(torch.arange(p), torch.arange(p), indexing="ij")
    X = torch.stack([a.reshape(-1), b.reshape(-1)], dim=1)
    y = (X[:, 0] + X[:, 1]) % p
    model = TinyTransformer(p)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1.0)
    lossf = nn.CrossEntropyLoss()
    for _ in range(steps):
        opt.zero_grad()
        loss = lossf(model(X), y)
        loss.backward()
        opt.step()
    model.eval()
    with torch.no_grad():
        acc = (model(X).argmax(1) == y).float().mean().item()
    model._train_acc = acc
    return model


def oracle_fn(model: TinyTransformer, b_fixed: int):
    """a (integers, any range) -> model's predicted (a + b_fixed) mod p, argmax."""
    def oracle(A):
        import numpy as np
        A = np.atleast_2d(np.asarray(A))
        a = torch.tensor([int(round(float(r[0]))) % model.p for r in A],
                         dtype=torch.long)
        b = torch.full_like(a, b_fixed % model.p)
        with torch.no_grad():
            pred = model(torch.stack([a, b], dim=1)).argmax(1)
        return pred.numpy().astype(float)
    return oracle


if __name__ == "__main__":
    for p, steps in [(7, 4000)]:
        m = train(p, steps)
        print(f"p={p}: train/exact accuracy = {m._train_acc:.4f}")
