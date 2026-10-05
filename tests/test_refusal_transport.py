"""Refusal transport witnesses; fake splits test plumbing, not discovery reach.

Falsifier: a structural first split followed by a rival-free split. The former
consumer returned only the last result and erased the evidence driving refusal.
"""
import numpy as np
import pytest
import sympy as sp

import lagh.mcp.core as core
import lagh.passive as passive
from lagh.acquisition import ActiveResult, Ledger
from lagh.certify import Certificate
from lagh.engine import Result

x = sp.Symbol('x_0')
X = np.linspace(1, 2, 40)[:, None]
Y = X[:, 0]


def refusal(rivals=()):
    return Result(Certificate(False, 0, 0, 8, [[1, 2]], '',
                              abstain='structural',
                              notes=['materially different classes']),
                  None, 1, 2, tuple(rivals))


def accepted():
    return Result(Certificate(True, 0, 0, 8, [[1, 2]], str(x)), x, 1, 1)


@pytest.mark.parametrize('last_passes', [False, True])
def test_rivals_survive_later_split_without_rivals(monkeypatch, last_passes):
    # The accepted later split triggers sticky demotion; neither exit may
    # discard the earlier pair, even though its own Result has no rivals.
    outcomes = iter([refusal((x, x + 1)), accepted() if last_passes else refusal()])
    monkeypatch.setattr(passive, 'discover', lambda *a, **kw: next(outcomes))
    out = passive.discover_passive(X, Y, n_resplits=2)
    assert not out.certified
    assert out.result.rivals == (x, x + 1)


def test_all_split_representatives_preserved_without_duplicate(monkeypatch):
    # Last-only retention loses x; naive concatenation duplicates x+1.
    outcomes = iter([refusal((x, x + 1)), refusal((x + 1, x + 2))])
    monkeypatch.setattr(passive, 'discover', lambda *a, **kw: next(outcomes))
    out = passive.discover_passive(X, Y, n_resplits=2)
    assert out.result.rivals == (x, x + 1, x + 2)


def test_mcp_passive_exposes_design_evidence_keeps_fixed_box(monkeypatch):
    # A rival-bearing refusal previously serialized only a class-count note.
    monkeypatch.setattr(core, 'discover_passive', lambda *a, **kw:
                        passive.PassiveResult(refusal((x, x + 1)), 1, None))
    monkeypatch.setattr(core, 'characterize', lambda *a, **kw:
                        {'research': {'move': 'acquire'}})
    out = core.recover(X, Y)
    assert out['tag'] == 'open' and not out['certified']
    assert 'law' not in out
    assert out['suggested_box'] == [[0.1], [20.0]]
    evidence = out['design_evidence']
    assert evidence['rivals'] == ['x_0', 'x_0 + 1']
    assert evidence['tag'] == 'empirical'
    assert evidence['requires_fresh_certification'] is True


def test_mcp_active_exposes_design_evidence(monkeypatch):
    # Both API paths must carry content, without turning rivals into laws.
    box = np.array([[1.0], [2.0]])
    active = ActiveResult(refusal((x, x + 1)), box, box, [], Ledger(), [], 0)
    monkeypatch.setattr(core, 'run_active', lambda *a, **kw: active)
    monkeypatch.setattr(core, '_characterize_oracle', lambda *a, **kw: None)
    out = core.recover(oracle=lambda _: pytest.fail('no oracle expected'),
                       box=box.tolist())
    assert out['tag'] == 'open' and 'law' not in out
    assert out['design_evidence']['rivals'] == ['x_0', 'x_0 + 1']
