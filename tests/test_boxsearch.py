"""Box-search's held-out guard (jascal/lagh#2): a result the independent box
REJECTS must leave the search demoted, on every exit path, and `recover` must
never translate it into `tag: proved`."""

import numpy as np

import lagh.acquisition as acq
from lagh.acquisition import run_active_boxsearch
from lagh.mcp.core import recover

ORACLE = lambda X: 3 * X[:, 0] ** 2                            # noqa: E731


def test_holdout_failure_on_the_last_box_demotes_the_result(monkeypatch):
    monkeypatch.setattr(acq, "_heldout_box_ok", lambda *a, **k: False)
    bs = run_active_boxsearch(ORACLE, [0.5], [4.0], max_boxes=1)
    c = bs.active.result.certificate
    assert bs.heldout_box_ok is False and bs.boxes_tried == 1
    assert c.certified is False and c.abstain == "heldout"
    assert bs.active.result.expr is None
    assert any("held-out box guard" in n for n in c.notes)


def test_holdout_failure_on_the_time_budget_exit_demotes_the_result(monkeypatch):
    class Clock:
        t = 0.0

        def time(self):
            return self.t

    clock = Clock()
    monkeypatch.setattr(acq, "_time", clock)

    def failing_holdout(*a, **k):
        clock.t = 1e9              # the budget expires right after the failure
        return False

    monkeypatch.setattr(acq, "_heldout_box_ok", failing_holdout)
    bs = run_active_boxsearch(ORACLE, [0.5], [4.0], time_budget_s=100.0)
    assert bs.boxes_tried == 1 and bs.heldout_box_ok is False
    c = bs.active.result.certificate
    assert c.certified is False and c.abstain == "heldout"


def test_recover_never_emits_proved_after_a_failed_holdout(monkeypatch):
    monkeypatch.setattr(acq, "_heldout_box_ok", lambda *a, **k: False)
    monkeypatch.setattr(acq, "_box_ladder", lambda lo, hi: iter(
        [("as-declared", np.asarray(lo, float), np.asarray(hi, float))]))
    r = recover(oracle=ORACLE, box=[[0.5], [4.0]], box_search=True)
    assert r["tag"] == "open" and r["certified"] is False
    assert r["abstain"] == "heldout"
    assert r["acquisition"]["heldout_box_ok"] is False
    assert "law" not in r
