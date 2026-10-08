"""The boundary audit rechecks evidence instead of trusting saved success flags."""
import json

import numpy as np
import pytest

from experiments.run_certification_boundaries import FAMILIES, REGIMES, TRUTHS, data, digest
from experiments.score_certification_boundaries import read_arm
from lagh.certify import claim_scope

CASE = "mono3__noise6__710"


def synthetic_arm(tmp_path):
    """Exercise the audit on current-runtime inputs, not cross-CPU float replay.

    These are unit-test records, not additional scored study attempts. The
    historical manifests and the audit's exact hash comparison stay untouched.
    """
    target = tmp_path/"arm"
    target.mkdir()
    for family in FAMILIES:
        for regime in REGIMES:
            for seed in (710, 711):
                name = f"{family}__{regime}__{seed}"
                X, y, _, _ = data(family, regime, seed)
                row = {"case": name, "family": family, "regime": regime,
                       "seed": seed, "status": "refused",
                       "input_hash": digest(np.column_stack([X, y]))}
                if name == CASE:
                    row.update(status="certified", law=str(TRUTHS[family]),
                        certificate={"claim": claim_scope(y, sigma=REGIMES[regime][0])},
                        score={"observed_band_ok": True, "interval_diagnostics": {
                            "slices": [{"atom": "13/10", "bounds": [1.3, 1.3]}]}})
                (target/(name+".json")).write_text(json.dumps(row))
    # Each mutation starts with a valid fixture; an always-rejecting audit fails.
    _, report = read_arm(target, (710, 711), require_scope=True)
    assert report["issues"] == []
    path = target/(CASE+".json")
    return target, path, json.loads(path.read_text())


def test_missing_scored_case_cannot_be_hidden_as_a_refusal(tmp_path):
    folder, path, _ = synthetic_arm(tmp_path)
    path.unlink()
    with pytest.raises(ValueError, match="incomplete/extra"):
        read_arm(folder, (710, 711), require_scope=True)


def test_audit_rechecks_law_even_when_all_recorded_checks_claim_success(tmp_path):
    folder, path, row = synthetic_arm(tmp_path)
    row["law"] = "x_0 + 100"
    path.write_text(json.dumps(row))
    _, report = read_arm(folder, (710, 711), require_scope=True)
    assert [CASE, "returned law fails supplied observations"] in report["issues"]


def test_uncertain_fit_cannot_be_relabelled_as_exact_identification(tmp_path):
    folder, path, row = synthetic_arm(tmp_path)
    row["certificate"]["claim"]["exact_form_identified"] = True
    path.write_text(json.dumps(row))
    _, report = read_arm(folder, (710, 711), require_scope=True)
    assert [CASE, "claim absent or stronger than evidence"] in report["issues"]


def test_changed_input_hash_is_rejected(tmp_path):
    folder, path, row = synthetic_arm(tmp_path)
    row["input_hash"] = "0"*64
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError, match="input hash differs"):
        read_arm(folder, (710, 711), require_scope=True)


def test_audit_rechecks_saved_interval_endpoints(tmp_path):
    folder, path, row = synthetic_arm(tmp_path)
    row["score"]["interval_diagnostics"]["slices"][0]["bounds"][1] = 100.
    path.write_text(json.dumps(row))
    _, report = read_arm(folder, (710, 711), require_scope=True)
    assert any(case == CASE and "infeasible slice endpoint" in issue
               for case, issue in report["issues"])
