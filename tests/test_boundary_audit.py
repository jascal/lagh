"""The boundary audit rechecks evidence instead of trusting saved success flags."""
import json
from pathlib import Path
import shutil

import pytest

from experiments.score_certification_boundaries import read_arm

ARM = Path(__file__).resolve().parents[1]/"experiments/results/certification_boundaries/repair"
CASE = "mono3__noise6__710"


def copied_arm(tmp_path):
    target = tmp_path/"arm"
    shutil.copytree(ARM, target)
    path = target/(CASE+".json")
    return target, path, json.loads(path.read_text())


def test_missing_scored_case_cannot_be_hidden_as_a_refusal(tmp_path):
    folder, path, _ = copied_arm(tmp_path)
    path.unlink()
    with pytest.raises(ValueError, match="incomplete/extra"):
        read_arm(folder, (710, 711), require_scope=True)


def test_audit_rechecks_law_even_when_all_recorded_checks_claim_success(tmp_path):
    folder, path, row = copied_arm(tmp_path)
    row["law"] = "x_0 + 100"
    path.write_text(json.dumps(row))
    _, report = read_arm(folder, (710, 711), require_scope=True)
    assert [CASE, "returned law fails supplied observations"] in report["issues"]


def test_noise_mismatch_cannot_be_relabelled_as_exact_identification(tmp_path):
    folder, path, row = copied_arm(tmp_path)
    row["certificate"]["claim"]["exact_form_identified"] = True
    path.write_text(json.dumps(row))
    _, report = read_arm(folder, (710, 711), require_scope=True)
    assert [CASE, "claim absent or stronger than evidence"] in report["issues"]
