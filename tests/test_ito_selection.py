"""Failing inputs executed first in selection_artifact_witnesses.json."""
import numpy as np
import pytest
from lagh.ito_selection import choose_test
import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "compare_artifacts", Path(__file__).resolve().parents[1]/"experiments/compare_artifacts.py")
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
differences = _module.differences


def test_selector_rejects_no_evidence_and_invalid_grid():
    with pytest.raises(ValueError):
        choose_test(np.arange(10.),np.ones((2,10)))
    with pytest.raises(ValueError):
        choose_test(np.arange(10.)**2,np.tile(np.arange(10.),(2,1)))


def test_artifact_tolerance_has_a_failing_input():
    assert differences({'a':1.},{'a':1.00001})
    assert not differences({'a':1.},{'a':float(np.nextafter(1.,2.))})
    assert differences({'a':1.},{'a':1})
    assert differences({'a':float('nan')},{'a':float('nan')})
