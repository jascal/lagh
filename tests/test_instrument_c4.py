"""C4: known causes, refusal controls and independent covariance truth."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest
from scipy.signal import lfilter

from experiments.tweezers.adapter import bfp_diode_at_power
from experiments.tweezers.run_c4 import covariance_measurement
from lagh.instrument import Retention, attribute_deviation


def test_slow_signature_predicts_diffusion_and_can_still_refuse():
    a = attribute_deviation(.019, 1.03, 3.797)
    assert a['verdict'] == 'slow-contaminant'
    assert a['hypotheses']['slow-contaminant']['value'] == pytest.approx(2.797)
    assert a['hypotheses']['slow-contaminant']['max_residual'] == pytest.approx(.03)
    assert attribute_deviation(.019, .885, 3.797)['verdict'] == 'unattributed'
    assert attribute_deviation(.019, 1.2, 3.797)['verdict'] == 'unattributed'
    assert attribute_deviation(.55, 1.02, 1.038)['verdict'] == 'unattributed'


@pytest.mark.parametrize('ratios,verdict', [
    ((1, 1, 1), 'consistent'), ((.7, .7, 1), 'drag'),
    ((1, 1.44, 1.44), 'scale'), ((.7, 1, 1/.7), 'unattributed'),
    ((.7, 1.4, .3), 'unattributed'), ((0, 1, 1), 'unreadable')])
def test_existing_signatures_and_ambiguity(ratios, verdict):
    assert attribute_deviation(*ratios)['verdict'] == verdict


def test_diode_uses_diagnostic_power_not_calibration_power():
    cal = {'Diode alpha delta': .245, 'Diode alpha max': .629,
           'Diode alpha rate': .85, 'Diode frequency delta': 2619,
           'Diode frequency max': 16607, 'Diode frequency rate': 1.07}
    # Producer tutorial's independent example at this recorded power.
    current = bfp_diode_at_power(cal, [.3622186360371995])
    assert current['alpha'] == pytest.approx(.4489251910346808)
    assert current['f_diode (Hz)'] == pytest.approx(14829.480905511606)
    old = bfp_diode_at_power(cal, [.11839534279659139])
    assert old['alpha'] == pytest.approx(.4074557886886479)
    with pytest.raises(ValueError):
        bfp_diode_at_power(cal, [np.nan])


def simulate_ou(n, fs, theta, variance, rng):
    decay = np.exp(-theta/fs)
    return lfilter([np.sqrt(variance*(1-decay**2))], [1, -decay],
                   rng.normal(size=n))


@pytest.mark.parametrize('kind', ['slow-ou', 'slow-line'])
def test_covariance_subtraction_measures_amplitude_without_assuming_decay_shape(kind):
    fs, theta, n = 10000., 300., 800000
    rng = np.random.default_rng(408)
    thermal = simulate_ou(n, fs, theta, 1., rng)
    tau = .04
    if kind == 'slow-ou':
        contaminant = simulate_ou(n, fs, 1/tau, 3., rng)
        expected_crossing = tau
    else:
        contaminant = np.sqrt(6)*np.sin(2*np.pi*5*np.arange(n)/fs)
        expected_crossing = np.arccos(1/np.e)/(2*np.pi*5)
    fc = theta/(2*np.pi)
    fr = np.linspace(0, fs/2, 32769)
    amplitude = 2*fc/np.pi  # known thermal variance=1, no data fit
    model = amplitude/(fc**2+fr**2)
    ret = Retention(fr, model, np.ones_like(fr), np.ones_like(fr), fs,
                    amplitude, fc, (100,2300), 1., 10000.)
    m = covariance_measurement(thermal+contaminant, fs, ret)
    assert abs(m['recorded_excess_variance_V2']/3-1) < .15
    assert abs(m['first_e_crossing_s']/expected_crossing-1) < .2
    if kind == 'slow-line':
        assert m['min_normalized_covariance'] < -.85
    assert 'theta' not in m and 'certified' not in m


def test_independent_slow_process_attributes_from_actual_observations():
    from lagh.instrument import acf_timescale, realized_diffusion
    fs, theta, n = 10000., 500., 800000
    rng = np.random.default_rng(409)
    x = (simulate_ou(n, fs, theta, 1., rng)
         + simulate_ou(n, fs, 2., 3., rng))
    fr = np.linspace(0, fs/2, 32769)
    fc, amplitude = theta/(2*np.pi), theta/np.pi**2
    model = amplitude/(fc**2+fr**2)
    # Exact sampled OU PSD includes aliases; no Nyquist-truncation fiction.
    decay = np.exp(-theta/fs)
    sampled_psd = 2/fs*(1-decay**2)/(1+decay**2-2*decay*np.cos(2*np.pi*fr/fs))
    ret = Retention(fr, model, np.ones_like(fr), sampled_psd/model, fs,
                    amplitude, fc, (100,2300), 1., 10000.)
    apparent = acf_timescale(x, 1/fs, min_lag=8, max_lag=200)['theta']
    d = realized_diffusion(x, 1/fs, theta, ret, strides=[2, 3, 4, 6])
    ratios = (apparent/theta, d['plateau']/(2*theta), float(x.var()))
    assert attribute_deviation(*ratios)['verdict'] == 'slow-contaminant'
    assert abs(ratios[1]-1) < .08
