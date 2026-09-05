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


@pytest.mark.parametrize('theta', [.999, 1., 1.001])
def test_free_theta_does_not_gate_slow_signature(theta):
    a = attribute_deviation(theta, 1., 3.)
    assert a['verdict'] == 'slow-contaminant'
    assert a['slow_contaminant_admission']['theta_used_for_admission'] is False


def test_materiality_is_reported_and_independent_of_residual_tolerance():
    a = attribute_deviation(.5, 1., 1.079, tol=.04)
    b = attribute_deviation(.5, 1., 1.081, tol=.04)
    assert not a['slow_contaminant_admission']['admitted']
    assert b['slow_contaminant_admission']['admitted']
    assert a['slow_contaminant_admission']['variance_floor'] == .08
    assert a['slow_contaminant_admission']['distance_to_floor'] < 0
    assert b['slow_contaminant_admission']['distance_to_floor'] > 0
    assert attribute_deviation(.5, 1., 1.06, contaminant_variance_floor=.02)['verdict'] == 'slow-contaminant'


def test_stiffness_increase_still_has_coverage():
    assert attribute_deviation(1.4, 1., 1/1.4)['verdict'] == 'stiffness'


def test_noise_training_repeatability_cannot_vote_in_validation():
    from experiments.tweezers.run_c4 import floor_comparison
    fr = np.arange(0, 50001, 100.)
    thermal = np.ones_like(fr)*1e-13
    train, held = thermal+2e-13, thermal+2e-13
    held[(fr >= 30000) & (fr < 40000)] += 1e-12
    # Simulate a stopband as well: neither it nor the training band may vote.
    held[fr >= 43000] = 0
    n = floor_comparison(fr, train, held, thermal)
    assert n['passband_model_passed']
    assert n['training_band_repeatability']['excess_over_train_floor'] > 2
    assert [b['band_hz'] for b in n['passband_validation']] == [[20000.,30000.], [40000.,43000.]]
    assert n['sigma_obs_V'] is None and n['detector_whiteness'] == 'unresolved'
    held[(fr >= 20000) & (fr < 30000)] += 1e-12
    assert not floor_comparison(fr, train, held, thermal)['passband_model_passed']


def test_unresolved_subtractions_are_counted_in_sensitivity():
    from experiments.tweezers.run_c4 import sensitivity_summary
    r = sensitivity_summary([{'rms_V': 2.}, {'verdict': 'unresolved', 'reason': 'negative excess'}])
    assert r['rms_V_range'] == [2., 2.] and r['n_unresolved'] == 1
    assert r['verdict'] == 'partly-unresolved'
    assert sensitivity_summary([{'reason': 'fit failed'}])['rms_V_range'] is None


def test_mean_removal_matches_the_thermal_estimator():
    fs, n = 100., 32
    fr = np.linspace(0, fs/2, 33)
    ret = Retention(fr, np.full_like(fr, 2/fs), np.ones_like(fr),
                    np.ones_like(fr), fs, 1., 1., (1., 20.), 1., 10.)
    x = np.tile([-1., 1.], n//2)
    a = covariance_measurement(x, fs, ret)
    b = covariance_measurement(x+100., fs, ret)
    assert a['thermal_mean_removal_V2'] == pytest.approx(1/n)
    assert a['thermal_retained_variance_V2'] == pytest.approx(1-1/n)
    assert a['recorded_excess_variance_V2'] == b['recorded_excess_variance_V2']
    assert a['population_variance'] is None
    assert 'unknown' in a['uncertainty']


@pytest.mark.parametrize('failure', ['acf', 'stride', 'fit'])
def test_axis_gate_survives_unresolved_paths(monkeypatch, failure):
    from experiments.tweezers import run_c4 as campaign
    x = np.arange(128.)
    de = {'D (V^2/s)': 1., 'fc (Hz)': 1.}
    fr = np.linspace(0, 50, 65)
    ret = Retention(fr, np.ones_like(fr), np.ones_like(fr), np.ones_like(fr),
                    100., 1., 1e9 if failure == 'stride' else 1., (1,20), 1., 10.)
    monkeypatch.setattr(campaign, 'read_axis', lambda *_: (x, de, {'fs_hz': 100.}))
    def fake_fit(*args):
        if failure == 'fit':
            raise ValueError('no spectral bins')
        return ret
    monkeypatch.setattr(campaign, 'fit', fake_fit)
    calls = []
    def gate(*args, **kwargs):
        calls.append(kwargs)
        return {'passed': False, 'theta_acf': None if failure == 'acf' else 1.,
                'ratio': None if failure == 'acf' else 1.}
    monkeypatch.setattr(campaign, 'axis_gate', gate)
    r = campaign.contaminated_axis(None, 'Force 1x')
    assert len(calls) == 2
    assert all(h['law_certificate'] is False and 'reason' in h for h in r['halves'])


def test_diode_provenance_is_checked_with_python_optimization():
    import subprocess
    script = '''
from experiments.tweezers.run_c4 import validate_diode_model
cal = {'Trap sum power (V)': 0., 'alpha': .9, 'f_diode (Hz)': 10.,
       'Diode alpha max': .5, 'Diode alpha delta': .1, 'Diode alpha rate': 1.,
       'Diode frequency max': 11., 'Diode frequency delta': 1., 'Diode frequency rate': 1.}
try:
    validate_diode_model(cal)
except ValueError:
    print('refused')
else:
    raise RuntimeError('provenance check vanished')
'''
    result = subprocess.run([sys.executable, '-O', '-c', script],
                            cwd=Path(__file__).resolve().parents[1],
                            capture_output=True, text=True, check=True)
    assert result.stdout.strip() == 'refused'


def test_plot_keeps_unresolved_covariance_visible(tmp_path):
    from experiments.tweezers.plot_c4 import plot
    r = {'passive': {'Force 1x': {'halves': [{}, {'reason': 'no ACF'}]}},
         'noise_floor': {'Force 1x': {'verdict': 'unresolved'}}}
    path = tmp_path/'unresolved.png'
    plot(r, path)
    assert path.stat().st_size > 0


def test_voltage_adapter_uses_applied_response_directly(monkeypatch):
    from experiments.tweezers import adapter
    class Dataset:
        attrs = {'Sample rate (Hz)': 10., 'Start time (ns)': 100, 'Stop time (ns)': 200}
        def __getitem__(self, key):
            return np.array([2., 4., 6.])
    cal = {'conversion_start': 0, 'voltage_start': 100,
           'Rf_transform': 2., 'Rd (um/V)': 3., 'item': 'test'}
    monkeypatch.setattr(adapter, 'bfp_calibrations', lambda *_: [cal])
    f = {'Force HF/Force 1x': Dataset()}
    _, volts, _, _, _ = adapter.bfp_voltage(f, 'Force 1x')
    _, nm, _, _, _ = adapter.bfp_position_nm(f, 'Force 1x')
    np.testing.assert_array_equal(volts, [1., 2., 3.])
    np.testing.assert_array_equal(nm, [3000., 6000., 9000.])
