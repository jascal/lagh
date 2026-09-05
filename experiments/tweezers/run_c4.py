"""Registered C4: additive contamination and held-out detector-noise tests."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.fft import next_fast_len, rfft, irfft
from scipy.signal import welch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.tweezers import adapter  # noqa: E402
from lagh.instrument import (axis_gate, attribute_deviation,  # noqa: E402
                             realized_diffusion, retention)

OUT = Path(__file__).resolve().parents[1] / 'results/tweezers_c4.json'
STRIDES = (2, 3, 4, 6, 8)


def read_axis(f, channel):
    _, volts, applied, de, meta = adapter.bfp_voltage(f, channel)
    return volts, de, {**meta,
        'applied': applied, 'reference_calibration': de,
        'record_start_ns': int(f[f'Force HF/{channel}'].attrs['Start time (ns)'])}


def fit(x, fs, de, fit_range=None):
    return retention(x, fs, alpha=de['alpha'], f_diode=de['f_diode (Hz)'],
                     fit_range=fit_range or (de['Fit range (min.) (Hz)'],
                                            de['Fit range (max.) (Hz)']))


def observed_covariance(x, fs):
    """Sample-mean-removed covariance, cached for repeated thermal subtractions."""
    x = np.asarray(x, float)
    if x.ndim != 1 or len(x) < 8 or not np.all(np.isfinite(x)) or fs <= 0:
        raise ValueError('finite one-dimensional record of at least eight samples required')
    x = x - x.mean()
    n = len(x)
    nfft = next_fast_len(2*n)
    maxlag = min(int(.5*fs), n//4)
    power = np.abs(rfft(x, nfft))**2
    cov = irfft(power, nfft)[:maxlag+1] / (n-np.arange(maxlag+1))
    return {'cov': cov, 'n': n, 'nfft': nfft, 'fs': fs, 'maxlag': maxlag}


def covariance_measurement(x, fs, thermal, *, observed=None, apply_measured=True):
    """Within-window statistic, NOT an estimate of stationary population variance.

    Match thermal-model mean removal to the empirical estimator. The unobserved
    slow component's removed sample mean cannot be corrected without a model;
    its population uncertainty remains unknown. Observed FFT is reused.
    """
    if observed is None:
        try:
            observed = observed_covariance(x, fs)
        except ValueError as exc:
            return {'verdict': 'unresolved', 'reason': str(exc)}
    n, nfft, maxlag = (observed[k] for k in ('n', 'nfft', 'maxlag'))
    cov = observed['cov']
    fr = np.fft.rfftfreq(nfft, 1/fs)
    response = thermal.measured if apply_measured else 1
    psd = np.interp(fr, thermal.freqs, thermal.model*thermal.diode*response)
    # Exact expectation of the same demeaning operation for the supplied
    # stationary thermal covariance, including edge-dependent sample means.
    tcov = irfft(psd, nfft)[:n] * fs/2
    csum = np.r_[0., np.cumsum(tcov[1:])]
    row_mean_cov = (tcov[0] + csum + csum[::-1])/n
    prefix = np.r_[0., np.cumsum(row_mean_cov)]
    lags = np.arange(maxlag+1)
    tcentered = (tcov[:maxlag+1]
                 - (prefix[n-lags]+prefix[n]-prefix[lags])/(n-lags)
                 + row_mean_cov.mean())
    excess = cov-tcentered
    out = {'evidence': 'empirical', 'recorded_excess_variance_V2': float(excess[0]),
           'thermal_retained_variance_V2': float(tcentered[0]),
           'thermal_mean_removal_V2': float(tcov[0]-tcentered[0]),
           'total_variance_V2': float(cov[0]), 'window_s': maxlag/fs,
           'statistic': 'within-window sample-mean-removed excess',
           'population_variance': None,
           'uncertainty': 'slow-process mean-removal bias and population uncertainty unknown; '
                          'thermal sensitivity is not total uncertainty',
           'response_assumption': 'clipped training PSD/model transfer' if apply_measured
                                  else 'declared diode only; clipping sensitivity alternative',
           'assumption': 'additive independent component; thermal PSD from other half'}
    if not np.all(np.isfinite(excess)) or excess[0] <= 0:
        return {**out, 'verdict': 'unresolved', 'reason': 'no positive finite excess variance'}
    c = excess/excess[0]
    hits = np.flatnonzero(c <= 1/np.e)
    tau = float(hits[0]/fs) if len(hits) else None
    picks = np.unique(np.r_[0, np.geomspace(1, maxlag, 100).astype(int)])
    return {**out, 'verdict': 'conditional-window-statistic',
            'rms_V': float(np.sqrt(excess[0])),
            'first_e_crossing_s': tau,
            'crossing_lower_bound_s': None if tau is not None else maxlag/fs,
            'tau_over_thermal': None if tau is None else tau*2*np.pi*thermal.fc_fit,
            'min_normalized_covariance': float(c.min()),
            'covariance_sample': {'lag_s': (picks/fs).tolist(),
                                  'normalized_excess': c[picks].tolist()}}


def sensitivity_summary(estimates):
    values = [e['rms_V'] for e in estimates if e.get('rms_V') is not None]
    return {'rms_V_range': [min(values), max(values)] if values else None,
            'n_resolved': len(values), 'n_unresolved': len(estimates)-len(values),
            'verdict': 'complete' if len(values) == len(estimates) else 'partly-unresolved',
            'note': 'training-block model sensitivity only; not total uncertainty',
            'unresolved_reasons': [e.get('reason') for e in estimates if 'rms_V' not in e]}


def contaminated_axis(f, channel):
    x, de, meta = read_axis(f, channel)
    fs = meta['fs_hz']
    halves = np.array_split(x, 2)
    ret, fit_errors = [], []
    for half in halves:
        try:
            ret.append(fit(half, fs, de))
            fit_errors.append(None)
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            ret.append(None)
            fit_errors.append(str(exc))
    b2ref = 2*de['D (V^2/s)']
    theta_ref = 2*np.pi*de['fc (Hz)']
    varref = b2ref/(2*theta_ref)
    out = {'provenance': meta, 'halves': [], 'eight_block_observations': []}
    start = 0
    for block in np.array_split(x, 8):
        out['eight_block_observations'].append({
            'window_s': [start/fs, (start+len(block))/fs],
            'mean_V': float(block.mean()), 'variance_V2': float(block.var())})
        start += len(block)
    for k, held in enumerate(halves):
        trained = ret[1-k]
        vr = float(np.var(held))/(varref*trained.variance) if trained else None
        gate = axis_gate(held, 1/fs, theta_ref, var_ratio=vr, min_lag=8, max_lag=200)
        apparent = gate.get('theta_acf')
        h = {'evaluation_half': k, 'training_half': 1-k,
             'window_s': [k*len(halves[0])/fs, (k*len(halves[0])+len(held))/fs],
             'axis_gate': gate, 'law_certificate': False,
             'attribution': {'verdict': 'unreadable'}}
        out['halves'].append(h)
        if trained is None:
            h['reason'] = 'thermal fit unavailable: '+str(fit_errors[1-k])
            continue
        thermal_theta = 2*np.pi*trained.fc_fit
        h['theta_thermal_per_s'] = thermal_theta
        h['theta_apparent_per_s'] = apparent
        h['theta_thermal_over_calibration'] = thermal_theta/theta_ref
        if apparent is None or not np.isfinite(apparent) or apparent <= 0:
            h['reason'] = 'apparent ACF timescale unavailable'
            continue
        h['theta_apparent_over_thermal'] = apparent/thermal_theta
        strides = [s for s in STRIDES if s*thermal_theta/fs <= .3]
        d = realized_diffusion(held, 1/fs, thermal_theta, trained, strides=strides)
        h['short_lag'] = d
        if d.get('plateau') is None or not np.isfinite(d['plateau']):
            h['reason'] = 'no usable short-lag diffusion estimate'
            continue
        old = realized_diffusion(held, 1/fs, apparent, trained)
        h['old_diffusion_ratio'] = old['plateau']/b2ref if old.get('plateau') is not None else None
        ratios = [apparent/theta_ref, d['plateau']/b2ref, vr]
        at = attribute_deviation(*ratios)
        h.update(ratios=ratios, attribution=at)
        if at['verdict'] == 'slow-contaminant':
            observed = observed_covariance(held, fs)
            h['contaminant'] = covariance_measurement(held, fs, trained, observed=observed)
            h['without_clipped_transfer'] = covariance_measurement(
                held, fs, trained, observed=observed, apply_measured=False)
            h['blocks'] = [covariance_measurement(b, fs, trained)
                           for b in np.array_split(held, 4)]
            estimates = []
            for block in np.array_split(halves[1-k], 4):
                try:
                    model = fit(block, fs, de)
                    estimates.append(covariance_measurement(held, fs, model, observed=observed))
                except (ValueError, RuntimeError, FloatingPointError) as exc:
                    estimates.append({'verdict': 'unresolved', 'reason': str(exc)})
            h['thermal_subtraction_sensitivity'] = sensitivity_summary(estimates)
    th = [h['axis_gate']['ratio'] for h in out['halves'] if h['axis_gate'].get('ratio')]
    spread = max(th)/min(th)-1 if len(th) == 2 else None
    out['theta_stability'] = {'half_ratio_spread': spread,
                             'stable_within_attribution_tol': spread is not None and spread <= .08,
                             'note': 'between-window sensitivity, not a physical drift uncertainty interval'}
    for h in out['halves']:
        h['interpretation'] = ('signature-compatible only; theta instability prevents physical-cause identification'
                               if spread is None or spread > .08 else 'empirical signature comparison only')
    return out


def moment(ret, stride, *, learned=False):
    response = ret.diode * (ret.measured if learned else 1)
    return float(np.trapezoid(ret._weight(stride)*ret.model*response, ret.freqs))


def validate_diode_model(de):
    at_calibration = adapter.bfp_diode_at_power(de, [de['Trap sum power (V)']])
    for key in ('alpha', 'f_diode (Hz)'):
        if not np.isclose(at_calibration[key], de[key], rtol=1e-12, atol=0):
            raise ValueError(f'diode provenance mismatch: {key} does not reproduce stored calibration')


def floor_comparison(fr, ptrain, pheld, thermal):
    """Separate correction training, temporal repeatability and new test bands.

    Transfer is not independently identified: even a pass cannot establish
    detector-only whiteness. Stopband values are diagnostics, never test inputs.
    """
    def mean_in(p, lo, hi):
        mask = (fr >= lo) & (fr < hi)
        if not mask.any():
            raise ValueError(f'no spectral bins in [{lo}, {hi})')
        return float(np.mean(p[mask]))
    N = mean_in(ptrain-thermal, 30000., 40000.)
    def summary(lo, hi):
        observed = mean_in(pheld, lo, hi)
        modeled = mean_in(thermal, lo, hi)
        excess = observed-modeled
        trained_excess = mean_in(ptrain-thermal, lo, hi)
        return {'band_hz': [lo, hi], 'heldout_PSD_V2_per_Hz': observed,
                'thermal_without_antialias_PSD_V2_per_Hz': modeled,
                'excess_PSD_V2_per_Hz': excess,
                'excess_over_train_floor': excess/N if N > 0 else None,
                'train_excess_PSD_V2_per_Hz': trained_excess,
                'half_difference_PSD_V2_per_Hz': abs(excess-trained_excess)/2}
    validation = [summary(20000., 30000.), summary(40000., 43000.)]
    passed = N > 0 and all(abs(b['excess_over_train_floor']-1) <= .10 for b in validation)
    return {'evidence': 'empirical', 'train_floor_PSD_V2_per_Hz': N,
            'training_band_hz': [30000., 40000.],
            'training_band_repeatability': summary(30000., 40000.),
            'passband_validation': validation, 'passband_model_passed': bool(passed),
            'rolloff_diagnostics': [summary(lo, hi) for lo, hi in
                                   ((43000.,45000.), (45000.,48000.), (48000.,50000.))],
            'detector_whiteness': 'unresolved', 'sigma_obs_V': None,
            'note': 'new exploratory passband protocol; unknown transfer and thermal-model error '
                    'prevent detector-only whiteness inference; stopband excluded from decision'}


def noise_axis(f, channel):
    x, de, meta = read_axis(f, channel)
    fs = meta['fs_hz']
    power_data = f[f'Diagnostics/Trap power {channel[6]}'][()]
    power_values = power_data['Value'] if power_data.dtype.names else power_data
    current = adapter.bfp_diode_at_power(de, power_values)
    validate_diode_model(de)
    de = {**de, 'alpha': current['alpha'], 'f_diode (Hz)': current['f_diode (Hz)']}
    train, held = np.array_split(x, 2)
    rt = fit(train, fs, de, (100., 2300.))
    fr, ptrain = welch(train, fs=fs, nperseg=65536)
    _, pheld = welch(held, fs=fs, nperseg=65536)
    thermal = rt.model*rt.diode
    noise = floor_comparison(fr, ptrain, pheld, thermal)
    N = noise['train_floor_PSD_V2_per_Hz']
    moments = []
    for stride in (1, 2, 4, 8, 16):
        observed = float(np.mean((held[stride:]-held[:-stride])**2))
        declared = moment(rt, stride)
        noise_corrected = declared + N*fs  # iid floor contributes 2 sigma²
        moments.append({'stride': stride, 'tau_s': stride/fs,
                        'measured_V2': observed,
                        'declared_thermal_V2': declared,
                        'declared_plus_floor_V2': noise_corrected,
                        'observed_over_thermal': observed/declared,
                        'observed_over_thermal_plus_floor': observed/noise_corrected
                        if noise_corrected > 0 else None,
                        'learned_response_thermal_V2': moment(rt, stride, learned=True)})
    dec = held[::2]
    decimation = []
    # Equal-lag decimated increments are a subsample of native increments.
    # Native-band PSD prediction measures temporal repeatability, not unknown
    # above-native-Nyquist power. Truncation is an arithmetic counterfactual.
    for stride in (1, 2, 4, 8):
        native = float(np.mean((held[2*stride:]-held[:-2*stride])**2))
        observed = float(np.mean((dec[stride:]-dec[:-stride])**2))
        mask = fr <= fs/4
        weight = 4*np.sin(np.pi*fr[mask]*(2*stride)/fs)**2
        truncated = float(np.trapezoid(weight*(thermal[mask]+N), fr[mask]))
        full_weight = 4*np.sin(np.pi*fr*(2*stride)/fs)**2
        native_band = float(np.trapezoid(full_weight*ptrain, fr))
        truncated_band = float(np.trapezoid(full_weight[mask]*ptrain[mask], fr[mask]))
        decimation.append({'decimated_stride': stride,
                           'observed_over_training_native_band_PSD': observed/native_band,
                           'observed_over_training_truncated_PSD': observed/truncated_band,
                           'counterfactual_discarded_fraction': 1-truncated_band/native_band, 'tau_s': 2*stride/fs,
                           'subsample_consistency_ratio': observed/native,
                           'observed_over_truncated_thermal_white_counterfactual': observed/truncated
                           if truncated > 0 else None})
    idx = np.unique(np.geomspace(1, len(fr)-1, 180).astype(int))
    return {'provenance': meta, 'diode_at_record_power': current, 'fit_half': 0, 'evaluation_half': 1,
            'thermal_fit': {'range_hz': [100,2300], 'fc_hz': rt.fc_fit,
                            'amplitude': rt.amplitude, 'alpha': rt.alpha,
                            'f_diode_hz': rt.f_diode},
            'noise': noise,
            'decimation_interpretation': 'subsample consistency and native-band temporal spectral repeatability only',
            'retention_factors': rt.fraction(1, factors=True),
            'heldout_moments': moments,
            'moment_model_note': 'thermal+white-floor counterfactual omits unknown anti-alias transfer',
            'decimation': decimation,
            'spectrum_sample': {'frequency_hz': fr[idx].tolist(),
                                'train_PSD': ptrain[idx].tolist(),
                                'heldout_PSD': pheld[idx].tolist(),
                                'frozen_thermal_PSD': thermal[idx].tolist()}}


def provenance(path):
    return {'file': path.name, 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    import h5py
    path = adapter.BFP_DIR/'passive_calibration.h5'
    out = {'stage': 'C4', 'tag': 'empirical',
           'registration': 'docs/TWEEZERS_C4_REGISTRATION.md',
           'inputs': [provenance(path)],
           'old_triple': attribute_deviation(.019, .885, 3.797),
           'R-C3': {'status': 'open', 'stage_length_unit': None,
                    'note': 'no new stage unit or drive frequency is supplied'}}
    with h5py.File(path) as f:
        out['passive'] = {ch: contaminated_axis(f, ch) for ch in ('Force 1x', 'Force 1y')}
    path = adapter.BFP_DIR/'noise_floor.h5'
    out['inputs'].append(provenance(path))
    with h5py.File(path) as f:
        out['noise_floor'] = {}
        for ch in f['Force HF']:
            try:
                out['noise_floor'][ch] = noise_axis(f, ch)
            except (ValueError, RuntimeError, KeyError, FloatingPointError) as exc:
                out['noise_floor'][ch] = {'verdict': 'unresolved', 'reason': str(exc)}
        out['noise_diagnostics'] = {}
        for key, d in f['Diagnostics'].items():
            values = d[()]
            if values.dtype.names:
                values = values['Value']
            out['noise_diagnostics'][key] = {'mean': float(np.mean(values)),
                                            'min': float(np.min(values)),
                                            'max': float(np.max(values))}
    c3path = OUT.parent/'tweezers_c3.json'
    c3 = json.loads(c3path.read_text())
    replay = {
        'evidence': 'empirical', 'source': provenance(c3path),
        'method': 'current classifier on frozen C3 ratios; no data or stage rerun',
        'axes': {record+'/'+ch: {'historical': a['attribution'],
                               'current': attribute_deviation(a['ratios']['theta'],
                                  a['ratios']['diffusion'], a['ratios']['variance'])}
                 for record, block in c3['records'].items() for ch, a in block['axes'].items()}}
    replay_path = OUT.parent/'tweezers_c3_attribution_replay.json'
    replay_path.write_text(json.dumps(replay, indent=2, allow_nan=False)+'\n')
    out['c3_classifier_replay'] = provenance(replay_path)
    out['c3_ratio_controls'] = {
        record+'/'+ch: {'ratios': a['ratios'],
                       'verdict': attribute_deviation(a['ratios']['theta'],
                           a['ratios']['diffusion'], a['ratios']['variance'])['verdict']}
        for record, block in c3['records'].items() for ch, a in block['axes'].items()}
    out['comparison_78125Hz'] = {
        'source_artifact': provenance(c3path),
        'axes': {record+'/'+ch: a['retention']
                 for record, block in c3['records'].items()
                 for ch, a in block['axes'].items()},
        'verdict': 'cross-record rate attribution unresolved',
        'note': 'different bead size, power, diode and acquisition; not a controlled fs-only experiment'}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
    for ch, a in out['noise_floor'].items():
        if 'noise' not in a:
            print('NOISE', ch, a)
            continue
        print('NOISE', ch, a['thermal_fit'], a['noise'])
        print('MOMENTS', [(m['stride'], m['observed_over_thermal'],
                           m['observed_over_thermal_plus_floor']) for m in a['heldout_moments']])
        print('DECIMATION', a['decimation'])
    for ch, a in out['passive'].items():
        print(ch, [(h.get('ratios'), h['attribution']['verdict']) for h in a['halves']])
        for h in a['halves']:
            if 'contaminant' in h:
                print({k:v for k,v in h['contaminant'].items() if k!='covariance_sample'})


if __name__ == '__main__':
    main()
