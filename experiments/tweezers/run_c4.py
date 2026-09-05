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
from lagh.instrument import (acf_timescale, attribute_deviation,  # noqa: E402
                             realized_diffusion, retention)

OUT = Path(__file__).resolve().parents[1] / 'results/tweezers_c4.json'
STRIDES = (2, 3, 4, 6, 8)


def read_axis(f, channel):
    _, x, applied, de, meta = adapter.bfp_position_nm(f, channel)
    return x / (de['Rd (um/V)'] * 1e3), de, {**meta,
        'applied': applied, 'reference_calibration': de,
        'record_start_ns': int(f[f'Force HF/{channel}'].attrs['Start time (ns)'])}


def fit(x, fs, de, fit_range=None):
    return retention(x, fs, alpha=de['alpha'], f_diode=de['f_diode (Hz)'],
                     fit_range=fit_range or (de['Fit range (min.) (Hz)'],
                                            de['Fit range (max.) (Hz)']))


def covariance_measurement(x, fs, thermal):
    """Recorded excess covariance conditional on the independent thermal PSD.

    No parametric contaminant model: report first 1/e crossing, not an OU rate.
    FFT integration of the retained thermal PSD uses the same sampling grid as
    the observed covariance. This is an empirical, finite-record subtraction.
    """
    x = np.asarray(x) - np.mean(x)
    n = len(x)
    nfft = next_fast_len(2*n)
    maxlag = min(int(.5*fs), n//4)
    power = np.abs(rfft(x, nfft))**2
    cov = irfft(power, nfft)[:maxlag+1] / (n-np.arange(maxlag+1))
    fr = np.fft.rfftfreq(nfft, 1/fs)
    psd = np.interp(fr, thermal.freqs, thermal.model*thermal.diode*thermal.measured)
    tcov = irfft(psd, nfft)[:maxlag+1] * fs/2
    excess = cov-tcov
    out = {'evidence': 'empirical', 'recorded_excess_variance_V2': float(excess[0]),
           'thermal_retained_variance_V2': float(tcov[0]),
           'total_variance_V2': float(cov[0]), 'window_s': maxlag/fs,
           'assumption': 'additive independent contaminant; thermal PSD from other half'}
    if excess[0] <= 0:
        return {**out, 'verdict': 'unresolved', 'reason': 'no positive excess variance'}
    c = excess/excess[0]
    hits = np.flatnonzero(c <= 1/np.e)
    tau = float(hits[0]/fs) if len(hits) else None
    picks = np.unique(np.r_[0, np.geomspace(1, maxlag, 100).astype(int)])
    return {**out, 'verdict': 'conditional-measurement',
            'rms_V': float(np.sqrt(excess[0])),
            'first_e_crossing_s': tau,
            'crossing_lower_bound_s': None if tau is not None else maxlag/fs,
            'tau_over_thermal': None if tau is None else tau*2*np.pi*thermal.fc_fit,
            'min_normalized_covariance': float(c.min()),
            'covariance_sample': {'lag_s': (picks/fs).tolist(),
                                  'normalized_excess': c[picks].tolist()}}


def contaminated_axis(f, channel):
    x, de, meta = read_axis(f, channel)
    fs = meta['fs_hz']
    halves = np.array_split(x, 2)
    ret = [fit(h, fs, de) for h in halves]
    b2ref = 2*de['D (V^2/s)']
    varref = b2ref/(4*np.pi*de['fc (Hz)'])
    out = {'provenance': meta, 'halves': [],
           'eight_block_observations': [
               {'window_s': [i*len(x)/8/fs, (i+1)*len(x)/8/fs],
                'mean_V': float(b.mean()), 'variance_V2': float(b.var())}
               for i, b in enumerate(np.array_split(x, 8))]}
    for k, held in enumerate(halves):
        trained = ret[1-k]
        apparent = acf_timescale(held, 1/fs, min_lag=8, max_lag=200)['theta']
        thermal_theta = 2*np.pi*trained.fc_fit
        strides = [s for s in STRIDES if s*thermal_theta/fs <= .3]
        d = realized_diffusion(held, 1/fs, thermal_theta, trained, strides=strides)
        old = realized_diffusion(held, 1/fs, apparent, trained)
        vr = float(np.var(held))/(varref*trained.variance)
        ratios = [apparent/(2*np.pi*de['fc (Hz)']), d['plateau']/b2ref, vr]
        at = attribute_deviation(*ratios)
        h = {'evaluation_half': k, 'window_s': [k*len(halves[0])/fs, (k*len(halves[0])+len(held))/fs], 'training_half': 1-k,
             'theta_thermal_per_s': thermal_theta, 'theta_apparent_per_s': apparent,
             'ratios': ratios, 'short_lag': d,
             'old_diffusion_ratio': old['plateau']/b2ref,
             'attribution': at}
        if at['verdict'] == 'slow-contaminant':
            h['contaminant'] = covariance_measurement(held, fs, trained)
            h['blocks'] = [covariance_measurement(b, fs, trained)
                           for b in np.array_split(held, 4)]
            alternative_models = [fit(b, fs, de) for b in np.array_split(halves[1-k], 4)]
            estimates = [covariance_measurement(held, fs, model)
                         for model in alternative_models]
            h['thermal_subtraction_sensitivity'] = {
                'rms_V_range': [min(e['rms_V'] for e in estimates),
                                max(e['rms_V'] for e in estimates)],
                'note': 'four training-block models on same evaluation half; not a confidence interval'}
        out['halves'].append(h)
    return out


def moment(ret, stride, *, learned=False):
    response = ret.diode * (ret.measured if learned else 1)
    return float(np.trapezoid(ret._weight(stride)*ret.model*response, ret.freqs))


def noise_axis(f, channel):
    x, de, meta = read_axis(f, channel)
    fs = meta['fs_hz']
    power_data = f[f'Diagnostics/Trap power {channel[6]}'][()]
    power_values = power_data['Value'] if power_data.dtype.names else power_data
    current = adapter.bfp_diode_at_power(de, power_values)
    at_calibration = adapter.bfp_diode_at_power(de, [de['Trap sum power (V)']])
    assert np.isclose(at_calibration['alpha'], de['alpha'], rtol=1e-12)
    assert np.isclose(at_calibration['f_diode (Hz)'], de['f_diode (Hz)'], rtol=1e-12)
    de = {**de, 'alpha': current['alpha'], 'f_diode (Hz)': current['f_diode (Hz)']}
    train, held = np.array_split(x, 2)
    rt = fit(train, fs, de, (100., 2300.))
    fr, ptrain = welch(train, fs=fs, nperseg=65536)
    _, pheld = welch(held, fs=fs, nperseg=65536)
    thermal = rt.model*rt.diode
    def band_mean(p, lo, hi):
        mask = (fr >= lo) & (fr < hi)
        return float(np.mean(p[mask]))
    N = band_mean(ptrain-thermal, 30000., 40000.)
    bands = []
    for lo, hi in ((20000., 30000.), (30000., 40000.), (40000., 48000.)):
        excess = band_mean(pheld-thermal, lo, hi)
        trained_excess = band_mean(ptrain-thermal, lo, hi)
        bands.append({'band_hz': [lo, hi], 'excess_PSD_V2_per_Hz': excess,
                      'excess_over_train_floor': excess/N if N > 0 else None,
                      'excess_band_power_V2': excess*(hi-lo),
                      'train_excess_PSD_V2_per_Hz': trained_excess,
                      'half_difference_PSD_V2_per_Hz': abs(excess-trained_excess)/2})
    white = N > 0 and all(abs(b['excess_over_train_floor']-1) <= .10 for b in bands)
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
    # The alias-free Retention counterfactual truncates the SAME frozen PSD at
    # the new Nyquist. Actual unfiltered decimation retains aliased power.
    for stride in (1, 2, 4, 8):
        native = float(np.mean((held[2*stride:]-held[:-2*stride])**2))
        observed = float(np.mean((dec[stride:]-dec[:-stride])**2))
        mask = fr <= fs/4
        weight = 4*np.sin(np.pi*fr[mask]*(2*stride)/fs)**2
        truncated = float(np.trapezoid(weight*(thermal[mask]+N), fr[mask]))
        full_weight = 4*np.sin(np.pi*fr*(2*stride)/fs)**2
        alias_retaining = float(np.trapezoid(full_weight*ptrain, fr))
        alias_discarding = float(np.trapezoid(full_weight[mask]*ptrain[mask], fr[mask]))
        decimation.append({'decimated_stride': stride,
                           'observed_over_train_PSD_alias_retaining': observed/alias_retaining,
                           'observed_over_train_PSD_alias_discarding': observed/alias_discarding,
                           'lost_tail_fraction': 1-alias_discarding/alias_retaining, 'tau_s': 2*stride/fs,
                           'decimated_over_native_same_lag': observed/native,
                           'decimated_over_unaliased_prediction': observed/truncated
                           if truncated > 0 else None})
    idx = np.unique(np.geomspace(1, len(fr)-1, 180).astype(int))
    return {'provenance': meta, 'diode_at_record_power': current, 'fit_half': 0, 'evaluation_half': 1,
            'thermal_fit': {'range_hz': [100,2300], 'fc_hz': rt.fc_fit,
                            'amplitude': rt.amplitude, 'alpha': rt.alpha,
                            'f_diode_hz': rt.f_diode},
            'noise': {'tag': 'empirical', 'train_floor_PSD_V2_per_Hz': N,
                      'white_floor_passed': bool(white), 'bands': bands,
                      'signed_band_20_48kHz_residual_power_V2': sum(b['excess_band_power_V2'] for b in bands),
                      'sigma_obs_V': float(np.sqrt(N*fs/2)) if white else None,
                      'conditional_white_equivalent_sigma_V': float(np.sqrt(N*fs/2))
                      if N > 0 else None,
                      'note': 'signed thermal-model residuals, not identified detector-only power; '
                              'iid sigma refused unless held-out bands agree within 10%'},
            'retention_factors': rt.fraction(1, factors=True),
            'heldout_moments': moments, 'decimation': decimation,
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
        out['noise_floor'] = {ch: noise_axis(f, ch) for ch in f['Force HF']}
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
        print('NOISE', ch, a['thermal_fit'], a['noise'])
        print('MOMENTS', [(m['stride'], m['observed_over_thermal'],
                           m['observed_over_thermal_plus_floor']) for m in a['heldout_moments']])
        print('DECIMATION', a['decimation'])
    for ch, a in out['passive'].items():
        print(ch, [(h['ratios'], h['attribution']['verdict']) for h in a['halves']])
        for h in a['halves']:
            if 'contaminant' in h:
                print({k:v for k,v in h['contaminant'].items() if k!='covariance_sample'})


if __name__ == '__main__':
    main()
