"""C4 figure from the numeric artifact; unresolved results remain visible."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / 'results'


def plot(r, output):
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    passive = r.get('passive', {}).get('Force 1x', {})
    b = passive.get('eight_block_observations', [])
    if b:
        means = np.array([v['mean_V'] for v in b])
        ax[0, 0].plot([np.mean(v['window_s']) for v in b], 1000*(means-means[0]), 'o-')
    ax[0, 0].set(title='Passive 1x: within-window level changes', xlabel='Record time (s)',
                 ylabel='Block mean relative to first block (mV)')
    halves = passive.get('halves', [])
    h = halves[1].get('contaminant', {}) if len(halves) > 1 else {}
    c = h.get('covariance_sample')
    if c:
        ax[0, 1].plot(c['lag_s'], c['normalized_excess'])
        ax[0, 1].axhline(1/np.e, color='gray', linestyle='--', label='1/e')
        ax[0, 1].legend()
    else:
        ax[0, 1].text(.5, .5, 'Excess covariance unresolved', ha='center', transform=ax[0, 1].transAxes)
    ax[0, 1].set(title='Demeaned excess: window statistic only',
                 xlabel='Lag (s)', ylabel='Normalized within-window excess')
    for ch, a in r.get('noise_floor', {}).items():
        if 'noise' not in a:
            continue
        bands = a['noise']['passband_validation']
        ax[1, 0].plot([np.mean(v['band_hz'])/1000 for v in bands],
                      [v['excess_over_train_floor'] for v in bands], 'o-', label=ch)
        ds = a['decimation']
        line, = ax[1, 1].plot([d['tau_s']*1e6 for d in ds],
                              [d['observed_over_training_native_band_PSD'] for d in ds], 'o-', label=ch)
        ax[1, 1].plot([d['tau_s']*1e6 for d in ds],
                      [d['observed_over_training_truncated_PSD'] for d in ds], '--', color=line.get_color())
    ax[1, 0].axhspan(.9, 1.1, alpha=.1, color='gray')
    ax[1, 0].set(title='Passband model check; whiteness unresolved',
                 xlabel='Validation band center (kHz)', ylabel='Excess / training 30–40 kHz excess')
    if ax[1, 0].lines:
        ax[1, 0].legend(fontsize=8)
    ax[1, 1].axhline(1, color='gray', linewidth=.7)
    ax[1, 1].set(title='Native-band repeatability\nDashed: counterfactual truncation',
                 xlabel='Physical lag (µs)', ylabel='Observed / predicted increment variance')
    fig.suptitle('Tweezers C4 — empirical window statistics; detector whiteness and R-C3 open')
    fig.savefig(output, dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    plot(json.loads((ROOT/'tweezers_c4.json').read_text()), ROOT/'tweezers_c4.png')
