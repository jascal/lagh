"""Standalone C4 research figure, generated from the committed numeric artifact."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / 'results'
r = json.loads((ROOT/'tweezers_c4.json').read_text())
fig, ax = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
b = r['passive']['Force 1x']['eight_block_observations']
t = [np.mean(v['window_s']) for v in b]
means = np.array([v['mean_V'] for v in b])
ax[0, 0].plot(t, 1000*(means-means[0]), 'o-')
ax[0, 0].axvspan(8.0169216, 16.0338432, alpha=.1, color='orange')
ax[0, 0].set(title='Passive 1x: late level shift', xlabel='Record time (s)',
             ylabel='Block mean relative to first block (mV)')
h = r['passive']['Force 1x']['halves'][1]['contaminant']
c = h['covariance_sample']
ax[0, 1].plot(c['lag_s'], c['normalized_excess'])
ax[0, 1].axhline(1/np.e, color='gray', linestyle='--', label='1/e')
ax[0, 1].set(title='Second-half excess covariance: no crossing by 0.5 s',
             xlabel='Lag (s)', ylabel='Normalized excess covariance', ylim=(0,1.1))
ax[0, 1].legend()
for ch, a in r['noise_floor'].items():
    ax[1, 0].plot([25, 35, 44], [v['excess_over_train_floor'] for v in a['noise']['bands']],
                  'o-', label=ch)
    ds = a['decimation']
    ax[1, 1].plot([d['tau_s']*1e6 for d in ds],
                  [d['observed_over_train_PSD_alias_retaining'] for d in ds], 'o-', label=ch)
    ax[1, 1].plot([d['tau_s']*1e6 for d in ds],
                  [d['observed_over_train_PSD_alias_discarding'] for d in ds], '--', color=ax[1,1].lines[-1].get_color())
ax[1, 0].axhspan(.9, 1.1, alpha=.1, color='gray')
ax[1, 0].set(title='Held-out spectral residual is not a white floor',
             xlabel='Band center (kHz)', ylabel='Excess / training 30–40 kHz excess')
ax[1, 0].legend(fontsize=8)
ax[1, 1].axhline(1, color='gray', linewidth=.7)
ax[1, 1].set(title='Decimation: retain aliases (solid), discard (dashed)',
             xlabel='Physical lag (µs)', ylabel='Observed / predicted increment variance')
fig.suptitle('Tweezers C4 — empirical, window- and model-conditional; R-C3 remains open')
fig.savefig(ROOT/'tweezers_c4.png', dpi=160)
