"""Pilot-only test choice. This proposes weights; it does not certify anything."""
from dataclasses import dataclass
import hashlib

import numpy as np


@dataclass(frozen=True)
class PilotTest:
    function: str
    pilot_sha256: str
    scores: tuple
    drift: tuple
    diffusion: tuple
    evidence: str = 'empirical'


def choose_test(t, pilot_paths) -> PilotTest:
    """Freeze a restricted signal/noise proxy before certification paths exist.

    Drift is fitted on (1,x); b² on (1,x²), nonnegatively. The only candidates
    are x, x²/2 and (on positive pilot paths) log(x). Certification must use
    independent paths and its original bands; a selected function with undefined
    values there must refuse, never trigger reselection. Selection makes no
    truth or coverage claim; the recorded hash audits the pilot input.
    """
    from scipy.optimize import nnls
    t = np.asarray(t, dtype=float)
    p = np.asarray(pilot_paths, dtype=float)
    if (t.ndim != 1 or p.ndim != 2 or p.shape[1] != len(t) or len(t) < 4
            or not np.all(np.isfinite(t)) or not np.all(np.isfinite(p))):
        raise ValueError('finite pilot trajectories aligned with t are required')
    dt = np.diff(t)
    if np.any(dt <= 0) or not np.allclose(dt, dt[0], rtol=1e-8, atol=0):
        raise ValueError('pilot grid must be uniform and increasing')
    x, dx = p[:,:-1].ravel(), np.diff(p,axis=1).ravel()
    if np.var(x) <= 0 or np.mean(dx**2) <= 0:
        raise ValueError('pilot contains no drift/diffusion evidence')
    A = np.column_stack([np.ones(len(x)), x])
    drift = np.linalg.lstsq(A, dx/dt[0], rcond=None)[0]
    B = np.column_stack([np.ones(len(x)), x*x])
    diffusion = nnls(B, dx*dx/dt[0])[0]
    a, b2 = A@drift, B@diffusion
    weights = [('x',np.ones(len(x))), ('x**2/2',x)]
    if np.all(p > 0):
        weights.append(('log(x)',1/x))
    scores = tuple((name,float(abs(np.mean(w*a))/np.sqrt(np.mean(w*w*b2))))
                   for name,w in weights)
    if not all(np.isfinite(value) for _,value in scores):
        raise ValueError('nonfinite pilot score')
    digest = hashlib.sha256(t.tobytes()+p.tobytes()).hexdigest()
    return PilotTest(max(scores,key=lambda item:item[1])[0], digest, scores,
                     tuple(drift.tolist()), tuple(diffusion.tolist()))
