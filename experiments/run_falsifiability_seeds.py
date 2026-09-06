"""Four registered seed interventions; synthetic evidence only."""
import json
from pathlib import Path
import sys

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tests'))
from test_instrument import bead, FS, RD_UM_PER_V
from lagh.instrument import drive_scale


def main():
    stage,volts,truth=bead(drag_over_bulk=2.3,n=400000,seed=4)
    given=drive_scale(stage,volts,FS,f_c=truth['fc_bulk_referenced'])
    read=drive_scale(stage,volts,FS)
    n=np.arange(10000)
    x=np.sin(2*np.pi*.4*n)+.1*np.sin(2*np.pi*.023*n)
    paired=x[2:]-x[:-2]
    identity_error=float(np.max(np.abs(np.diff(x[::2])-paired[::2])))
    frequency=np.arange(0.,50001.,10.)
    white=np.ones(len(frequency))
    transfer=np.where(frequency>43000.,.01,1.)
    observed=white*transfer
    band=(frequency>=40000)&(frequency<48000)
    old_whiteness_error=float(abs(observed[band].mean()-1))
    # Training-band repeatability fails, disjoint validation bands remain flat.
    training_repeatability=[2.,1.,1.]
    disjoint=[1.,1.]
    bumped_disjoint=[2.,1.]
    out=dict(evidence='empirical',
      scope='computed synthetic responses to the declared inputs below',
      declarations=dict(drag_ratio=2.3, seed=4, n=400000,
        detector_psd='constant 1', transfer='H²=.01 above43kHz, 1 elsewhere',
        training_repeatability=training_repeatability, disjoint_validation=disjoint,
        perturbed_validation=bumped_disjoint,
        note='Constructed inputs, not empirical findings; flags below are computed predicates.'),
      circular_corner=dict(input='bead drag=2.3, seed=4, n=400000, declared bulk corner',
        true_scale=RD_UM_PER_V,assumed_corner_scale=given['scale'],record_corner_scale=read['scale'],
        assumed_relative_error=abs(given['scale']/RD_UM_PER_V-1),
        note='Changing true drag defeats independence while supplied-corner inversion remains evaluable.'),
      decimation=dict(input='sin(2*pi*.4*n)+.1*sin(2*pi*.023*n), n=0..9999',
        paired_identity_error=identity_error,scope="paired-sample identity; no independent response validation"),
      whiteness=dict(input='white PSD=1; H²=.01 above43kHz; subtracted thermal PSD=0',
        old_band_relative_error=old_whiteness_error,declared_detector_psd='constant one-sided PSD of 1',
        construction_identity_equal=bool(np.array_equal(observed,transfer)),
        scope="identical constructed spectra do not independently test whiteness"),
      training_band_vote=dict(input='training-repeat ratio2; disjoint-band ratios1,1',
        old_all_bands_pass=all(abs(r-1)<=.1 for r in training_repeatability),
        disjoint_bands_pass=all(abs(r-1)<=.1 for r in disjoint),
        independent_failure_input='change first disjoint-band ratio to2',
        independent_failure_pass=all(abs(r-1)<=.1 for r in bumped_disjoint)))
    Path('experiments/results/falsifiability_seeds.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
