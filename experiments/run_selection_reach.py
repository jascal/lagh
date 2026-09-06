"""Unchanged deterministic reach consumer; pilot state tests do not apply."""
import json
from pathlib import Path
from experiments.reach import audit


def main():
    audit.OUT=Path('experiments/results/test_selection_reach.json')
    audit.main()
    data=json.loads(audit.OUT.read_text())
    # Runtime is telemetry, not a reproducible scientific observable.
    for value in data.values():
        value.pop('seconds',None)
    audit.OUT.write_text(json.dumps(dict(evidence='empirical',
        comparison='Unchanged deterministic consumer; no state-test selection applies.',
        cells=data),indent=2)+'\n')


if __name__=='__main__':
    main()
