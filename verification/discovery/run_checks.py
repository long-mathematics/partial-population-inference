#!/usr/bin/env python3
"""Run the retained and independent verification programs; write separate logs."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
TASKS = [
    ('sampling', ROOT/'originals/sampling_calculus/check_calculus.py'),
    ('recovery', ROOT/'originals/recovery_geometry/recovery_geometry.py'),
    ('calibration', ROOT/'originals/statistical_calibration/statistical_calibration.py'),
    ('independent', ROOT/'checks/independent_audit.py'),
    ('gt_variance', ROOT/'checks/gt_variance_check.py'),
    ('revision_v02', ROOT/'checks/revision_v02.py'),
]

def main() -> int:
    target = ROOT/'checks'
    target.mkdir(exist_ok=True)
    for name, script in TASKS:
        if not script.is_file():
            raise FileNotFoundError(script)
        completed = subprocess.run([sys.executable, str(script)], cwd=ROOT,
                                   capture_output=True, text=True, check=False)
        (target/f'{name}_rerun.log').write_text(completed.stdout+completed.stderr)
        if completed.returncode:
            print(f'FAILED: {name}; see checks/{name}_rerun.log', file=sys.stderr)
            return completed.returncode
        print(f'PASSED: {name}')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
