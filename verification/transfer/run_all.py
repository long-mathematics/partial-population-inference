#!/usr/bin/env python3
"""Run both verification suites; fail immediately on any failed assertion."""
from pathlib import Path
import json
import platform
import subprocess
import sys
import numpy
import scipy
import sympy

ROOT=Path(__file__).resolve().parent
for name in ['check_finite_original.py','check_extensions.py']:
    print(f'Running {name}',flush=True)
    result=subprocess.run([sys.executable,str(ROOT/name)],cwd=ROOT,
                          text=True,capture_output=True)
    (ROOT/(name.removesuffix('.py')+'.log')).write_text(result.stdout+result.stderr,
                                                       encoding='utf-8')
    if result.returncode:
        print(result.stdout); print(result.stderr,file=sys.stderr)
        raise SystemExit(result.returncode)
original=json.loads((ROOT/'verification.json').read_text())
extensions=json.loads((ROOT/'extensions_results.json').read_text())
summary={
    'all_passed':original['all_passed'] and extensions['all_passed'],
    'retained_exact_checks':original['total_checks'],
    'extension_checks':extensions['total_checks'],
    'total_checks':original['total_checks']+extensions['total_checks'],
    'environment':{'python':platform.python_version(),'numpy':numpy.__version__,
                   'scipy':scipy.__version__,'sympy':sympy.__version__},
    'scope':'Exact finite algebra and explicit numerical corroboration, not formal verification or a novelty certificate.'
}
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
