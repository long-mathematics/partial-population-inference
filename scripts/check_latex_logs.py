"""Reject compilation failures, unresolved references, and layout warnings."""
from pathlib import Path
import re
from pdfs import ROOT, SOURCES

problems=[]
for source in SOURCES:
    log=ROOT/'build'/(Path(source).stem+'.log')
    text=log.read_text(errors='replace')
    if 'Output written on' not in text:
        problems.append(f'{log.name}: no completed PDF')
    for pattern in (r'^!.*',r'^.*Warning:.*',r'^Overfull .*',r'^Underfull .*'):
        problems.extend(f'{log.name}: {line}' for line in re.findall(pattern,text,re.MULTILINE))
if problems:
    raise SystemExit('\n'.join(problems))
print(f'{len(SOURCES)} LaTeX logs: no errors, warnings, or overfull/underfull boxes.')
