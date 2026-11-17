"""Publish built PDFs, or check that committed downloads match their sources."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
import shutil

ROOT=Path(__file__).resolve().parents[1]
SOURCES=(
    'manuscripts/discovery_rate_inversion.tex',
    'manuscripts/transfer_generation.tex',
    'research/master_identity_and_deletion.tex',
    'manuscripts/capture_recapture.tex',
    'manuscripts/classical_data.tex',
)

def dependencies(source):
    if source!='manuscripts/classical_data.tex':return []
    return sorted(str(p.relative_to(ROOT)) for p in (ROOT/'manuscripts/generated').glob('*'))+[
        'verification/empirical/classical_results.json','verification/empirical/coverage_width.json',
        'data/processed/wartime_comparison.json','data/processed/frequency_variants.json']

def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args()
    manifest=ROOT/'pdf/source_manifest.json'
    if args.publish:
        (ROOT/'pdf').mkdir(exist_ok=True)
        old={e['source']:e for e in json.loads(manifest.read_text())['papers']} if manifest.exists() else {}
        entries=[]
        for source in SOURCES:
            filename=Path(source).stem+'.pdf'
            built=ROOT/'build'/filename
            data=built.read_bytes()
            if not data.startswith(b'%PDF-') or b'%%EOF' not in data[-1024:]:
                raise ValueError(f'Not a complete PDF: {built}')
            target='pdf/'+filename
            deps={p:digest(p) for p in dependencies(source)}
            prior=old.get(source,{})
            unchanged=(prior.get('source_sha256')==digest(source) and prior.get('dependencies',{})==deps
                       and (ROOT/target).exists() and prior.get('pdf_sha256')==digest(target))
            if not unchanged:shutil.copyfile(built,ROOT/target)
            entries.append({'source':source,'source_sha256':digest(source),
                            'pdf':target,'pdf_sha256':digest(target),
                            'dependencies':deps})
        manifest.write_text(json.dumps({'papers':entries},indent=2)+'\n')
    entries=json.loads(manifest.read_text())['papers']
    if [entry['source'] for entry in entries]!=list(SOURCES):
        raise ValueError('PDF source inventory changed; run make publish-pdfs')
    for entry in entries:
        if sorted(entry.get('dependencies',{}))!=sorted(dependencies(entry['source'])):
            raise ValueError('PDF dependency inventory changed; run make publish-pdfs')
        for path,expected in entry.get('dependencies',{}).items():
            if digest(path)!=expected:raise ValueError(f'Stale PDF input: {path}; run make publish-pdfs')
        for kind in ('source','pdf'):
            if digest(entry[kind])!=entry[kind+'_sha256']:
                raise ValueError(f"Stale {kind}: {entry[kind]}; run make publish-pdfs")
    print(f'{len(entries)} downloadable PDFs match their recorded source and PDF hashes.')

if __name__=='__main__':
    main()
