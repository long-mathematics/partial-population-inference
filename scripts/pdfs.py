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
)

def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--publish',action='store_true')
    args=parser.parse_args()
    manifest=ROOT/'pdf/source_manifest.json'
    if args.publish:
        (ROOT/'pdf').mkdir(exist_ok=True)
        entries=[]
        for source in SOURCES:
            filename=Path(source).stem+'.pdf'
            built=ROOT/'build'/filename
            data=built.read_bytes()
            if not data.startswith(b'%PDF-') or b'%%EOF' not in data[-1024:]:
                raise ValueError(f'Not a complete PDF: {built}')
            target='pdf/'+filename
            shutil.copyfile(built,ROOT/target)
            entries.append({'source':source,'source_sha256':digest(source),
                            'pdf':target,'pdf_sha256':digest(target)})
        manifest.write_text(json.dumps({'papers':entries},indent=2)+'\n')
    entries=json.loads(manifest.read_text())['papers']
    if [entry['source'] for entry in entries]!=list(SOURCES):
        raise ValueError('PDF source inventory changed; run make publish-pdfs')
    for entry in entries:
        for kind in ('source','pdf'):
            if digest(entry[kind])!=entry[kind+'_sha256']:
                raise ValueError(f"Stale {kind}: {entry[kind]}; run make publish-pdfs")
    print(f'{len(entries)} downloadable PDFs match their recorded source and PDF hashes.')

if __name__=='__main__':
    main()
