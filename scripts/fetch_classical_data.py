"""Fetch pinned CRAN mirror data and documentation; retain hashes and licenses."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,urllib.request
ROOT=Path(__file__).resolve().parents[1]
PACKAGES={'VGAM':['corbet'],'sads':['moths','birds'],'vcdExtra':['ShakeWords'],
          'vegan':['BCI'],'Rcapture':['hare'],'RMark':['edwards.eberhardt'],
          'Petersen':['data_rodli'],'fishmethods':['Gerking']}

def fetch_package(item):
 pkg,names=item
 meta=json.load(urllib.request.urlopen(f'https://api.github.com/repos/cran/{pkg}/commits/master'))
 rev=meta['sha']
 tree=json.load(urllib.request.urlopen(f'https://api.github.com/repos/cran/{pkg}/git/trees/{rev}?recursive=1'))
 entries=[]
 for item in tree['tree']:
  path=item['path']; stem=Path(path).stem
  if path!='DESCRIPTION' and not (stem in names and (path.startswith('data/') or path.startswith('man/'))):continue
  url=f'https://raw.githubusercontent.com/cran/{pkg}/{rev}/{path}'
  data=urllib.request.urlopen(url).read()
  local=Path('data')/('raw' if path.startswith('data/') else 'documentation')/pkg/Path(path).name
  (ROOT/local).parent.mkdir(parents=True,exist_ok=True);(ROOT/local).write_bytes(data)
  entries.append({'package':pkg,'revision':rev,'upstream_path':path,'url':url,'path':str(local),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
 return entries

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--refresh',action='store_true',help='Resolve new upstream revisions and replace the recorded manifest')
 args=parser.parse_args()
 manifest=ROOT/'data/source_manifest.json'
 if manifest.exists() and not args.refresh:
  sources=json.loads(manifest.read_text())['sources']
  for entry in sources:
   data=urllib.request.urlopen(entry['url']).read()
   if hashlib.sha256(data).hexdigest()!=entry['sha256']:
    raise ValueError(f"Upstream content changed: {entry['url']}")
   dest=ROOT/entry['path'];dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
  print(f'Restored {len(sources)} recorded source files with matching SHA-256 hashes.')
  return
 out=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
  for entries in pool.map(fetch_package,PACKAGES.items()):out.extend(entries)
 url='https://www.capturerecapture.co.uk/datasets/schnabel_Tab3.1.xlsx'
 data=urllib.request.urlopen(url).read();local=Path('data/raw/McCreaMorgan/schnabel_Tab3.1.xlsx')
 (ROOT/local).parent.mkdir(parents=True,exist_ok=True);(ROOT/local).write_bytes(data)
 out.append({'package':'McCreaMorgan','url':url,'path':str(local),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'source':'McCrea and Morgan (2014), Table 3.1; attribute original collectors.'})
 (ROOT/'data/source_manifest.json').write_text(json.dumps({'sources':out},indent=2)+'\n')
 print(f'Preserved {len(out)} pinned source files, {sum(x["bytes"] for x in out)} bytes.')
if __name__=='__main__':main()
