"""Create or independently replay a certified missing-mass/novelty interval."""
import argparse,gzip,json
from pathlib import Path
from partial_population.projection import Design,Projector,verify_kernel_table,verify_interval
from partial_population.approximations import good_turing_normal

def read(path):
 with gzip.open(path,'rt') as f:return json.load(f)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--n',type=int)
 p.add_argument('--future',type=int,default=0,help='fixed future sample size; zero selects missing mass')
 p.add_argument('--fingerprint',type=int,nargs='+',help='Phi_1 for mass; Phi_1 Phi_2 Phi_3 for discoveries')
 p.add_argument('--alpha',default='1/20',help='rational significance level')
 p.add_argument('--normal',choices=('fixed_n','poisson'),help='approximate Good--Turing calibration; requires Phi_1 Phi_2')
 p.add_argument('--score-weights',nargs='+',default=(),help='predeclared rational weights for a linear fingerprint score')
 p.add_argument('--output',type=Path,default=Path('build/projection'))
 p.add_argument('--verify',type=Path,help='replay an output prefix without running the optimizer')
 a=p.parse_args()
 if a.normal:
  from fractions import Fraction
  if a.future or a.verify or a.score_weights:p.error('--normal is for mass and cannot be combined with --verify or --score-weights')
  if a.n is None or a.fingerprint is None or len(a.fingerprint)!=2:p.error('--normal requires --n and two fingerprint values')
  result=good_turing_normal(a.n,*a.fingerprint,variance_model=a.normal,alpha=float(Fraction(a.alpha)))
  a.output.parent.mkdir(parents=True,exist_ok=True)
  Path(str(a.output)+'.normal.json').write_text(json.dumps(result,indent=2)+'\n')
  print(json.dumps(result,indent=2));return
 prefix=a.verify or a.output
 certpath=Path(str(prefix)+'.certificate.json.gz');tablepath=Path(str(prefix)+'.kernel.json.gz')
 if a.verify:
  cert,table=read(certpath),read(tablepath)
 else:
  if a.n is None or a.fingerprint is None:p.error('--n and --fingerprint are required')
  d=Design(a.n,'discoveries' if a.future else 'missing_mass',a.future,a.alpha,score_weights=tuple(a.score_weights))
  model=Projector(d);cert,table=model.interval(a.fingerprint),model.table
 verify_kernel_table(table);verify_interval(cert,table)
 if not a.verify:
  prefix.parent.mkdir(parents=True,exist_ok=True)
  for path,obj in ((certpath,cert),(tablepath,table)):
   with gzip.open(path,'wt') as f:json.dump(obj,f,separators=(',',':'))
 print(json.dumps({'interval':cert['interval'],'target':cert['design']['target'],
                   'alpha':cert['design']['alpha'],'verified':True,'certificate':str(certpath)},indent=2))

if __name__=='__main__':main()
