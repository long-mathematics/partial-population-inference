"""Predeclared, paired fixed-sample coverage/width experiment (seed 20261117)."""
import gzip,json,math
from pathlib import Path
from collections import defaultdict
from fractions import Fraction
import numpy as np
from scipy.stats import norm,beta
from partial_population.projection import Design,Projector,design_key,dual_endpoint,outward,verify_interval,verify_kernel_table
from partial_population.approximations import good_turing_normal
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'verification/empirical'

def savegz(path,obj):
 with gzip.open(path,'wt') as out:json.dump(obj,out,separators=(',',':'))

def run(reps=200):
 rng=np.random.default_rng(20261117)
 pops={'uniform100':np.full(100,.01),'geometric100':.95**np.arange(100),
       'rare_tail':np.r_[.8,np.full(200,.001)],'Rothamsted_reference':None}
 data=json.loads((OUT/'classical_results.json').read_text())
 pops['Rothamsted_reference']=np.array(data['abundance']['moths']['counts'],float)
 rows=[];certificates={};cache={};case_rows=[]
 if (OUT/'simulation-certificates.json.gz').exists():
  with gzip.open(OUT/'simulation-certificates.json.gz','rt') as f:certificates=json.load(f)
 for n in (100,500):
  models={target:Projector(Design(n,target,n if target=='discoveries' else 0)) for target in ('missing_mass','discoveries')}
  for model in models.values():
   assert verify_kernel_table(model.table)
   savegz(OUT/f'kernel-{design_key(model.design)}.json.gz',model.table)
  for name,p in pops.items():
   p=p/p.sum();results=defaultdict(list)
   for rep in range(reps):
    x=rng.multinomial(n,p);y=rng.multinomial(n,p)
    fp=[int(sum(x==r)) for r in range(1,4)]
    truths={'missing_mass':float(sum(p[x==0])),'discoveries':int(sum((x==0)&(y>0)))}
    for target,model in models.items():
     phi=fp[:1] if target=='missing_mass' else fp
     key=f'{n}-{target}-'+','.join(map(str,phi))
     if key not in certificates:
      cert=model.interval(phi);assert verify_interval(cert,model.table)
      cert['table']=f'kernel-{design_key(model.design)}.json.gz';certificates[key]=cert
     interval=certificates[key]['interval']
     choices={'certified_joint':interval,'full_support':[0,model.design.target_scale]}
     # A separate equal-budget tail construction on the same directions.
     splitkey=(n,target)
     if splitkey not in cache:
      offsets=[]
      for j in range(len(model.design.scales)):
       qi=[0]*len(model.design.scales);qi[j]=10**12
       v,positive=dual_endpoint(model.design,[0]*len(phi),qi,model.table)
       offsets.append((outward(v,positive),positive))
      cache[splitkey]=offsets
     # To preserve rigorous arithmetic, compute scalar shifts with Arb and
     # outward rounding; use the cached bound for the zero fingerprint.
     from flint import arb
     lows=[0.];ups=[1.]
     for (off,pos),(lam,eta) in zip(cache[splitkey],model.design.directions()):
      shift=sum((arb(str(l))*v for l,v in zip(lam,phi)),arb(0))/(arb(str(eta))*model.design.target_scale)
      v=outward(arb(off)-shift,pos)
      (ups if pos else lows).append(v)
     lo,hi=max(lows),min(ups)
     if lo>hi:lo,hi=0.,1.
     scale=model.design.target_scale
     lo=max(0,math.nextafter(lo*scale,-math.inf));hi=min(scale,math.nextafter(hi*scale,math.inf))
     if target=='discoveries':
      lo,hi=math.ceil(lo),math.floor(hi)
      if lo>hi:lo,hi=0,scale
     choices['certified_split']=[lo,hi]
     if target=='missing_mass':
      choices['GT_normal']=good_turing_normal(n,fp[0],fp[1])['interval']
      choices['GT_poisson_normal']=good_turing_normal(n,fp[0],fp[1],variance_model='poisson')['interval']
     for method,(lo,hi) in choices.items():
      covered=lo<=truths[target]<=hi
      results[(target,method)].append((covered,hi-lo))
     case_rows.append({'population':name,'n':n,'replicate':rep,'target':target,'truth':truths[target],
                       'certificate':key,'intervals':choices})
   for (target,method),values in results.items():
    hits=sum(int(v[0]) for v in values);widths=[v[1] for v in values]
    ci=[0. if hits==0 else float(beta.ppf(.025,hits,reps-hits+1)),
        1. if hits==reps else float(beta.ppf(.975,hits+1,reps-hits))]
    rows.append({'population':name,'n':n,'m':n if target=='discoveries' else 0,'target':target,'method':method,'replicates':reps,
                 'hits':hits,'coverage':hits/reps,'coverage_mc95':ci,'mean_width':float(np.mean(widths)),
                 'mean_width_mc_se':float(np.std(widths,ddof=1)/math.sqrt(reps))})
   print(n,name,'completed',len(certificates),'unique certified inputs',flush=True)
 savegz(OUT/'simulation-certificates.json.gz',certificates)
 savegz(OUT/'simulation-cases.json.gz',case_rows)
 (OUT/'coverage_width.json').write_text(json.dumps({'seed':20261117,'replicates_per_scenario':reps,'rows':rows,
  'scope':'Fixed n iid initial draws and independent fixed m iid future draws; paired methods; Rothamsted counts normalized to a specified empirical reference population, not the unknown natural community.'},indent=2)+'\n')
 print('Saved paired study and certificates.',flush=True)
if __name__=='__main__':run()
