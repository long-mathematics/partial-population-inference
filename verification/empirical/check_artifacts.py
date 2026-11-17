"""Replay continuum certificates, source hashes, simulated truths, and summary statistics.

Does not call a conic optimizer. Raw-data transformations are also reproduced by
scripts/classical_analysis.py; this check checks their essential sufficient statistics.
"""
import gzip,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from statistics import NormalDist
import numpy as np
import pyreadr
from scipy.stats import beta
from flint import arb
from partial_population.projection import (Design,verify_kernel_table,verify_interval,
                                            dual_endpoint,outward)
from partial_population.approximations import good_turing_normal
from partial_population.removals import removal_estimates,removal_interval
from partial_population.capture_recapture import chao_incidence

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'verification/empirical'
def gz(path):
 with gzip.open(path,'rt') as f:return json.load(f)
def read(pkg,name):
 return next(iter(pyreadr.read_r(str(next((ROOT/'data/raw'/pkg).glob(name+'.*')))).values()))
def close(a,b):assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),(a,b)

def main():
 sources=json.loads((ROOT/'data/source_manifest.json').read_text())['sources']
 for s in sources:
  b=(ROOT/s['path']).read_bytes()
  assert len(b)==s['bytes'] and hashlib.sha256(b).hexdigest()==s['sha256']
 classical=json.loads((OUT/'classical_results.json').read_text())
 g=read('fishmethods','Gerking');s=classical['sunfish']
 C,R,net=[g[col].astype(int).tolist() for col in ('C','R','nM')]
 d=[c-r-m for c,r,m in zip(C[:-1],R[:-1],net[:-1])]
 assert d==s['removed_after_except_last'] and C==s['catches'] and R==s['recaptures']
 assert removal_interval(C,R,d)==s['initial_population_95_e_interval']
 assert removal_estimates(C,R,d)==s['estimates']
 variants=json.loads((ROOT/'data/processed/frequency_variants.json').read_text())['variants']
 for v,row in zip(variants,classical['frequency_version_sensitivity']):
  f=v['alternative_frequency'];base=next(r for r in classical['capture_collection'] if r['dataset']==v['dataset'])
  assert row['alternative_observed']==sum(f) and row['alternative_singletons']==f[0]
  close(row['alternative_point'],float(chao_incidence(base['occasions'],f)))
 for a in classical['abundance'].values():
  for model in ('fixed_n','poisson'):
   assert a['normal_approximations'][model]==good_turing_normal(a['n'],a['fingerprint'].get('1',0),a['fingerprint'].get('2',0),variance_model=model)
 for name,pkg in [('moths','sads'),('birds','sads'),('BCI','vegan')]:
  df=read(pkg,name);counts=(df.sum(axis=0) if name=='BCI' else df.iloc[:,0]).astype(int).tolist()
  a=classical['abundance'][name]
  assert counts==a['counts'] and sum(counts)==a['n'] and len(counts)==a['observed']
  assert {str(k):v for k,v in Counter(counts).items()}==a['fingerprint']
 for pkg,name,r,f,n,D in [('VGAM','corbet','ofreq','species',9031,620),('vcdExtra','ShakeWords','count','Freq',884647,31534)]:
  df=read(pkg,name);a=classical['abundance'][name]
  assert a['fingerprint']=={str(int(x)):int(y) for x,y in zip(df[r],df[f])}
  assert (a['n'],a['observed'])==(n,D)
 real=gz(OUT/'real-data-certificates.json.gz');sim=gz(OUT/'simulation-certificates.json.gz')
 tables={}
 for name in sorted({c['table'] for c in real+list(sim.values())}):
  t=gz(OUT/name);verify_kernel_table(t);tables[name]=t
 print(f'{len(sources)} source hashes and {len(tables)} continuum tables verified.',flush=True)
 statuses=Counter();witnesses=0
 for i,c in enumerate(real+list(sim.values())):
  verify_interval(c,tables[c['table']])
  for ep in c['endpoints']:
   statuses[ep['solver_status']]+=1;witnesses+=bool(ep.get('witness'))
  if i and i%250==0:print(f'{i} certificates replayed.',flush=True)
 for c,row in zip(real,classical['projections']):
  assert c['dataset']==row['dataset'] and c['interval']==row['interval']
  d=Design(**c['design']);fp=classical['abundance'][c['dataset']]['fingerprint']
  assert c['fingerprint']==[fp.get(str(r),0) for r in range(1,len(c['fingerprint'])+1)]
 cases=gz(OUT/'simulation-cases.json.gz');report=json.loads((OUT/'coverage_width.json').read_text())
 pops={'uniform100':np.full(100,.01),'geometric100':.95**np.arange(100),
       'rare_tail':np.r_[.8,np.full(200,.001)],
       'Rothamsted_reference':np.array(classical['abundance']['moths']['counts'],float)}
 rng=np.random.default_rng(report['seed']);idx=0;groups=defaultdict(list);offsets={}
 for n in (100,500):
  for name,p in pops.items():
   p=p/p.sum()
   for rep in range(report['replicates_per_scenario']):
    x=rng.multinomial(n,p);y=rng.multinomial(n,p);fp=[int(sum(x==r)) for r in (1,2,3)]
    for target in ('missing_mass','discoveries'):
     case=cases[idx];idx+=1
     assert (case['population'],case['n'],case['replicate'],case['target'])==(name,n,rep,target)
     truth=float(sum(p[x==0])) if target=='missing_mass' else int(sum((x==0)&(y>0)))
     # NumPy powers and floating sums can differ by a few ulps across CPUs.
     # Integer discoveries, fingerprints, and coverage decisions stay exact.
     if target=='missing_mass':close(case['truth'],truth)
     else:assert case['truth']==truth
     c=sim[case['certificate']];d=Design(**c['design']);phi=fp[:1] if target=='missing_mass' else fp
     assert c['fingerprint']==phi and c['design']['n']==n and c['design']['target']==target
     assert case['intervals']['certified_joint']==c['interval']
     assert case['intervals']['full_support']==[0,d.target_scale]
     if (n,target) not in offsets:
      values=[]
      for j in range(len(d.scales)):
       q=[0]*len(d.scales);q[j]=10**12
       v,pos=dual_endpoint(d,[0]*len(phi),q,tables[c['table']]);values.append((outward(v,pos),pos))
      offsets[n,target]=values
     low,high=[0.],[1.]
     for (v,pos),(lam,eta) in zip(offsets[n,target],d.directions()):
      shift=sum((arb(str(l))*f for l,f in zip(lam,phi)),arb(0))/(arb(str(eta))*d.target_scale)
      (high if pos else low).append(outward(arb(v)-shift,pos))
     lo,hi=max(low),min(high)
     if lo>hi:lo,hi=0.,1.
     lo=max(0,math.nextafter(lo*d.target_scale,-math.inf));hi=min(d.target_scale,math.nextafter(hi*d.target_scale,math.inf))
     if target=='discoveries':
      lo,hi=math.ceil(lo),math.floor(hi)
      if lo>hi:lo,hi=0,d.target_scale
     assert case['intervals']['certified_split']==[lo,hi]
     if target=='missing_mass':
      radius=NormalDist().inv_cdf(.975)*math.sqrt(max(0,(fp[0]+2*fp[1]-fp[0]**2/n)/n**2))
      expected=[max(0,fp[0]/n-radius),min(1,fp[0]/n+radius)]
      for a,b in zip(case['intervals']['GT_normal'],expected):close(a,b)
      radius=NormalDist().inv_cdf(.975)*math.sqrt((fp[0]+2*fp[1])/n**2)
      expected=[max(0,fp[0]/n-radius),min(1,fp[0]/n+radius)]
      for a,b in zip(case['intervals']['GT_poisson_normal'],expected):close(a,b)
     for method,(lo,hi) in case['intervals'].items():
      groups[name,n,target,method].append((lo<=truth<=hi,hi-lo))
 assert idx==len(cases) and len(groups)==len(report['rows'])
 for row in report['rows']:
  values=groups[row['population'],row['n'],row['target'],row['method']]
  hits=sum(v[0] for v in values);R=len(values);widths=[v[1] for v in values]
  assert hits==row['hits'] and R==row['replicates'];close(hits/R,row['coverage'])
  close(float(np.mean(widths)),row['mean_width']);close(float(np.std(widths,ddof=1)/math.sqrt(R)),row['mean_width_mc_se'])
  ci=[0 if hits==0 else beta.ppf(.025,hits,R-hits+1),1 if hits==R else beta.ppf(.975,hits+1,R-hits)]
  for a,b in zip(ci,row['coverage_mc95']):close(a,b)
 result={'source_files':len(sources),'continuum_tables':len(tables),'real_data_certificates':len(real),
         'simulation_certificates':len(sim),'verified_primal_witnesses':witnesses,
         'endpoint_solver_statuses':dict(statuses),'paired_cases':len(cases),'summary_rows':len(groups),
         'status':'all source hashes, continuum certificates, simulated truths and reported summaries verified'}
 (OUT/'replay_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
