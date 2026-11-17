"""Reproduce classical-data summaries and certified fixed-count projections."""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import gzip,hashlib,json,math
import numpy as np
import pyreadr
from partial_population.projection import Design,Projector,design_key,verify_interval,verify_kernel_table
from partial_population.capture_recapture import lincoln_petersen,chapman,fixed_interval,chao_incidence,repeated_estimate
from partial_population.approximations import good_turing_normal
from partial_population.removals import removal_states,removal_estimates,removal_interval
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/empirical'

def read(pkg,name):
 f=next((ROOT/'data/raw'/pkg).glob(name+'.*'))
 return next(iter(pyreadr.read_r(str(f)).values()))

def savegz(path,obj):
 with gzip.open(path,'wt') as out:json.dump(obj,out,separators=(',',':'))

def run():
 OUT.mkdir(exist_ok=True)
 abundance={}
 for name in ('moths','birds'):
  a=read('sads',name).iloc[:,0].astype(int).tolist()
  abundance[name]={'counts':a,'n':sum(a),'observed':len(a),'fingerprint':dict(Counter(a)),'complete':True}
 a=read('vegan','BCI').sum(axis=0).astype(int).tolist()
 abundance['BCI']={'counts':a,'n':sum(a),'observed':len(a),'fingerprint':dict(Counter(a)),'complete':True,'sampling_note':'pooled tree counts; spatial census is not an iid token sample'}
 for pkg,name,rcol,fcol,n,D in [('VGAM','corbet','ofreq','species',9031,620),('vcdExtra','ShakeWords','count','Freq',884647,31534)]:
  df=read(pkg,name);fp={int(r):int(f) for r,f in zip(df[rcol],df[fcol])}
  abundance[name]={'n':n,'observed':D,'fingerprint':fp,'complete':False,'omitted_observed_types':D-sum(fp.values()),'recorded_frequency_max':max(fp)}
 rows=[];certs=[];tables=[]
 cached={}
 if (OUT/'real-data-certificates.json.gz').exists():
  with gzip.open(OUT/'real-data-certificates.json.gz','rt') as f:
   cached={(c['dataset'],c['design']['target'],c['design']['m']):c for c in json.load(f)}
 for name,data in abundance.items():
  fp=data['fingerprint'];n=data['n'];f1,f2=fp.get(1,0),fp.get(2,0)
  data['good_turing']=f1/n
  data['normal_approximations']={model:good_turing_normal(n,f1,f2,variance_model=model) for model in ('fixed_n','poisson')}
  data['chao1_point']=data['observed']+(f1*f1/(2*f2) if f2 else f1*(f1-1)/2)
  horizons=[10000,100000,884647] if name=='ShakeWords' else [n]
  for target,m in [('missing_mass',0)]+[('discoveries',m) for m in horizons]:
   design=Design(n,target,m)
   phi=[fp.get(r,0) for r in range(1,len(design.directions()[0][0])+1)]
   old=cached.get((name,target,m))
   if old and Design(**old['design'])==design and old['fingerprint']==phi and all('solver_status' in e for e in old['endpoints']):
    with gzip.open(OUT/old['table'],'rt') as f:table=json.load(f)
    model=Projector(design,table=table);cert=old
   else:
    model=Projector(design);cert=model.interval(phi)
   assert verify_kernel_table(model.table) and verify_interval(cert,model.table)
   key=design_key(design);savegz(OUT/f'kernel-{key}.json.gz',model.table)
   cert['dataset']=name;cert['table']=f'kernel-{key}.json.gz';certs.append(cert)
   gaps=[]
   for ep in cert['endpoints']:
    if ep.get('witness'):
     gaps.append(abs(ep['bound_unit']-ep['witness']['z'])*design.target_scale)
    else:gaps.append(None)
   rows.append({'dataset':name,'n':n,'target':target,'m':m,'interval':cert['interval'],'relaxation_endpoint_gap_bounds':gaps,'fallback':cert['fallback']})
   print(name,target,m,cert['interval'],flush=True)
 savegz(OUT/'real-data-certificates.json.gz',certs)
 # Full histories retain occasion order for repeated-capture estimates.
 histories={'hares':read('Rcapture','hare').to_numpy(dtype=int),
            'cottontails':np.array([[int(v) for v in h] for h in read('RMark','edwards.eberhardt')['ch']])}
 capture=[]
 for name,H in histories.items():
  K=H.shape[1];freq=dict(Counter(map(int,H.sum(axis=1))))
  seen=np.zeros(len(H),dtype=bool);marked=[];sizes=[];overlaps=[]
  for column in H.T:
   current=column.astype(bool)
   marked.append(int(seen.sum()));sizes.append(int(current.sum()));overlaps.append(int((seen&current).sum()));seen|=current
  # Zero catches and an empty marked pool contribute zero to both ratios.
  usable=[(a,b,j) for a,b,j in zip(marked,sizes,overlaps) if a>0 and b>0]
  marked,sizes,overlaps=map(list,zip(*usable))
  a,b,j=map(int,(H[:,0].sum(),H[:,1].sum(),(H[:,0]*H[:,1]).sum()))
  row={'dataset':name,'observed':len(H),'occasions':K,'frequency':freq,
       'schnabel':float(repeated_estimate(marked,sizes,overlaps)),
       'schumacher_eschmeyer':float(repeated_estimate(marked,sizes,overlaps,weight_marked=True)),
       'chao_incidence_point':float(chao_incidence(K,[freq.get(r,0) for r in range(1,K+1)])),
       'first_two_samples':[a,b,j],'first_two_exact_interval':fixed_interval(a,b,j),
       'histories':H.tolist()}
  if name=='cottontails':row['external_population_total']=135
  capture.append(row)
 rodli=read('Petersen','data_rodli');hist={str(h):int(f) for h,f in zip(rodli['cap_hist'],rodli['freq'])}
 a,b,j=hist['10']+hist['11'],hist['01']+hist['11'],hist['11']
 trout={'sample_sizes':[a,b],'recaptures':j,'lincoln_petersen':float(lincoln_petersen(a,b,j)),'chapman':float(chapman(a,b,j)),'exact_95_interval':fixed_interval(a,b,j)}
 gerking=read('fishmethods','Gerking')
 C,R,net=[gerking[col].astype(int).tolist() for col in ('C','R','nM')]
 removed=[c-r-v for c,r,v in zip(C[:-1],R[:-1],net[:-1])]
 assert removed==[0]*10+[2,0,1]
 states=removal_states(C,R,removed)
 assert [v[0] for v in states]==[0]+list(np.cumsum(net[:-1]))
 sunfish={'catches':C,'recaptures':R,'net_marked_column':net,'removed_after_except_last':removed,
          'marked_before':[s[0] for s in states],'cumulative_removals_before':[s[1] for s in states],
          'estimates':removal_estimates(C,R,removed),'initial_population_95_e_interval':removal_interval(C,R,removed),
          'note':'Final nM=0 is an unused package placeholder (Krebs table has a dash), not 16 inferred deaths. Only inter-occasion removals are reconstructed.'}
 raw=json.loads((ROOT/'data/processed/schnabel_table_raw.json').read_text())['Sheet1']
 collection=[]
 for row in raw[1:]:
  name=row[0];T=int(row[-1]);freq={i+1:int(v) for i,v in enumerate(row[1:-1]) if v is not None}
  # Wood mice have T=21 but columns 19..21 absent: preserve the omission.
  item={'dataset':name,'occasions':T,'frequency':freq,'recorded_observed':sum(freq.values()),'complete_frequency_table':T<=18,
        'chao_incidence_point':float(chao_incidence(T,[freq.get(r,0) for r in range(1,T+1)])) if T<=18 else None}
  if name.startswith('Taxicabs'):item['external_population_total']=420
  collection.append(item)
 sensitivity=[]
 for variant in json.loads((ROOT/'data/processed/frequency_variants.json').read_text())['variants']:
  base=next(r for r in collection if r['dataset']==variant['dataset']);f=variant['alternative_frequency']
  sensitivity.append({'dataset':variant['dataset'],'workbook_singletons':base['frequency'][1],
   'alternative_singletons':f[0],'workbook_observed':base['recorded_observed'],'alternative_observed':sum(f),
   'workbook_point':base['chao_incidence_point'],'alternative_point':float(chao_incidence(base['occasions'],f)),
   'source':variant['source'],'url':variant['url'],'reference_population':variant['reference_population']})
 # Independent cross-source agreement prevents orientation/transcription mistakes.
 for name,key in [('hares','Hares'),('cottontails','Rabbits')]:
  a=next(x for x in capture if x['dataset']==name)['frequency'];b=next(x for x in collection if x['dataset']==key)['frequency']
  assert {k:v for k,v in a.items() if v}=={k:v for k,v in b.items() if v}
 result={'abundance':abundance,'projections':rows,'capture_histories':capture,'trout':trout,'sunfish':sunfish,'capture_collection':collection,'frequency_version_sensitivity':sensitivity,
         'assumptions':'Projection coverage is unconditional under iid sampling with fixed known counts. Historical observational designs are not asserted iid; estimates are model-based illustrations.'}
 (OUT/'classical_results.json').write_text(json.dumps(result,indent=2)+'\n')
 print('Saved classical analysis and independently replayed all certificates.',flush=True)
if __name__=='__main__':run()
