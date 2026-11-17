"""Generate manuscript tables and a scientific coverage/width figure from checked JSON."""
import argparse,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'manuscripts/generated'
NAMES={'moths':'Rothamsted moths','birds':'Preston birds','BCI':'BCI trees','corbet':'Corbet butterflies','ShakeWords':'Shakespeare',
       'uniform100':'Uniform','geometric100':'Geometric','rare_tail':'Rare tail','Rothamsted_reference':'Moth reference'}
METHODS={'certified_joint':'Joint','certified_split':'Split','GT_normal':'GT fixed','GT_poisson_normal':'GT Poisson','full_support':'Support'}
def line(items):return ' & '.join(map(str,items))+r' \\'+'\n'
def gap(v):return '--' if v is None else f'{v:.4g}'
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args()
 a=json.loads((ROOT/'verification/empirical/classical_results.json').read_text())
 s=json.loads((ROOT/'verification/empirical/coverage_width.json').read_text());texts={}
 texts['abundance.tex']=''.join(line([NAMES[k],f"{v['n']:,}",f"{v['observed']:,}",v['fingerprint'].get('1',0),v['fingerprint'].get('2',0),f"{100*v['good_turing']:.3f}",f"{v['chao1_point']:.1f}"]) for k,v in a['abundance'].items())
 mass=[];future=[];gaps=[]
 for row in a['projections']:
  lo,hi=row['interval'];name=NAMES[row['dataset']]
  if row['target']=='missing_mass':
   # Display percentages with outward rounding, preserving the interval claim.
   interval=f'$[{math.floor(lo*10000)/100:.2f}, {math.ceil(hi*10000)/100:.2f}]$'
   mass.append(line([name,interval,*[gap(v) for v in row['relaxation_endpoint_gap_bounds']]]))
  else:
   future.append(line([name,f"{row['m']:,}",f'$[{lo:,}, {hi:,}]$',*[gap(v) for v in row['relaxation_endpoint_gap_bounds']]]))
 texts['missing_mass.tex']=''.join(mass);texts['future.tex']=''.join(future)
 texts['histories.tex']=''.join(line([r['dataset'].capitalize(),r['observed'],r['occasions'],f"{r['schnabel']:.3f}",f"{r['schumacher_eschmeyer']:.3f}",f"{r['chao_incidence_point']:.3f}"]) for r in a['capture_histories'])
 def percent_interval(bounds):
  lo,hi=bounds
  return f'$[{math.floor(lo*100000)/1000:.3f}, {math.ceil(hi*100000)/1000:.3f}]$'
 texts['normal_choices.tex']=''.join(line([NAMES[k],percent_interval(v['normal_approximations']['fixed_n']['interval']),percent_interval(v['normal_approximations']['poisson']['interval'])]) for k,v in a['abundance'].items())
 texts['source_sensitivity.tex']=''.join(line([r['dataset'],f"{r['workbook_singletons']} / {r['alternative_singletons']}",f"{r['workbook_observed']} / {r['alternative_observed']}",f"{r['workbook_point']:.3f} / {r['alternative_point']:.3f}"]) for r in a['frequency_version_sensitivity'])
 war=json.loads((ROOT/'data/processed/wartime_comparison.json').read_text())
 texts['wartime.tex']=''.join(line([r['month'],r['intelligence_estimate'],r['statistical_estimate'],r['german_records'],f"{100*(r['statistical_estimate']/r['german_records']-1):+.1f}"]) for r in war['rows'])
 texts['collection.tex']=''.join(line([r['dataset'],r['occasions'],str(r['recorded_observed'])+('$^*$' if not r['complete_frequency_table'] else ''),r['frequency'].get('1',0),r['frequency'].get('2',0),gap(r['chao_incidence_point']),r.get('external_population_total','--')]) for r in a['capture_collection'])
 for target in ('missing_mass','discoveries'):
  rows=[]
  for r in s['rows']:
   if r['target']!=target:continue
   low,high=r['coverage_mc95']
   width=(f"{r['mean_width']:.4f} ({r['mean_width_mc_se']:.4f})" if target=='missing_mass' else f"{r['mean_width']:.2f} ({r['mean_width_mc_se']:.2f})")
   rows.append(line([r['n'],NAMES[r['population']],METHODS[r['method']],f"{r['coverage']:.3f} [{low:.3f}, {high:.3f}]",width]))
  texts[f'simulation_{target}.tex']=''.join(rows)
 OUT.mkdir(exist_ok=True)
 for name,content in texts.items():
  if args.check:assert (OUT/name).read_text()==content,f'Stale generated table: {name}'
  else:(OUT/name).write_text(content)
 if args.check:
  print(f'{len(texts)} manuscript tables match the checked JSON.');return
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 import numpy as np
 plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 fig,axes=plt.subplots(2,2,figsize=(7.3,5.7),layout='constrained')
 pops=list(NAMES)[5:];colors={'certified_joint':'#126879','certified_split':'#ce702d','GT_normal':'#77459a','GT_poisson_normal':'#537528'}
 for i,target in enumerate(('missing_mass','discoveries')):
  methods=['certified_joint','certified_split']+(['GT_normal','GT_poisson_normal'] if i==0 else [])
  for j,method in enumerate(methods):
   rows=[next(r for r in s['rows'] if r['n']==500 and r['target']==target and r['method']==method and r['population']==pop) for pop in pops]
   x=np.arange(4)+(j-(len(methods)-1)/2)*.16
   cov=np.array([r['coverage'] for r in rows]);err=np.array([[r['coverage']-r['coverage_mc95'][0] for r in rows],[r['coverage_mc95'][1]-r['coverage'] for r in rows]])
   axes[i,0].errorbar(x,cov,yerr=err,fmt='o',markersize=4,capsize=2,color=colors[method],label=METHODS[method])
   scale=1 if i==0 else 500
   axes[i,1].errorbar(x,[r['mean_width']/scale for r in rows],yerr=[1.96*r['mean_width_mc_se']/scale for r in rows],fmt='o',markersize=4,capsize=2,color=colors[method])
  axes[i,0].axhline(.95,color='gray',linestyle='--',linewidth=.8);axes[i,0].set_ylim(.85,1.025)
  axes[i,0].set_ylabel(('Missing mass' if i==0 else 'Future discoveries')+'\nCoverage')
  axes[i,1].set_ylabel('Mean width'+(' / m' if i else ''));axes[i,1].set_ylim(bottom=0)
  for ax in axes[i]:
   ax.set_xticks(range(4),['Uniform','Geometric','Rare tail','Moths']);ax.grid(axis='y',alpha=.2)
 axes[0,0].set_title('Coverage with exact 95% MC intervals');axes[0,1].set_title('Width with 1.96 MC standard errors')
 axes[0,0].legend(loc='lower left',fontsize=7,frameon=False,ncol=2)
 fig.savefig(OUT/'coverage_width.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
 print(f'Generated {len(texts)} tables and coverage_width.pdf.')
if __name__=='__main__':main()
