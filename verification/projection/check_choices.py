"""Exact removal-law checks and reproduction of approximate calibration choices."""
import json,math
from pathlib import Path
from fractions import Fraction
from scipy.stats import norm
from partial_population.approximations import good_turing_normal,normal_interval
from partial_population.removals import removal_states,removal_evalue,removal_interval
from partial_population.capture_recapture import hypergeom_pmf
from partial_population.projection import Design,Projector,verify_kernel_table,verify_interval

def main():
 checks=0;coverages=[]
 for N in (5,8,12):
  sizes=[2,2,2];deaths=[1,1];leaves=[]
  def walk(t,M,D,hist,prob):
   if t==len(sizes):leaves.append((hist,prob));return
   for r,p in enumerate(hypergeom_pmf(N-D,M,sizes[t])):
    if p:
     d=deaths[t] if t<len(deaths) else 0
     walk(t+1,M+sizes[t]-r-d,D+d,hist+[r],prob*p)
  walk(0,0,0,[],Fraction(1));assert sum(p for _,p in leaves)==1
  for z in (Fraction(1,2),Fraction(2),Fraction(3,2)):
   total=sum((p*removal_evalue(N,removal_states(sizes,r,deaths),z) for r,p in leaves),Fraction(0))
   assert total==1;checks+=1
  coverage=Fraction(0)
  for r,p in leaves:
   lo,hi=removal_interval(sizes,r,deaths)
   coverage+=p*int(lo<=N and (hi is None or N<=hi))
  assert coverage>=Fraction(19,20);coverages.append(str(coverage));checks+=1
 # A known law with changes in available N and the marked pool.
 assert removal_states([3,2,2],[0,1,1],[1,1])==((0,0,3,0),(2,1,2,1),(2,2,2,1));checks+=1
 # Two normal variance models reproduce the analytic formulas, not certificates.
 n,f1,f2=100,20,5;z=norm.ppf(.975)
 fixed=good_turing_normal(n,f1,f2,clip=False);pois=good_turing_normal(n,f1,f2,variance_model='poisson',clip=False)
 assert fixed['variance']==26/10000 and pois['variance']==30/10000;checks+=1
 assert fixed['interval']==[.2-z*math.sqrt(26/10000),.2+z*math.sqrt(26/10000)];checks+=1
 assert good_turing_normal(10,0,0)['degenerate_variance'];checks+=1
 custom=normal_interval(.1,.0004,target='declared custom target',assumptions='user-supplied CLT and variance consistency',bias_bound=.01)
 assert math.isclose(custom['interval'][1],.11+z*.02);checks+=1
 # A different predeclared score uses the same continuum certificate machinery.
 d=Design(100,score_weights=('1','1/10'));model=Projector(d);cert=model.interval([20,5])
 assert verify_kernel_table(model.table) and verify_interval(cert,model.table);checks+=1
 result={'checks':checks,'exact_removal_coverages':coverages,
         'scope':'known-removal martingale expectations and interval inversion; normal-formula reproduction; custom-score continuum replay'}
 (Path(__file__).resolve().parent/'choice_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
