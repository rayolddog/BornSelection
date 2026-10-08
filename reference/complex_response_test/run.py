"""Transient linear complex response, explicitly assumed rate readout and threshold comparator.
No sampled input site. Closed-form response; only hazard integrals use a time grid.
"""
from pathlib import Path
import sys,json,csv,gzip
import numpy as np
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'local_coupling_sweep'))
from run import spatial,x,SITES
sys.path.insert(0,str(P.parent/'microscopic_detector'))
from detector_model import goodness
T=6.;DT=.05;N=6000;RATE=.6;BASE=2026100500

def profile(kind):
 phase=1.4 if kind=='shifted' else .4
 ratio=0. if kind=='single' else (.4 if kind=='unequal' else 1.)
 a=np.exp(-x*x/(4*.85**2))*(np.exp(2j*np.pi*np.sin(np.pi/3)*x)+ratio*np.exp(1j*phase-2j*np.pi*np.sin(np.pi/3)*x))
 return a/np.linalg.norm(a)
CASES=[dict(name='uniform_'+k,profile=k) for k in ['double','shifted','single','unequal']]+[
 dict(name='slow_response',gamma=.25),dict(name='fast_response',gamma=2.),dict(name='fixed_electronic_phases',phase=True),
 dict(name='independent_weak',kind='independent',noise=.25),dict(name='independent_strong',kind='independent',noise=1.),
 dict(name='free_correlated',kind='free',noise=1.),dict(name='mirror_parallel',kind='parallel',noise=1.),dict(name='mirror_normal',kind='normal',noise=1.),
 dict(name='detuning_step',detuning='step'),dict(name='detuning_disorder',detuning='disorder')]

def response(t,l,omega=0.):
 return (np.exp(1j*omega*t)-np.exp(-l*t))/(l+1j*omega)

def simulate(case,idx,n=N,dt=DT,suffix=''):
 psi=profile(case.get('profile','double'));w=abs(psi)**2;gamma=case.get('gamma',1.);b=case.get('noise',0.)
 delta=np.zeros(SITES)
 if case.get('detuning')=='step':delta[x>=0]=2.
 if case.get('detuning')=='disorder':delta=np.random.default_rng(BASE+700).normal(0,.8,SITES)
 l=gamma+1j*delta
 L,var,check=spatial(case.get('kind','independent')) if b else (None,0.,{})
 rank=SITES if L is None else L.shape[1];rng=np.random.default_rng(BASE+idx)
 counts=np.zeros((2,SITES),int);eventrows=[];mean_probability=np.zeros(SITES);meansurvival=0.;timesall=[];thresholdtimes=[]
 checkpoints=[1.,3.,6.];responsecurve=np.zeros((round(T/dt)+1,));gatecounts=np.zeros((2,3),int)
 for start in range(0,n,100):
  size=min(100,n-start)
  def field():
   z=(rng.normal(size=(size,rank))+1j*rng.normal(size=(size,rank)))/np.sqrt(2)
   return z if L is None else z@L.T
  u,v=field(),field();cp=(u-1j*v)/2;cm=(u+1j*v)/2
  ephase=np.exp(1j*np.random.default_rng(BASE+701).uniform(-np.pi,np.pi,SITES)) if case.get('phase') else np.ones(SITES)
  clocks=rng.exponential(size=(size,SITES));hazard=np.zeros((size,SITES));oldpower=np.zeros((size,SITES));totalprob=np.zeros((size,SITES))
  sites=np.full((size,2),-1,int);ts=np.full((size,2),np.inf)
  for it in range(1,round(T/dt)+1):
   t=it*dt
   z=psi*ephase*(response(t,l)+b*(cp*response(t,l,1.)+cm*response(t,l,-1.)))
   power=abs(z)**2;responsecurve[it]+=power.sum()
   dh=RATE*(oldpower+power)*dt/2;newh=hazard+dh
   survival=np.exp(-hazard.sum(axis=1));dhall=dh.sum(axis=1)
   fractions=np.divide(dh,dhall[:,None],out=np.zeros_like(dh),where=dhall[:,None]>0)
   totalprob+=survival[:,None]*(-np.expm1(-dhall))[:,None]*fractions
   crossing=(hazard<clocks)&(newh>=clocks)&(sites[:,0,None]<0)
   frac=np.divide(clocks-hazard,dh,out=np.full_like(dh,np.inf),where=crossing)
   j=frac.argmin(axis=1);f=frac[np.arange(size),j];hit=np.isfinite(f)
   sites[hit,0]=j[hit];ts[hit,0]=t-dt+dt*f[hit]
   # Alternative rule with no separate clock RNG: first power threshold crossing.
   threshold=1/SITES;cross=(oldpower<threshold)&(power>=threshold)&(sites[:,1,None]<0)
   frac=np.divide(threshold-oldpower,power-oldpower,out=np.full_like(power,np.inf),where=cross)
   j=frac.argmin(axis=1);f=frac[np.arange(size),j];hit=np.isfinite(f)
   sites[hit,1]=j[hit];ts[hit,1]=t-dt+dt*f[hit]
   hazard=newh;oldpower=power
  mean_probability+=totalprob.sum(axis=0);meansurvival+=np.exp(-hazard.sum(axis=1)).sum()
  assert np.max(abs(totalprob.sum(axis=1)+np.exp(-hazard.sum(axis=1))-1))<2e-12
  for j in range(2):
   hit=sites[:,j]>=0;counts[j]+=np.bincount(sites[hit,j],minlength=SITES)
   for gi,gate in enumerate(checkpoints):gatecounts[j,gi]+=np.count_nonzero(ts[:,j]<=gate)
  for i in range(size):eventrows.append([start+i]+[v for j in range(2) for v in (int(sites[i,j]),float(ts[i,j]) if np.isfinite(ts[i,j]) else '')])
  timesall.extend(ts[sites[:,0]>=0,0]);thresholdtimes.extend(ts[sites[:,1]>=0,1])
 name=case['name']+suffix
 with gzip.open(P/(name+'_events.csv.gz'),'wt') as f:
  wr=csv.writer(f);wr.writerow(['trial','clock_site','clock_time','threshold_site','threshold_time']);wr.writerows(eventrows)
 mean_probability/=n;meansurvival/=n;responsecurve/=n
 result=dict(name=name,settings=case,incident=n,dt=dt,gamma=gamma,noise=b,spatial_variance=var,weights=w.tolist(),delta=delta.tolist(),clock_captures=int(counts[0].sum()),threshold_captures=int(counts[1].sum()),clock_counts=counts[0].tolist(),threshold_counts=counts[1].tolist(),clock_fit=goodness(counts[0],w,seed=BASE+idx+100,replicas=10000) if counts[0].sum()>=1000 else None,threshold_fit=goodness(counts[1],w,seed=BASE+idx+200,replicas=10000) if counts[1].sum()>=1000 else None,gate_counts=gatecounts.tolist(),mean_record_weights=mean_probability.tolist(),mean_survival=float(meansurvival),mean_response_curve=responsecurve.tolist(),clock_median=float(np.median(timesall)) if timesall else None,expected_captures_from_hazards=float(n*mean_probability.sum()))
 if not b and not case.get('detuning'):
  H=RATE/gamma**2*(T-2*(-np.expm1(-gamma*T))/gamma+(-np.expm1(-2*gamma*T))/(2*gamma))
  result['exact_capture_probability']=float(-np.expm1(-H));result['hazard_quadrature_capture_error']=abs(mean_probability.sum()+np.expm1(-H))
 print(name,result['clock_captures'],result['clock_fit']['p'] if result['clock_fit'] else None,result['threshold_captures'],flush=True)
 return result
if __name__=='__main__':
 cases=[simulate(c,i) for i,c in enumerate(CASES)]
 fits=[c[k] for c in cases for k in ['clock_fit','threshold_fit'] if c[k]];bound=0
 for rank,f in enumerate(sorted(fits,key=lambda f:f['p'])):
  bound=max(bound,min(1,(len(fits)-rank)*f['p']));f['holm_p']=bound
 (P/'results.json').write_text(json.dumps(dict(rate=RATE,gate=T,incident_total=N*len(cases),cases=cases),indent=2))
 # Identical seeds and batches preserve input fields and clocks for paired refinement.
 refs=[]
 for idx in [0,8,11,12]:
  runs=[simulate(CASES[idx],idx,n=200,dt=d,suffix=f'_check_{d}') for d in [.05,.025,.0125]]
  changes=[]
  for aa,bb in zip(runs,runs[1:]):
   logs=[]
   for rr in [aa,bb]:
    with gzip.open(P/(rr['name']+'_events.csv.gz'),'rt') as f:logs.append(list(csv.DictReader(f)))
   changes.append(dict(max_record_weight_change=float(np.max(abs(np.array(aa['mean_record_weights'])-bb['mean_record_weights']))),winner_changes=[sum(a[key]!=b[key] for a,b in zip(*logs)) for key in ['clock_site','threshold_site']]))
  refs.append(dict(case=CASES[idx]['name'],paired_histories=200,changes=changes))
 (P/'refinement.json').write_text(json.dumps(refs,indent=2))
