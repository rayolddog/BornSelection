"""Finite microscopic-site parameter study; phenomenological selection laws, not full QED.
Run from any directory. Requires NumPy and earlier benchmark/model helpers.
"""
from pathlib import Path
import sys,json,gzip,csv,time
import numpy as np
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'qed_vacuum_benchmark'))
sys.path.insert(0,str(P.parent/'microscopic_detector'))
from benchmark import kernel
from detector_model import goodness
SITES=256;TRIALS=4000;DT=.2;T=10.;THRESHOLDS=[.5,2.,8.];GATES=[2.,5.,10.]
x=(np.arange(SITES)-(SITES-1)/2)*.01
# Illustrative two crossing coherent beams, 50 nm wavelength, .5 nm site spacing.
psi=np.exp(-x*x/(4*.85**2))*(np.exp(2j*np.pi*np.sin(np.pi/3)*x)+np.exp(.4j)*np.exp(-2j*np.pi*np.sin(np.pi/3)*x))
psi/=np.linalg.norm(psi);w=abs(psi)**2
configs=[
 ('no_fluctuations','independent',0.,0.,1.),
 ('electronic_only','independent',1.,0.,1.),
 ('independent_fast','independent',.5,1.,.1),
 ('independent_medium','independent',.5,1.,1.),
 ('independent_slow','independent',.5,1.,5.),
 ('free_weak','free',.5,.25,1.),
 ('free_medium','free',.5,1.,1.),
 ('free_strong','free',.5,2.,1.),
 ('free_fast','free',.5,1.,.1),
 ('free_slow','free',.5,1.,5.),
 ('mirror_parallel','parallel',.5,1.,1.),
 ('mirror_normal','normal',.5,1.,1.)]
def spatial(kind):
 if kind=='independent':return None,1.,dict(rank=SITES,min_eigenvalue=1.,reconstruction_error=0.)
 h=.05 if kind in ['parallel','normal'] else None
 K=kernel(x,height=h,orientation='z' if kind=='normal' else 'y')
 ev,V=np.linalg.eigh(K);keep=ev>1e-10*ev.max();L=V[:,keep]*np.sqrt(ev[keep])
 return L,float(K[0,0]),dict(rank=int(keep.sum()),min_eigenvalue=float(ev.min()),reconstruction_error=float(abs(L@L.T-K).max()))
def simulate(cfg,number,trials=TRIALS,dt=DT,suffix=''):
 name,kind,e,v,tau=cfg;L,var,check=spatial(kind);rank=SITES if L is None else L.shape[1]
 rng=np.random.default_rng(2026100200+number); calibration=1+e*e+v*v*var
 def noise(n):return (rng.normal(size=(n,rank))+1j*rng.normal(size=(n,rank)))/np.sqrt(2)
 total_counts=np.zeros((4,3,SITES),int);events=[];examples=[]
 # Four counterfactual selectors use the same wave/field histories. Not four photons.
 for start in range(0,trials,200):
  n=min(200,trials-start);z=noise(n);phase=rng.uniform(-np.pi,np.pi,(n,SITES))
  found=np.full((4,n),-1,int);times=np.full((4,n),np.inf);previous=np.zeros((n,SITES));hazard=np.zeros((n,SITES))
  clocks=rng.exponential(size=(n,SITES)) # extra stochastic absorption assumption
  rho=np.exp(-dt/tau)
  for step in range(1,round(T/dt)+1):
   t=step*dt;z=rho*z+np.sqrt(1-rho*rho)*noise(n)
   phase+=np.sqrt(2*dt/tau)*rng.normal(size=(n,SITES))
   xi=z if L is None else z@L.T
   # Relative electronic pathways, not an overall random electronic state phase.
   g=(1+e*np.exp(1j*phase)+v*xi)/np.sqrt(calibration)
   # Full complex incident field retained locally. Optical paths were added before modulus.
   amplitude=psi[None,:]*g
   score=abs(amplitude)**2*np.sin(np.pi*t/T)**2
   for rule,theta in enumerate(THRESHOLDS):
    target=theta/SITES
    crossing=(previous<target)&(score>=target)
    frac=np.divide(target-previous,score-previous,out=np.ones_like(score),where=crossing)
    ct=np.where(crossing,t-dt+dt*frac,np.inf);j=ct.argmin(axis=1);first=ct[np.arange(n),j]
    hit=(found[rule]<0)&np.isfinite(first);found[rule,hit]=j[hit];times[rule,hit]=first[hit]
   # Rate=3 * normalized pulse intensity * local coupling weight.
   increment=3*(previous+score)/T*dt # trapezoid of 2/T sin^2 pulse
   new=hazard+increment;crossing=(hazard<clocks)&(new>=clocks)
   frac=np.divide(clocks-hazard,increment,out=np.ones_like(score),where=crossing)
   ct=np.where(crossing,t-dt+dt*frac,np.inf);j=ct.argmin(axis=1);first=ct[np.arange(n),j]
   hit=(found[3]<0)&np.isfinite(first);found[3,hit]=j[hit];times[3,hit]=first[hit]
   hazard=new;previous=score
   if start==0 and step in [5,25,50]:examples.append(dict(time=t,score=score[0].tolist(),field_real=xi[0].real.tolist(),field_imag=xi[0].imag.tolist()))
  for rule in range(4):
   for gi,gate in enumerate(GATES):
    select=(found[rule]>=0)&(times[rule]<=gate+1e-10)
    total_counts[rule,gi]+=np.bincount(found[rule,select],minlength=SITES)
  for i in range(n):events.append([start+i]+[value for rule in range(4) for value in (int(found[rule,i]),float(times[rule,i]) if np.isfinite(times[rule,i]) else '')])
 with gzip.open(P/(name+suffix+'_events.csv.gz'),'wt') as f:
  wr=csv.writer(f);wr.writerow(['history']+[a for s in ['threshold_0.5','threshold_2','threshold_8','absorption_clock'] for a in [s+'_site',s+'_time']]);wr.writerows(events)
 results=[]
 for rule in range(4):
  counts=total_counts[rule,-1];captured=int(counts.sum())
  fit=goodness(counts,w,seed=number*10+rule+90000,replicas=10000) if captured>=1000 else None
  results.append(dict(rule=['threshold_0.5','threshold_2','threshold_8','absorption_clock'][rule],recorded=captured,unresolved=trials-captured,counts=counts.tolist(),gate_counts=[int(a.sum()) for a in total_counts[rule]],gate_site_counts=total_counts[rule].tolist(),fit=fit))
 result=dict(name=name+suffix,kind=kind,electronic_strength=e,vacuum_strength=v,correlation_time=tau,dt=dt,trials=trials,spatial=check,vacuum_variance=var,normalization=calibration,results=results,examples=examples)
 print(name+suffix,[(a['rule'],a['recorded'],round(a['fit']['tv'],4) if a['fit'] else None) for a in results],flush=True)
 return result
if __name__=='__main__':
 results=[]
 for k,c in enumerate(configs):
  results.append(simulate(c,k));(P/'progress.json').write_text(json.dumps({'completed':len(results),'total':len(configs)}))
 # All 48 final-gate tests form one family; only >=1000-record tests are formal.
 tests=[s for r in results for s in r['results'] if s['fit']]
 bound=0
 for rank,s in enumerate(sorted(tests,key=lambda a:a['fit']['p'])):
  bound=max(bound,min(1,(len(tests)-rank)*s['fit']['p']));s['fit']['holm_p']=bound
 (P/'results.json').write_text(json.dumps(dict(sites=SITES,x_lambda=x.tolist(),weights=w.tolist(),psi_real=psi.real.tolist(),psi_imag=psi.imag.tolist(),gates=GATES,thresholds=THRESHOLDS,pulse_duration=T,dt=DT,results=results),indent=2))
 with (P/'counts.csv').open('w') as f:
  wr=csv.writer(f);wr.writerow(['case','rule','gate','site','count','born_weight'])
  for r in results:
   for s in r['results']:
    for gi,gate in enumerate(GATES):
     for j,c in enumerate(s['gate_site_counts'][gi]):wr.writerow([r['name'],s['rule'],gate,j,c,w[j]])
