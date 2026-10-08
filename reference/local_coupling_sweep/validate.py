from pathlib import Path
import json,csv,gzip
import numpy as np
from run import spatial,w,psi,SITES,T,GATES
P=Path(__file__).resolve().parent
r=json.loads((P/'results.json').read_text());checks=[]
for c in r['results']:
 with gzip.open(P/(c['name']+'_events.csv.gz'),'rt') as f: rows=list(csv.reader(f))[1:]
 assert len(rows)==c['trials']
 for rule,s in enumerate(c['results']):
  for gi,gate in enumerate(GATES):
   sites=[int(a[1+2*rule]) for a in rows if a[2+2*rule] and float(a[2+2*rule])<=gate+1e-10]
   assert all(0<=j<SITES for j in sites)
   assert np.array_equal(np.bincount(sites,minlength=SITES),s['gate_site_counts'][gi])
  assert s['recorded']+s['unresolved']==c['trials']
base=r['results'][0];N=base['trials'];clock=base['results'][3]
for gate,n in zip(GATES,clock['gate_counts']):
 F=gate/T-np.sin(2*np.pi*gate/T)/(2*np.pi);p=1-np.exp(-3*F)
 checks.append(dict(gate=gate,clock_captures=n,expected=N*p,z=(n-N*p)/np.sqrt(N*p*(1-p))))
assert max(abs(c['z']) for c in checks)<4
peak=int(w.argmax());assert all(c==0 for j,c in enumerate(base['results'][0]['counts']) if j!=peak)
# Single global phase and exact darkness must not change scalar local transition weights.
rng=np.random.default_rng(101011);g=1+(rng.normal(size=(100,SITES))+1j*rng.normal(size=(100,SITES)))/np.sqrt(2)
phase_error=float(abs(abs(psi*g)**2-abs(psi*np.exp(.713j)*g)**2).max());assert phase_error<1e-14
psi0=psi.copy();psi0[13]=0;assert np.all(abs(psi0[None,:]*g)[:,13]==0)
# Check expected spatial covariance on a manageable subset, independent of full trajectories.
fieldchecks=[]
for kind in ['free','parallel','normal']:
 L,var,kcheck=spatial(kind);ids=np.array([0,40,80,120,160,200,255]);root=L[ids]
 z=(rng.normal(size=(50000,L.shape[1]))+1j*rng.normal(size=(50000,L.shape[1])))/np.sqrt(2)
 f=z@root.T;expect=root@root.T;emp=f.T@f.conj()/len(f)
 err=float(abs(emp-expect).max()/var);assert err<.025
 fieldchecks.append(dict(kind=kind,relative_covariance_max_error=err,**kcheck))
validation=dict(event_logs_verified=len(r['results']),event_histories=sum(c['trials'] for c in r['results']),normalization_error=abs(float(w.sum())-1),global_phase_error=phase_error,exact_dark_site_preserved=True,deterministic_brightest_site=peak,analytical_clock_capture_checks=checks,field_checks=fieldchecks)
(P/'validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation,indent=2))
