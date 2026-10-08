from pathlib import Path
import json,gzip,csv,sys
import numpy as np
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parent/'microscopic_detector'))
from detector_model import goodness
r=json.loads((P/'results.json').read_text());counts=[]
for c in r['cases']:
 with gzip.open(P/(c['name']+'_events.csv.gz'),'rt') as f:rows=list(csv.DictReader(f))
 assert len(rows)==c['incident']
 for ri,kind in enumerate(['clock','threshold']):
  loc=[int(a[kind+'_site']) for a in rows if int(a[kind+'_site'])>=0]
  assert np.array_equal(np.bincount(loc,minlength=256),c[kind+'_counts'])
  for gi,gate in enumerate([1.,3.,6.]):assert sum(bool(a[kind+'_time']) and float(a[kind+'_time'])<=gate for a in rows)==c['gate_counts'][ri][gi]
 assert abs(sum(c['mean_record_weights'])+c['mean_survival']-1)<1e-12
for c in r['cases']:
 if not c['settings'].get('detuning'):continue
 w=np.array(c['weights']);l=c['gamma']+1j*np.array(c['delta'])
 # Independent dense integration of deterministic instantaneous density lambda_j exp(-sum H).
 estimates=[]
 for steps in [6000,12000]:
  dt=6/steps;t=np.arange(steps+1)*dt;F=(1-np.exp(-t[:,None]*l))/l
  rates=.6*abs(F)**2*w
  H=np.zeros(steps+1);H[1:]=np.cumsum((rates[:-1].sum(axis=1)+rates[1:].sum(axis=1))*dt/2)
  density=rates*np.exp(-H)[:,None];pr=np.sum((density[:-1]+density[1:])*dt/2,axis=0);estimates.append(pr)
 ref=estimates[-1]/estimates[-1].sum()
 c['response_corrected_fit']=goodness(c['clock_counts'],ref,seed=82711,replicas=10000)
 c['response_reference_quadrature_change']=float(abs(estimates[0]-estimates[1]).max())
 c['coarse_hazard_reference_max_error']=float(abs(np.array(c['mean_record_weights'])-estimates[-1]).max())
 c['independent_response_prediction']=ref.tolist()
(P/'results.json').write_text(json.dumps(r,indent=2))
validation=dict(verified_main_histories=r['incident_total'],verified_main_logs=len(r['cases']),probability_conservation=True,max_equal_response_capture_error=max(c.get('hazard_quadrature_capture_error',0) for c in r['cases']),detuned_checks=[dict(case=c['name'],fit=c['response_corrected_fit'],quadrature_change=c['response_reference_quadrature_change'],coarse_error=c['coarse_hazard_reference_max_error']) for c in r['cases'] if 'response_corrected_fit' in c])
(P/'validation.json').write_text(json.dumps(validation,indent=2))
with (P/'counts.csv').open('w') as f:
 wr=csv.writer(f);wr.writerow(['case','site','incident_born_weight','clock_count','threshold_count'])
 for c in r['cases']:
  for j,w in enumerate(c['weights']):wr.writerow([c['name'],j,w,c['clock_counts'][j],c['threshold_counts'][j]])
print(validation)
