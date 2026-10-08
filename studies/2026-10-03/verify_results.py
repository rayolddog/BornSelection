"""Verify saved results and independent analytic checks; run after simulate.py."""
from pathlib import Path
import csv, json
import numpy as np
from scipy.special import ndtr
from scipy.stats import chi2
D=Path(__file__).resolve().parent
rows=list(csv.DictReader(open(D/'basin_measures.csv')))
checks=list(csv.DictReader(open(D/'sensitivity.csv')))
assert len(rows)==287 and len(checks)==50
assert all(int(r['unresolved'])==0 for r in rows+checks)
# Hoeffding + union bound for 3*41 baseline cells and 2*41 analytic controls.
m=205; eps=np.sqrt(np.log(2*m/.05)/(2*8000))
report={'simultaneous_95pct_absolute_error_bound':float(eps),'gof':{}}
for name in dict.fromkeys(r['model'] for r in rows):
    a=[r for r in rows if r['model']==name]
    p=np.array([float(r['p']) for r in a]); f=np.array([float(r['frequency']) for r in a])
    assert np.all((f>=0)&(f<=1))
    pred=None
    if name.startswith('Unbiased'): pred=p
    if name.startswith('Amplitude'): pred=np.sqrt(p)/(np.sqrt(p)+np.sqrt(1-p))
    if name.startswith('Vacuum'): pred=ndtr((2*p-1)/.6)
    if pred is not None: assert max(abs(f-pred))<eps
    mask=(p>0)&(p<1); statistic=8000*np.sum((f[mask]-p[mask])**2/(p[mask]*(1-p[mask])))
    report['gof'][name]={'chi_square':float(statistic),'df':39,'p_value':float(chi2.sf(statistic,39))}
# Independent antithetic check of the exact boundary regulator's zero drift.
rng=np.random.default_rng(11); x=rng.uniform(size=10000); d=rng.normal(size=10000)*.07
step=lambda v:np.sign(v)*np.minimum(abs(v),np.minimum(x,1-x))
assert np.max(abs(step(d)+step(-d)))==0
assert all(0<=float(r['frequency'])<=1 for r in checks)
report['all_checks_passed']=True
(D/'validation.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
