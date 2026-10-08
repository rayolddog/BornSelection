"""Reproduce the fixed experiments and write all events and summary results."""
import csv,json,time
from pathlib import Path
from dataclasses import asdict,replace
import numpy as np
from detector_model import *
from test_model import run as validate

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
OUT.mkdir(exist_ok=True)


def save_events(name,e):
    keys=list(e)
    with (OUT/(name+'_events.csv')).open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['trial']+keys)
        for i in range(len(e['site'])):
            w.writerow([i+1]+[e[k][i] for k in keys])


def run(count=12000):
    start=time.time()
    checks=validate()
    (OUT/'validation.json').write_text(json.dumps(checks,indent=2))
    print('Validation passed',flush=True)
    p=Parameters()
    scenarios=[('single', 'single',p),('double','double',p),
      ('shifted','shifted',p),('unequal','unequal',p),
      ('incoherent','incoherent',p),
      ('coupling_phases','double',replace(p,random_coupling_phase=True)),
      ('thermal','double',replace(p,thermal_sigma=2.)),
      ('disorder','double',replace(p,disorder=True))]
    summary={'trials_per_case':count,'seed_policy':'Fixed seeds declared in source before results; no parameter fitting',
      'validation':checks,'quantum':{},'threshold':{},'playback':{}}
    for k,(name,kind,param) in enumerate(scenarios):
        e=quantum_trials(kind,param,count=count,seed=721000+k)
        counts=histogram(e,p.n); target=intensity(kind,p.n)
        reference=None
        if not param.thermal_sigma:
            reference=efficiency(param)
            target=target*reference['capture']; target/=target.sum()
        fit=goodness(counts,target,seed=84100+k)
        result={'kind':kind,'parameters':asdict(param),'counts':counts.tolist(),
          'reference':target.tolist(),'fit':fit,
          'outcomes':{str(i):int(np.count_nonzero(e['outcome']==i)) for i in range(4)},
          'mean_capture_time':float(e['time'][e['outcome']==1].mean()),
          'mean_avalanche_latency':float(np.nanmean(e['latency']))}
        if reference is not None:
            result['predicted_capture_efficiency']=float(np.sum(intensity(kind,p.n)*reference['capture']))
        if name=='disorder':
            result['uncorrected_born_fit']=goodness(counts,intensity(kind,p.n),seed=84999)
        summary['quantum'][name]=result
        summary['playback'][name]={'site':e['site'][:2000].tolist(),'outcome':e['outcome'][:2000].tolist(),
          'time':e['time'][:2000].round(3).tolist()}
        save_events(name,e)
        print(name, result['outcomes'], 'p=',fit['p'],flush=True)
    # Holm adjusted p across the eight prespecified quantum spatial cases.
    ordered=sorted(summary['quantum'],key=lambda key:summary['quantum'][key]['fit']['p'])
    running=0.
    for k,key in enumerate(ordered):
        running=max(running,min(1.,(len(ordered)-k)*summary['quantum'][key]['fit']['p']))
        summary['quantum'][key]['fit']['holm_p']=running
    tests=[('fixed','fixed',0.,False),('uniform','uniform',0.,False),
      ('exponential','exponential',0.,False),('additive_noise','exponential',.3,False),
      ('noise_dark','exponential',.3,True)]
    for k,(name,law,noise,vacuum) in enumerate(tests):
        e=threshold_trials(count=count,seed=93100+k,law=law,noise=noise,vacuum=vacuum)
        counts=histogram(e,p.n)
        result={'law':law,'threshold_mean_scale':.02,'noise':noise,'vacuum':vacuum,
          'counts':counts.tolist(),'captured':int(e['outcome'].sum()),
          'fraction_multiple_without_stop':float(np.mean(e['potential_count']>=2)),
          'mean_potential_crossings':float(np.mean(e['potential_count'])),
          'fit':goodness(counts,intensity('double',p.n),seed=89100+k) if not vacuum else None}
        summary['threshold'][name]=result
        save_events('threshold_'+name,e)
        print('threshold',name,'captures',result['captured'],'multiple',result['fraction_multiple_without_stop'],flush=True)
    # Independent-seed replication of the key race, plus a resolution check.
    summary['exponential_replications']=[]
    for seed in [61001,61002,61003]:
        e=threshold_trials(count=count,seed=seed)
        summary['exponential_replications'].append(goodness(histogram(e,p.n),intensity(),seed=seed+100))
    summary['elapsed_seconds']=time.time()-start
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    print('Saved',OUT/'summary.json','seconds',summary['elapsed_seconds'],flush=True)

if __name__=='__main__':
    run()
