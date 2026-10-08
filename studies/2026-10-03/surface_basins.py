"""Phenomenological surface-field/thermal threshold basins, not a QED event model.
Run: python3 surface_basins.py (numpy scipy matplotlib). No Born draws.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/born-mpl')
import numpy as np
from scipy.signal import lfilter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv,json
D=Path(__file__).resolve().parent
P=np.linspace(0,1,21); N=4000; dt=.04; times=np.arange(0,12+dt/2,dt)
GATES=[2.,4.,6.,8.,12.]; QS=[.15,.35,.7]; SEED=20261004
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})

def surface(n,ts,rng,res=2.,common=0.):
    # Positive-frequency classical analytic Gaussian surrogate, two contacts.
    # A passive Lorentz dielectric + quasistatic reflected contribution.
    nf=2048; period=81.92; w=2*np.pi*np.arange(nf)/period
    eps=1+2*res**2/(res**2-w*w-1j*.3*w)
    rp=(eps-1)/(eps+1)
    S=(.02*w**3+np.maximum(rp.imag,0))*np.exp(-w/6);S[0]=0
    S[w>12]=0; S/=S.sum()
    # Coupled absorber response to each Fourier frequency; carrier = 2.
    H=-1j/(.5+1j*(w-2))
    idx=np.rint(ts/(period/nf)).astype(int)
    # dt=.02 refinement evaluated by direct half-grid interpolation using nf=4096.
    if abs(ts[1]-.02)<1e-8:
        nf=4096; w=2*np.pi*np.arange(nf)/period
        eps=1+2*res**2/(res**2-w*w-1j*.3*w);rp=(eps-1)/(eps+1)
        S=(.02*w**3+np.maximum(rp.imag,0))*np.exp(-w/6);S[0]=0;S[w>12]=0;S/=S.sum()
        H=-1j/(.5+1j*(w-2));idx=np.rint(ts/(period/nf)).astype(int)
    values=[]
    for channel in range(2):
        result=np.empty((len(ts),n),complex)
        for start in range(0,n,200):
            m=min(200,n-start)
            coeff=(rng.normal(size=(m,nf))+1j*rng.normal(size=(m,nf)))/np.sqrt(2)
            z=np.fft.ifft(coeff*np.sqrt(S)*H,axis=1)*nf
            result[:,start:start+m]=z[:,idx].T*np.exp(-2j*ts[:,None])
        values.append(result)
    if common: values[1]=common*values[0]+np.sqrt(1-common**2)*values[1]
    return .12*np.array(values),w,S,H

def environment(n,ts,seed,res=2.,common=0.):
    rng=np.random.default_rng(seed);d=ts[1]-ts[0]
    vac,w,S,H=surface(n,ts,rng,res,common)
    a=np.exp(-.5*d); thermal=np.empty((2,len(ts),n),complex)
    for j in range(2):
        eta=(rng.normal(size=(len(ts),n))+1j*rng.normal(size=(len(ts),n)))/np.sqrt(2)
        eta*=.08*np.sqrt(1-a*a)
        eta[0]=(rng.normal(size=n)+1j*rng.normal(size=n))*.08/np.sqrt(2)
        thermal[j]=lfilter([1],[1,-a],eta,axis=0)
    delay=rng.normal(0,.18,n);angle=rng.normal(0,.05,n)
    return vac,thermal,delay,angle,(w,S,H)

def signal(ts,delay):
    d=ts[1]-ts[0];f=np.exp(-(ts[:,None]-3-delay[None,:])**2/(4*.45**2))
    f/=np.sqrt(np.trapezoid(abs(f)**2,ts,axis=0))[None,:]
    a=np.exp(-.5*d)
    # Piecewise-constant midpoint drive; exponential update, second-order global.
    z=np.zeros_like(f,dtype=complex)
    for k in range(1,len(ts)):
        z[k]=a*z[k-1]-1j*np.sqrt(.8)*(1-a)/.5*(f[k]+f[k-1])/2
    return z

def crossings(z,q,ts):
    intensity=abs(z)**2
    accum=np.vstack([np.zeros(z.shape[1]),np.cumsum((intensity[:-1]+intensity[1:])*.5*(ts[1]-ts[0]),axis=0)])
    mask=accum>=q;has=mask.any(axis=0);idx=mask.argmax(axis=0);idx=np.maximum(idx,1)
    cols=np.arange(z.shape[1]);prev=accum[idx-1,cols];cur=accum[idx,cols]
    frac=np.divide(q-prev,cur-prev,out=np.zeros_like(prev),where=cur>prev)
    return np.where(has,ts[idx-1]+(ts[1]-ts[0])*frac,np.inf)

rows=[];maps={};prepared={}; diagnostics={}
configs=[('No fluctuations',False,False,False,0,2.,0.),('Source variations',True,False,False,0,2.,0.),('Surface vacuum surrogate',False,True,False,0,2.,0.),('Thermal equilibrium noise',False,False,True,0,2.,0.),('Combined',True,True,True,0,2.,0.),('Combined correlated surface',True,True,True,.6,2.,0.),('Combined shifted surface',True,True,True,0,3.,0.),('Combined biased thresholds',True,True,True,0,2.,.12)]
base=environment(N,times,SEED)

def run_case(name,src,usevac,useth,common,res,bias,n=N,ts=times,env=None,ps=P,thresholds=QS,tag='primary',phi=0.):
    if env is None:env=environment(n,ts,SEED+17,res,common)
    vac,thermal,delay,angle,spectral=env
    if not src:delay=np.zeros(n);angle=np.zeros(n)
    sig=signal(ts,delay)
    noise=(vac if usevac else np.zeros_like(vac))+(thermal if useth else 0)
    # Record background field is retained: no post-hoc vacuum subtraction.
    curves=[]
    for p in list(ps)+[None]:
        if p is None:c1=np.zeros(n);c2=np.zeros(n);peff=0
        else:
            c1=np.cos(angle)*np.sqrt(p)+np.sin(angle)*np.sqrt(1-p)*np.exp(1j*phi)
            c2=-np.sin(angle)*np.sqrt(p)+np.cos(angle)*np.sqrt(1-p)*np.exp(1j*phi)
            peff=np.mean(abs(c1)**2)
        z1=sig*c1[None,:]+noise[0];z2=sig*c2[None,:]+noise[1]
        for q in thresholds:
            t1=crossings(z1,q*(1-bias),ts);t2=crossings(z2,q*(1+bias),ts)
            for gate in GATES:
                a=t1<=gate;b=t2<=gate
                # All categories retained. Near-simultaneous crossings are ties.
                tie=a&b&np.isclose(t1,t2,atol=1e-10,rtol=0)
                first1=a&(~b|(t1<t2))&~tie
                first2=b&(~a|(t2<t1))&~tie
                none=~(a|b)
                assert np.all(first1.astype(int)+first2+tie+none==1)
                row=dict(run=tag,model=name,p='dark' if p is None else float(p),prepared_p=peff,n=n,q=q,gate=gate,first1=int(first1.sum()),first2=int(first2.sum()),tie=int(tie.sum()),none=int(none.sum()),only1=int((a&~b).sum()),only2=int((b&~a).sum()),both=int((a&b).sum()))
                rows.append(row)
                if tag=='primary' and name=='Combined' and q==.35 and p is not None:
                    maps[(float(p),gate)]=np.where(first1,1,np.where(first2,2,np.where(tie,3,0)))
    print(name,tag,'done',flush=True)

for cfg in configs:
    env=base if cfg[4]==0 and cfg[5]==2 else None
    run_case(*cfg,env=env)
# Refinement with independent ensembles at the same physical spectra/rates.
for d in [.04,.02]:
    ts=np.arange(0,12+d/2,d); env=environment(2000,ts,SEED+91)
    run_case('Combined',True,True,True,0,2.,0.,n=2000,ts=ts,env=env,ps=[.1,.25,.5,.75,.9],thresholds=[.35],tag=f'dt={d}')
phase_env=environment(2000,times,SEED+123,common=.6)
for phi in [0,np.pi/2,np.pi]:
    run_case('Combined correlated surface',True,True,True,.6,2.,0.,n=2000,env=phase_env,ps=[.1,.25,.5,.75,.9],thresholds=[.35],tag=f'phase={phi}',phi=phi)
# Pointwise Wilson intervals; shared histories correlate different rows.
for r in rows:
    f=r['first1']/r['n']; z=1.95996398454; den=1+z*z/r['n']
    center=(f+z*z/(2*r['n']))/den; half=z*np.sqrt(f*(1-f)/r['n']+z*z/(4*r['n']**2))/den
    r['first1_ci_low']=center-half; r['first1_ci_high']=center+half
with open(D/'surface_basin_counts.csv','w') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)

def select(name,q=.35,gate=12,tag='primary',phi=0.):
    return [r for r in rows if r['model']==name and r['q']==q and r['gate']==gate and r['run']==tag and r['p']!='dark']
def summary(name):
    a=select(name);p=np.array([r['p'] for r in a]);f=np.array([r['first1']/r['n'] for r in a]);eff=np.array([1-r['none']/r['n'] for r in a]);den=np.array([r['first1']+r['first2'] for r in a]);cond=np.divide([r['first1'] for r in a],den,out=np.full(len(a),np.nan),where=den>0)
    dark=next(r for r in rows if r['model']==name and r['p']=='dark' and r['q']==.35 and r['gate']==12 and r['run']=='primary')
    return dict(unconditional_rmse=float(np.sqrt(np.mean((f-p)**2))),conditional_rmse=float(np.sqrt(np.nanmean((cond-p)**2))),mean_any_record=float(eff.mean()),mean_both_if_no_latch=float(np.mean([r['both']/r['n'] for r in a])),mean_tie=float(np.mean([r['tie']/r['n'] for r in a])),dark_any_record=1-dark['none']/dark['n'])
summ={cfg[0]:summary(cfg[0]) for cfg in configs}
(D/'surface_basin_summary.json').write_text(json.dumps(summ,indent=2));print(json.dumps(summ,indent=2),flush=True)
fig,axs=plt.subplots(2,2,figsize=(12,9),layout='constrained')
for name in [c[0] for c in configs[:5]]:
    a=select(name);axs[0,0].plot(P,[r['first1']/N for r in a],'.-',label=name)
axs[0,0].plot(P,P,'k--',label='Born target');axs[0,0].set(xlabel='Nominal input |c₁|²',ylabel='Unconditional first-record-1 frequency',title='A. Noise changes the basins; no universal Born curve');axs[0,0].legend(fontsize=7)
for gate in GATES:
    a=select('Combined',gate=gate);axs[0,1].plot(P,[r['first1']/N for r in a],'.-',label=f'T={gate:g}')
axs[0,1].plot(P,P,'k--');axs[0,1].set(xlabel='Nominal input |c₁|²',ylabel='Unconditional first-record-1 frequency',title='B. The observation gate changes the measured basins');axs[0,1].legend(fontsize=8)
a=select('Combined');x=P
for key,label in [('only1','Only channel 1'),('only2','Only channel 2'),('both','Both cross (without latch)'),('none','Neither crosses')]:axs[1,0].plot(x,[r[key]/N for r in a],'.-',label=label)
axs[1,0].set(xlabel='Nominal input |c₁|²',ylabel='Fraction of all incident trials',title='C. Independent detector outcomes, T=12');axs[1,0].legend(fontsize=8)
for q in QS:
    a=select('Combined',q=q);den=np.array([r['first1']+r['first2'] for r in a]);cond=np.divide([r['first1'] for r in a],den,out=np.full(len(a),np.nan),where=den>0);axs[1,1].plot(x,cond,'.-',label=f'Threshold {q}')
axs[1,1].plot(P,P,'k--');axs[1,1].set(xlabel='Nominal input |c₁|²',ylabel='Channel-1 fraction among unique first records',title='D. Conditional curves depend on threshold');axs[1,1].legend(fontsize=8)
fig.suptitle('Complex excitation + surface surrogate + thermal bath → threshold basins\n4,000 noise histories per state; all misses, ties and double crossings retained',fontsize=14)
fig.savefig(D/'surface_basins.png',dpi=160);plt.close(fig)
from matplotlib.colors import ListedColormap,BoundaryNorm
cmap=ListedColormap(['#e4e8ec','#2189b8','#e68b39','#9a59b5'])
fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained')
for ax,gate in zip(axs,[4.,6.,12.]):
    mat=np.array([maps[(float(p),gate)][:180] for p in P]);ax.imshow(mat,origin='lower',aspect='auto',extent=[0,180,0,1],cmap=cmap,norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4),interpolation='nearest');ax.set(title=f'Observation gate T={gate:g}',xlabel='Same microscopic history',ylabel='Nominal |c₁|²')
fig.suptitle('Time-dependent basins: blue = first 1; orange = first 2; gray = no record; purple = tie')
fig.savefig(D/'surface_basin_maps.png',dpi=160);plt.close(fig)
# Finite-mode covariance, equilibrium variance, count arithmetic checks.
w,S,H=base[-1];pred=.12**2*np.sum(S*abs(H)**2);meas=float(np.mean(abs(base[0])**2))
validation={'vacuum_filtered_variance_predicted':float(pred),'vacuum_filtered_variance_measured':meas,'thermal_variance_target':.08**2,'thermal_variance_measured':float(np.mean(abs(base[1])**2)),'primary_histories_per_case':N,'primary_state_values':21,'primary_cases':8,'dt':dt,'seed':SEED,'note':'Common random histories across input states and thresholds; pointwise trials independent. Refinement uses independent ensembles.'}
assert abs(meas/pred-1)<.1
assert all(r['only1']+r['only2']+r['both']+r['none']==r['n'] for r in rows)
(D/'surface_basin_validation.json').write_text(json.dumps(validation,indent=2))
