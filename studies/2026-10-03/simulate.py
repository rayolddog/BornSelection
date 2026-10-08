"""Reproducible two-outcome basin experiments. Python + numpy/scipy/matplotlib.
Run: MPLCONFIGDIR=/tmp/born-mpl python3 simulate.py
All outputs are written beside this file. No outcome is sampled from Bernoulli(p).
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/born-mpl')
import numpy as np
from scipy.special import ndtr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv, json
OUT=Path(__file__).resolve().parent
SEED=20261003
rng=np.random.default_rng(SEED)
P=np.linspace(0,1,41); N=8000
plt.rcParams.update({'figure.dpi':150,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})

def walk(p,n,noise='gaussian',h=.07,rho=0.,gain=0.,amplitude=False,seed=None,shared=False):
    r=np.random.default_rng(seed) if seed is not None else rng
    p=np.asarray(p); initial=np.sqrt(p)/(np.sqrt(p)+np.sqrt(1-p)) if amplitude else p
    x=np.repeat(initial,n).copy(); memory=r.normal(size=x.size)
    steps=np.zeros(x.size,dtype=int)
    for t in range(20000):
        ids=np.flatnonzero((x>0)&(x<1))
        if not len(ids): break
        if shared:
            eta=r.normal(size=n)[ids % n]
        elif noise=='laplace': eta=r.laplace(size=len(ids))/np.sqrt(2)
        elif noise=='binary': eta=r.choice([-1.,1.],size=len(ids))
        else: eta=r.normal(size=len(ids))
        if rho:
            memory[ids]=rho*memory[ids]+np.sqrt(1-rho*rho)*eta
            eta=memory[ids]
        old=x[ids]; proposal=h*eta+gain*(2*old-1)*h*h
        # Symmetric distance-to-boundary regulator: preserves zero conditional
        # drift for independent symmetric noise, unlike ordinary clipping.
        magnitude=np.minimum(np.abs(proposal),np.minimum(old,1-old))
        new=old+np.sign(proposal)*magnitude
        # Snap only arithmetic roundoff, never a finite epsilon interval.
        hit=magnitude==np.minimum(old,1-old)
        new=np.where(hit & (proposal>0) & (old>=.5),1.,new)
        new=np.where(hit & (proposal<0) & (old<=.5),0.,new)
        x[ids]=new; steps[ids]+=1
    x=x.reshape(len(p),n)
    return (x==1).sum(axis=1),((x>0)&(x<1)).sum(axis=1),steps.reshape(len(p),n).mean(axis=1),x

def wilson(k,n):
    z=1.95996398454; f=k/n; d=1+z*z/n
    a=(f+z*z/(2*n))/d; b=z*np.sqrt(f*(1-f)/n+z*z/(4*n*n))/d
    return a-b,a+b

rows=[]; results={}
def add(name,k,u,steps,n=N,p=P):
    lo,hi=wilson(k,n); f=k/n
    for j,v in enumerate(p): rows.append(dict(model=name,p=float(v),n=n,upper=int(k[j]),unresolved=int(u[j]),frequency=float(f[j]),ci_low=float(lo[j]),ci_high=float(hi[j]),mean_steps=float(steps[j])))
    results[name]=(np.array(p),f,lo,hi)
    print(name,'RMSE',round(float(np.sqrt(np.mean((f-p)**2))),5),'max',round(float(np.max(abs(f-p))),5),'unresolved',int(sum(u)),flush=True)

models=[('Unbiased / Gaussian',dict()),('Unbiased / Laplace',dict(noise='laplace')),('Unbiased / binary',dict(noise='binary')),('Memory / Gaussian rho=.9',dict(rho=.9)),('Amplifier / Gaussian gain=4',dict(gain=4)),('Amplitude coordinate / Gaussian',dict(amplitude=True))]
for name,kw in models:
    k,u,s,_=walk(P,N,**kw); add(name,k,u,s)

# Free-space dipole-inspired symmetrized spectrum with explicit UV cutoff.
# Independent Gaussian quadratures are sampled before detector evolution.
# eta(t)=sum sqrt(w_j) [A_j cos(omega_j t)+B_j sin(omega_j t)].
# A leaky integrator z'= -lambda*z + eta(t), z(0)=0, is solved exactly.
w=np.linspace(.05,10,96); dw=w[1]-w[0]; cutoff=2.; lam=1.; T=3.
spec=w**3*np.exp(-w/cutoff); weights=spec*dw; weights/=weights.sum()
H=(np.exp(1j*w*T)-np.exp(-lam*T))/(lam+1j*w)
ca=np.sqrt(weights)*H.real; cb=np.sqrt(weights)*H.imag
sd=np.sqrt(np.sum(ca*ca+cb*cb))
# Fixed sigma, not fitted to Born frequencies. Sign selects deterministic basin
# of dy/dt=y-y^3 after the preparation pulse.
sigma=.6
ks=[]
for p in P:
    a=rng.normal(size=(N,len(w))); b=rng.normal(size=(N,len(w)))
    kick=(a@ca+b@cb)/sd
    ks.append(np.count_nonzero((2*p-1)+sigma*kick>0))
add('Vacuum-spectrum kick / bistable',np.array(ks),np.zeros(len(P),int),np.ones(len(P)))

# Resolution and parameter checks, using independent ensembles.
checks=[]
for h in [.14,.07,.035]:
    for label,kw in [('Unbiased',{}),('Amplifier',{'gain':4})]:
        ps=np.array([.1,.25,.5,.75,.9]); k,u,s,_=walk(ps,12000,h=h,**kw)
        for j,p in enumerate(ps): checks.append(dict(model=label,h=h,p=p,n=12000,frequency=k[j]/12000,unresolved=int(u[j])))
for rho in [0.,.3,.6,.9]:
    ps=np.array([.1,.25,.5,.75,.9]); k,u,s,_=walk(ps,12000,rho=rho)
    for j,p in enumerate(ps): checks.append(dict(model='Memory',h=.07,rho=rho,p=p,n=12000,frequency=k[j]/12000,unresolved=int(u[j])))
for filename,data in [('basin_measures.csv',rows),('sensitivity.csv',checks)]:
    keys=list(dict.fromkeys(k for row in data for k in row))
    with open(OUT/filename,'w') as f:
        writer=csv.DictWriter(f,fieldnames=keys); writer.writeheader(); writer.writerows(data)
summary={name:dict(rmse=float(np.sqrt(np.mean((v[1]-v[0])**2))),max_error=float(np.max(abs(v[1]-v[0])))) for name,v in results.items()}
(OUT/'summary.json').write_text(json.dumps({'seed':SEED,'n_per_state':N,'states':len(P),'models':summary},indent=2))
fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
for ax,indices,title in [(axs[0,0],[0,1,2],'Independent symmetric noise'),(axs[0,1],[3,4],'Memory and amplification'),(axs[1,0],[5,6],'Other detector couplings')]:
    ax.plot(P,P,'k--',lw=1,label='Born target')
    for i in indices:
        name=list(results)[i]; p,f,lo,hi=results[name]
        ax.errorbar(p,f,yerr=[f-lo,hi-f],fmt='.-',lw=1,ms=3,capsize=1,label=name)
    if ax is axs[1,0]:
        ax.plot(P,ndtr((2*P-1)/sigma),color='gray',ls=':',label='Gaussian threshold prediction')
    ax.set(xlabel='Initial squared amplitude |c₁|²',ylabel='Fraction reaching outcome 1',title=title,xlim=(0,1),ylim=(-.03,1.03));ax.legend(fontsize=7)
ax=axs[1,1]
for name in list(results)[:3]:
    p,f,lo,hi=results[name];ax.errorbar(p,f-p,yerr=[f-lo,hi-f],fmt='.',capsize=2,label=name)
ax.axhline(0,c='k',ls='--');ax.set(title='Benchmark residuals (95% pointwise intervals)',xlabel='Initial squared amplitude |c₁|²',ylabel='Measured frequency − Born target');ax.legend(fontsize=7)
fig.suptitle('Two-outcome selection: detector dynamics decides the basin measure\n8,000 independent trials per state; 41 states per model',fontsize=14)
fig.savefig(OUT/'frequencies.png');fig.savefig(OUT/'frequencies.svg');plt.close(fig)
# A single noise tape per column, reused across initial states; no sorting.
ps=np.linspace(0,1,101); _,_,_,basins=walk(ps,160,seed=4001,shared=True)
fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
axs[0].imshow(basins,origin='lower',aspect='auto',extent=[0,160,0,1],cmap='coolwarm',vmin=0,vmax=1,interpolation='nearest');axs[0].set(xlabel='Microscopic noise tape (unsorted)',ylabel='Initial |c₁|²',title='Basin slice: unbiased Gaussian competition')
axs[1].plot(ps,basins.mean(axis=1),label='160 displayed tapes');axs[1].plot(ps,ps,'k--',label='Born target');axs[1].set(xlabel='Initial |c₁|²',ylabel='Fraction in outcome-1 basin',title='Slice measure (illustration only)');axs[1].legend()
fig.savefig(OUT/'basins.png');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for model in ['Unbiased','Amplifier']:
    for h in [.14,.07,.035]:
        z=[q for q in checks if q['model']==model and q['h']==h];axs[0].plot([q['p'] for q in z],[q['frequency'] for q in z],'.-',label=f'{model}, h={h}')
axs[0].plot(P,P,'k--');axs[0].set(title='Step-size check',xlabel='Initial |c₁|²',ylabel='Outcome-1 frequency');axs[0].legend(fontsize=7)
for rho in [0,.3,.6,.9]:
    z=[q for q in checks if q['model']=='Memory' and q['rho']==rho];axs[1].plot([q['p'] for q in z],[q['frequency'] for q in z],'.-',label=f'rho={rho}')
axs[1].plot(P,P,'k--');axs[1].set(title='Dependence on temporal memory',xlabel='Initial |c₁|²',ylabel='Outcome-1 frequency');axs[1].legend()
fig.savefig(OUT/'sensitivity.png');plt.close(fig)
print('DONE',flush=True)
