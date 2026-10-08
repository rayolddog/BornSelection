"""Numerical and physical consistency checks; run with Python + NumPy."""
import json
from dataclasses import replace
import numpy as np
from detector_model import *


def density_reference(p, psi):
    """Independent full density matrix with three sets of absorbing sinks."""
    n=p.n; dim=5*n
    g, delta=couplings(p)
    h=np.zeros((dim,dim),complex)
    for j in range(n):
        h[j,n+j]=g[j]; h[n+j,j]=g[j]; h[n+j,n+j]=delta[j]
    eig,v=np.linalg.eigh(h)
    u=(v*np.exp(-1j*eig*p.dt/2))@v.conj().T
    diag=np.ones(dim)
    diag[:n]=np.exp(-p.escape*p.dt/2)
    diag[n:2*n]=np.exp(-(p.capture+p.radiative)*p.dt/2)
    k0=np.diag(diag)
    initial=np.zeros(dim,complex); initial[:n]=psi
    rho=np.outer(initial,initial.conj())
    for step in range(round(p.gate/p.dt)):
        rho=u@rho@u.conj().T
        d=rho.diagonal().real.copy()
        updated=k0@rho@k0
        for j in range(n):
            lost=d[n+j]*(-np.expm1(-(p.capture+p.radiative)*p.dt))
            updated[2*n+j,2*n+j]+=lost*p.capture/(p.capture+p.radiative)
            updated[3*n+j,3*n+j]+=lost*p.radiative/(p.capture+p.radiative)
            updated[4*n+j,4*n+j]+=d[j]*(-np.expm1(-p.escape*p.dt))
        rho=u@updated@u.conj().T
    assert abs(np.trace(rho)-1)<1e-11
    return rho.diagonal().real


def run():
    rng=np.random.default_rng(55199); results={}
    for kind in ['single','double','shifted','unequal','incoherent']:
        psi,_=screen_state(kind,rng=rng,count=5)
        assert np.max(abs(np.sum(abs(psi)**2,axis=1)-1))<1e-14
    results['normalized_unitary_propagation']='PASS'
    error=0.
    for k in range(100):
        g=rng.normal()+1j*rng.normal(); d=rng.normal(); t=rng.uniform(0,3)
        v0=rng.normal(size=2)+1j*rng.normal(size=2)
        h=np.array([[0,np.conj(g)],[g,d]])
        ev,vec=np.linalg.eigh(h)
        reference=(vec*np.exp(-1j*ev*t))@vec.conj().T@v0
        a,b=unitary(v0[0],v0[1],g,d,t)
        error=max(error,float(np.max(abs(reference-[a,b]))))
    assert error<1e-12
    results['local_unitary_max_error']=error
    p=Parameters(disorder=True)
    e=efficiency(p)
    err=float(np.max(abs(sum(e.values())-1)))
    assert err<1e-11
    results['probability_conservation_max_error']=err
    # Coupling phase is removable when each independent electron starts unexcited.
    a=rng.normal(size=64)+1j*rng.normal(size=64); a/=np.linalg.norm(a)
    phase=np.exp(1j*rng.uniform(-np.pi,np.pi,64)); g,d=couplings(p)
    aa,bb=unitary(a,np.zeros(64),g,d,.3)
    ap,bp=unitary(a,np.zeros(64),g*phase,d,.3)
    assert np.max(abs(abs(bb)**2-abs(bp)**2))<1e-14
    ag,bg=unitary(a*np.exp(1.23j),np.zeros(64),g,d,.3)
    assert np.max(abs(abs(ag)**2-abs(aa)**2))<1e-14
    results['local_coupling_and_global_optical_phase_invariance']='PASS'
    a,b=unitary(np.zeros(64,complex),np.zeros(64,complex),g,d,.3)
    cp,rp,ep=damping_probabilities(a,b,p)
    assert (cp+rp+ep).sum()==0
    results['vacuum_has_no_absorption']='PASS'
    # Compare discrete splitting with an independent continuous calculation.
    continuous=continuous_efficiency(p)
    errs=[]
    for dt in [.1,.05,.025]:
        errs.append(float(np.max(abs(efficiency(replace(p,dt=dt))['capture']-continuous))))
    assert errs[1]<errs[0]/3.5 and errs[2]<errs[1]/3.5
    results['continuous_capture_max_errors_dt_010_005_0025']=errs
    # Full density operator comparison: n=16 still contains the specified apertures.
    small=Parameters(n=16,gate=5,dt=.1,disorder=True)
    psi=screen_state('unequal',16)[0][0]
    rho=density_reference(small,psi)
    response=efficiency(small)
    assert np.max(abs(rho[32:48]-abs(psi)**2*response['capture']))<1e-11
    event=quantum_trials('unequal',small,count=10000,seed=2026121)
    observed=np.concatenate([histogram(event,16),
        np.bincount(event['site'][event['outcome']==2],minlength=16),
        np.bincount(event['site'][event['outcome']==3],minlength=16),
        [np.count_nonzero(event['outcome']==0)]])
    predicted=np.concatenate([rho[32:80],[rho[:32].sum()]])
    fit=goodness(observed,predicted,seed=2231)
    assert fit['p']>.001
    results['trajectory_vs_full_density_matrix']=fit
    rates=np.array([.2,.8])
    first=np.argmin(rng.exponential(size=(100000,2))/rates,axis=1)
    frequency=float(np.mean(first==0))
    assert abs(frequency-.2)<.005
    results['two_channel_exponential_race_expected_020_observed']=frequency
    no=threshold_trials(vacuum=True,count=100)
    assert no['outcome'].sum()==0
    results['zero_signal_zero_noise_threshold']='PASS'
    return results

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
