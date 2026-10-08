"""Single-photon to two-channel material excitation; no event sampling.
Run with numpy, scipy, matplotlib. Files saved beside this script.
"""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/born-mpl')
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json,csv
D=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
T=18.; beta=.8; Gamma=1.; ts=np.linspace(-T,0,9001)
# Markov single excitation sector. f has units 1/sqrt(time), e is dimensionless.
# Model matched input: normalized rising exponential stopped at t=0.
def rising(t): return np.exp(t/2)/np.sqrt(1-np.exp(-T))
def gaussian(t): return np.exp(-(t+3)**2/4)/(2*np.pi)**.25/gnorm
# Normalize the Gaussian on the same finite preparation interval.
gnorm=np.sqrt(np.trapezoid(np.exp(-(ts+3)**2/2)/np.sqrt(2*np.pi),ts))
def solve(pulse,delta=0.,b=beta):
    sol=solve_ivp(lambda t,e:-(Gamma/2+1j*delta)*e-1j*np.sqrt(b*Gamma)*pulse(t),[-T,0],[0j],t_eval=ts,rtol=2e-10,atol=2e-12)
    assert sol.success
    e=sol.y[0]; f=pulse(ts); fout=f-1j*np.sqrt(b*Gamma)*e
    incoming=cumulative_trapezoid(abs(f)**2,ts,initial=0)
    outgoing=cumulative_trapezoid(abs(fout)**2,ts,initial=0)
    loss=cumulative_trapezoid((1-b)*Gamma*abs(e)**2,ts,initial=0)
    return e,incoming,outgoing,loss
r,inc,out,loss=solve(rising); g,_,_,_=solve(gaussian)
# Exact finite-pulse solution at capture time.
def analytic_capture(delta,b=beta):
    return b/(1-np.exp(-T))*abs((1-np.exp(-(1+1j*delta)*T))/(1+1j*delta))**2
assert abs(abs(r[-1])**2-analytic_capture(0))<1e-8
balance=np.max(abs(inc-out-loss-abs(r)**2)); assert balance<3e-6
# Channel 2 has optional detuning; these states are orthogonal material modes.
e2,*_=solve(rising,delta=.7)
transfer=np.array([r[-1],e2[-1]])
c=np.array([np.sqrt(.3),np.sqrt(.7)*np.exp(1j*.8)])
v=transfer*c; rho=np.outer(v,v.conj()); eta=np.trace(rho).real
# Matched channels preserve any input qubit conditional on successful capture.
max_infidelity=0
for p in np.linspace(0,1,21):
    for phi in np.linspace(0,2*np.pi,25,endpoint=False):
        u=np.array([np.sqrt(p),np.sqrt(1-p)*np.exp(1j*phi)])
        z=r[-1]*u; normalized=z/np.linalg.norm(z)
        max_infidelity=max(max_infidelity,abs(1-abs(np.vdot(u,normalized))**2))
assert max_infidelity<1e-12
# Capture into stable storage is an ideal instantaneous coherent swap at t=0.
# Separate phenomenological dephasing + irreversible recording thereafter.
t=np.linspace(0,8,801); k=.6; gamma=.8
# Equal-superposition input, phi=.8, matched absorbers, unnormalized rho.
eta0=abs(r[-1])**2; phi=.8
pop=.5*eta0*np.exp(-k*t)
coh=.5*eta0*np.exp(-1j*phi)*np.exp(-(k+gamma)*t)
records=eta0*(1-np.exp(-k*t))
assert np.max(abs(2*pop+records-eta0))<1e-12
# Complementary, destructive readout after a storage delay, k=0 during delay.
# pplus is unconditional weight per incident photon, not conditional fidelity.
delta=np.linspace(0,2*np.pi,241)
fringes={delay:.5*eta0*(1+np.exp(-gamma*delay)*np.cos(delta-phi)) for delay in [0,1,3]}
# Detector orientation: real rotation acting on input polarization.
theta=np.linspace(0,np.pi,241); u=np.array([np.sqrt(.3),np.sqrt(.7)*np.exp(1j*.8)])
channel1=eta0*abs(np.cos(theta)*u[0]+np.sin(theta)*u[1])**2
channel2=eta0*abs(-np.sin(theta)*u[0]+np.cos(theta)*u[1])**2
assert max(abs(channel1+channel2-eta0))<1e-12
summary={'time_units':'Gamma inverse; dimensionless Gamma=1','beta':beta,'pulse_start':-T,'matched_capture':float(eta0),'gaussian_capture_at_0':float(abs(g[-1])**2),'gaussian_peak_excitation':float(max(abs(g)**2)),'flux_balance_max_error':float(balance),'matched_qubit_max_infidelity_525_states':float(max_infidelity),'detuned_channel_capture':float(abs(e2[-1])**2),'detuned_relative_transfer_phase':float(np.angle(e2[-1]/r[-1])),'unmatched_example_capture':float(eta),'record_rate':k,'coherence_decay_rate':gamma,'fringe_visibility_delays':{str(q):float(np.exp(-gamma*q)) for q in fringes}}
(D/'microscopic_summary.json').write_text(json.dumps(summary,indent=2))
with open(D/'microscopic_absorption.csv','w') as f:
 w=csv.writer(f);w.writerow(['t','rising_excited_weight','gaussian_excited_weight','input_cumulative','output_cumulative','loss_cumulative']); w.writerows(zip(ts,abs(r)**2,abs(g)**2,inc,out,loss))
with open(D/'microscopic_record.csv','w') as f:
 w=csv.writer(f);w.writerow(['t','remaining_excited_weight','coherence_magnitude','recorded_weight']);w.writerows(zip(t,2*pop,abs(coh),records))
fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
ax=axs[0,0];ax.plot(ts,abs(r)**2,label='Matched rising pulse');ax.plot(ts,abs(g)**2,label='Earlier Gaussian pulse');ax.plot(ts,inc,'k:',label='Cumulative incident photon weight');ax.set(xlim=(-10,0),xlabel='Time (1/Γ)',ylabel='Excitation / incident weight',title='A. Photon wavepacket → material excitation');ax.legend(fontsize=8)
ax=axs[0,1];ax.plot(theta*180/np.pi,channel1,label='Material channel 1');ax.plot(theta*180/np.pi,channel2,label='Material channel 2');ax.axhline(eta0,c='gray',ls=':',label='Total capture');ax.set(xlabel='Apparatus orientation (degrees)',ylabel='Excited-state weight',title='B. Apparatus orientation defines the mapping');ax.legend(fontsize=8)
ax=axs[1,0]
for delay,curve in fringes.items():ax.plot(delta,curve,label=f'Storage delay {delay}; visibility {np.exp(-gamma*delay):.2f}')
ax.axhline(eta0/2,c='gray',ls=':',label='Direct channel count (no phase analyzer)');ax.set(xlabel='Analyzer reference phase δ (radians)',ylabel='Outcome + weight per incident photon',title='C. Stored relative phase becomes measurable');ax.legend(fontsize=7)
ax=axs[1,1];ax.plot(t,2*pop,label='Unrecorded excitation');ax.plot(t,2*abs(coh),label='Twice coherence magnitude');ax.plot(t,records,label='Transferred to records');ax.set(xlabel='Time after ideal storage transfer (1/Γ)',ylabel='Weight / coherence',title='D. Dephasing and recording are separate processes');ax.legend(fontsize=8)
fig.suptitle('Minimal microscopic transfer model: absorption, coherence and readout\nβ = 0.8; matched orthogonal material channels; no sampled detection events',fontsize=14)
fig.savefig(D/'microscopic_transfer.png',dpi=170);fig.savefig(D/'microscopic_transfer.svg');plt.close(fig)
print(json.dumps(summary,indent=2))
