from pathlib import Path
import json
import numpy as np
from PIL import Image,ImageDraw,ImageFont
P=Path(__file__).resolve().parent;r=json.loads((P/'results.json').read_text());font='/System/Library/Fonts/Supplemental/Arial.ttf'
im=Image.new('RGB',(1450,920),'white');d=ImageDraw.Draw(im)
def text(x,y,s,size=22):d.text((x,y),s,font=ImageFont.truetype(font,size),fill='#19364d')
text(50,25,'Transient complex response and a recording race',35)
text(50,80,'The complex response builds continuously; the recording rule remains a separate assumption.',22)
for k,idx in enumerate([0,12]):
 c=r['cases'][idx];x0=80+k*710;y0=200;W=570;H=300
 counts=np.array(c['clock_counts']).reshape(32,8).sum(axis=1);born=np.array(c['weights']).reshape(32,8).sum(axis=1)*c['clock_captures'];mx=max(counts.max(),born.max())*1.12
 text(x0,y0-62,'Identical local resonances' if k==0 else 'Unequal local resonances',26);text(x0,y0-26,f"{c['clock_captures']:,} records / 6,000 histories; count scale 0–{mx:.0f}",17)
 d.line((x0,y0,x0,y0+H,x0+W,y0+H),fill='#7b929f',width=2)
 for j,v in enumerate(counts):
  x=x0+(j+.5)*W/32;d.rectangle((x-5,y0+H-H*v/mx,x+5,y0+H),fill='#388ca6')
 d.line([(x0+(j+.5)*W/32,y0+H-H*v/mx) for j,v in enumerate(born)],fill='#c36c35',width=3)
 if 'independent_response_prediction' in c:
  pred=np.array(c['independent_response_prediction']).reshape(32,8).sum(axis=1)*c['clock_captures']
  d.line([(x0+(j+.5)*W/32,y0+H-H*v/mx) for j,v in enumerate(pred)],fill='#7052a6',width=3)
 text(x0,y0+H+12,'Site 0',17);text(x0+W-67,y0+H+12,'Site 255',17)
text(80,570,'Bars: counts in groups of 8 sites. Orange: incident Born weights. Purple: detector prediction.',21)
text(80,635,'Uniform linear response + squared-response Poisson rate preserves Born weights exactly.',23)
text(80,680,'That conclusion holds during the transient, not just after a steady response is reached.',23)
text(80,735,'A common response threshold still selects bright sites too often.',23)
text(80,790,'The response scaling is established within this model. Vacuum-generated recording',22)
text(80,825,'clocks, photon depletion and irreversible detection have not been derived.',22)
im.save(P/'summary.png')
rows=''
for c in r['cases']:
 fit=c['clock_fit'];rows+=f"<tr><td>{c['name'].replace('_',' ')}</td><td>{c['clock_captures']:,}</td><td>{fit['tv']:.4f}</td><td>{fit['holm_p']:.4f}</td><td>{c['threshold_captures']:,}</td></tr>"
body=f'''<h1>Testing the retrieved complex response mechanism</h1><p>1 October 2026</p>
<p><strong>The retrieved linear complex-response model preserves Born spatial weights exactly during its transient evolution when detector sites have identical responses and recording rates are proportional to squared response.</strong> This is distinct from the earlier sine-based phase-arrival race, which showed a finite bias. Fluctuation tests are compatible with Born counts at the tested sample size; unequal resonances produce predicted detector-response distortions.</p>
<p>Fourteen cases each used 6,000 independent incident histories: 84,000 total. Four paired numerical refinements reused 200 histories each. Both a Poisson recording race and a threshold comparator were evaluated on each history. No target site was drawn from the Born distribution beforehand.</p>
<img src="summary.png" alt="Uniform detector records follow Born weights; unequal resonances follow the detector-corrected prediction">
<h2>Connection to the Claude notes</h2>
<p>The read-only source is section 12 of Superposition_Equations (5).md in the user's Downloads. That model deliberately replaces Adler sine coupling with a damped complex oscillator. Linear input-amplitude coupling produces squared-amplitude response strengths, and an explicitly separate Poisson race generates outcomes. The present test implements the local-site generalization and follows the build-up from zero response instead of assuming a steady state.</p>
<h2>Local dynamics</h2><pre>dzⱼ/dt = −(γⱼ + iΔⱼ)zⱼ + ψⱼ exp(iθⱼ)[1 + b ξⱼ(t)]
zⱼ(0) = 0
Recording rate λⱼ(t) = 0.6 |zⱼ(t)|²
Alternative readout: first |zⱼ(t)|² to cross 1/256</pre>
<p>The shared photon profile ψ has normalized coherent amplitudes at 256 sites in the previous finite one-dimensional array. The incident drive is held constant during a dimensionless observation window 0≤t≤6. This is a linear-response representation of the incident excitation, not a normalized evolving photon–detector wavefunction. z is a response amplitude; |z|² is not automatically a probability or energy. Coupling/readout units and the finite drive duration are illustrative, with no fitted conversion to physical seconds.</p>
<p>Baseline damping is γ=1, detuning Δ=0, and θ=0. Slow and fast cases use γ=0.25 and 2 without changing the drive coefficient. A fixed-electronic-phase control assigns each site a random θ; this multiplies its whole response by a phase and cannot change its rate. It tests the distinction between arbitrary local phase and physically relevant interference between separate pathways.</p>
<p>The complex fluctuation is ξ(t)=u cos t+v sin t, with independent proper complex Gaussian spatial vectors u and v. Spatial covariance is independent-site, free-space Green-tensor, or ideal-mirror covariance, as in the earlier benchmark. Noise strengths b=0.25 and 1 are tested. Unlike the earlier normalized-coupling sweep, <strong>the mean drive strength is not rescaled to compensate for added fluctuations</strong>; both constructive/destructive interference and the increased mean squared response are retained.</p>
<p>Only spatial covariance is QED-inspired. The selected temporal spectrum and its modulation of a material transition remain assumptions, not a derived vacuum reservoir. Multiplication by ψ makes the drive vanish at exact optical zeros and when no photon drive is supplied. This suppresses vacuum-only excitation by construction; it does not model dark counts.</p>
<h2>Exact transient response</h2>
<p>For no fluctuations and lⱼ=γⱼ+iΔⱼ, the solution is zⱼ(t)=ψⱼ exp(iθⱼ)[1−exp(−lⱼt)]/lⱼ. For each harmonic component of the colored drive, the exact transfer function is [exp(iΩt)−exp(−lⱼt)]/(lⱼ+iΩ). These are evaluated directly, avoiding integration error in the complex oscillator itself. The time grid is used only for integrated recording rates and threshold crossings.</p>
<p>When the response parameters are identical, zⱼ(t)=ψⱼF(t) up to rate-irrelevant fixed phases. Therefore λⱼ(t)=wⱼh(t), where wⱼ=|ψⱼ|² and h(t)=0.6|F(t)|². The probability density for the first record is</p>
<pre>p(j,t) = wⱼ h(t) exp[−∫₀ᵗ h(s) ds].
P(j | a record by T) = wⱼ.</pre>
<p>This factorization is exact for any finite observation window with nonzero capture probability. It does not require steady state or a specially chosen phase-locking tolerance. The squared dependence follows from the linear response plus the assumed squared-response recording law. The exponential survival and individual event are supplied by the Poisson assumption, not by complex geometry alone.</p>
<h2>Parameter and generalization checks</h2>
<p>The same response and readout parameters are used for four input patterns: the original interference profile, a shifted relative optical phase, unequal beam amplitudes, and a single-beam envelope. We also vary damping, fixed electronic phases, absolute fluctuation strength and spatial field correlations. A detuning step sets Δ=0 on one half and 2 on the other; a separate fixed disordered detector uses Gaussian site detunings of standard deviation 0.8.</p>
<table><tr><th>Case</th><th>Poisson records</th><th>TV from incident Born weights</th><th>Adjusted p</th><th>Threshold records</th></tr>{rows}</table>
<p>All 12 clock cases without fixed site-response disorder are compatible with the incident Born reference after Holm correction of the 27 primary comparisons. For the uniform deterministic cases, agreement is also an exact analytic result. For fluctuating cases, this is a finite-sample compatibility statement, not proof of exact equality: correlated rate fluctuations and survival weighting can in principle alter the averaged distribution.</p>
<p>The two detuned cases disagree with the uncorrected incident reference. Their records agree with an independent calculation including the site-specific transient detector response: p≈0.643 for the step and p≈0.170 for fixed disorder. These distortions are expected detector response effects, not violations of quantum mechanics. A physically appropriate Born prediction for a nonuniform detector includes its measurement response.</p>
<p>Thirteen threshold cases produce 6,000 records each and all disagree strongly with the incident reference. The fast-response case never reaches its threshold and produces no threshold records. Thus linear complex response alone does not rescue the tested first-threshold rule. A record is still imposed at the first event; no microscopic circuit latch or cascade is derived here.</p>
<h2>Verification</h2>
<p>All 84,000 event histories were read back and matched against site counts and gate counts at times 1, 3 and 6. Integrated event probabilities plus no-event survival sum to one. Exact analytic capture fractions for uniform deterministic response agree with numerical hazard integration to below 4×10⁻⁸ in probability.</p>
<p>For four paired 200-history samples, time steps 0.05, 0.025 and 0.0125 use identical field coefficients and recording clocks. No clock winner changes in either halving. One threshold winner changes in the first halving of the strong independent-noise case, and none changes in the final halving. Maximum changes in mean local record weights decrease by approximately a factor of four, reaching below 5.6×10⁻⁸ in the final halving.</p>
<p>Two deterministic detuned predictions were independently integrated on 6,000- and 12,000-step grids. Their maximum local probability change is below 7.4×10⁻¹¹, and the main-grid predictions differ by below 1.9×10⁻⁷. Statistical comparisons use multinomial Monte Carlo reference distributions, accommodating sparse site counts.</p>
<h2>What has advanced</h2>
<p>The original complex-response idea is now tested in its own form, rather than replaced by an Adler phase-arrival law. Its amplitude-squared scaling and uniform-site spatial statistics are supported analytically, through time-dependent simulations and across changed incident profiles without retuning. Site response differences have a predictable effect.</p>
<p>The unresolved physical step is the origin of the recording process. Independent unit-exponential thresholds for integrated rates represent assumed Poisson clocks. Assigning those random numbers to vacuum fluctuations by name would not derive them. The linear response also omits photon depletion, a full single-excitation entangled state, material band structure and irreversible amplification. A detector-specific open-system derivation is needed to connect its damping, fluctuations and capture channel consistently.</p>
<p>The next test should therefore start from an explicit electronic transition and reservoir, derive the response and its event law together, and determine which stochastic assumptions are genuinely necessary. The earlier reversible single-excitation model supplies a useful comparison for excitation conservation; the present susceptibility model should not be mistaken for that complete quantum state.</p>
<p><a href="counts.csv">All raw site counts</a> · <a href="results.json">Parameters and results</a> · <a href="refinement.json">Paired numerical refinements</a> · <a href="validation.json">Verification record</a>. Event logs and reproducible scripts are saved alongside this report.</p>'''
(P/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Transient complex response test</title><style>body{max-width:1080px;margin:40px auto;padding:0 24px;font:18px/1.65 system-ui;color:#19364d}h1,h2{line-height:1.2}h2{margin-top:2em}img{width:100%}table{width:100%;border-collapse:collapse;font-size:14px}td,th{padding:9px;border-bottom:1px solid #cad6df;text-align:right}td:first-child{text-align:left}pre{white-space:pre-wrap;background:#edf4f8;padding:20px}a{color:#126a92}</style>'+body+'</html>')
