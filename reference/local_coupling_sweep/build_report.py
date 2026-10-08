from pathlib import Path
import json,html
import numpy as np
from PIL import Image,ImageDraw,ImageFont
P=Path(__file__).resolve().parent
r=json.loads((P/'results.json').read_text());ref=json.loads((P/'refinement.json').read_text());val=json.loads((P/'validation.json').read_text());w=np.array(r['weights'])
font='/System/Library/Fonts/Supplemental/Arial.ttf'
im=Image.new('RGB',(1450,1080),'white');d=ImageDraw.Draw(im)
def text(x,y,s,size=22,col='#18344b'):d.text((x,y),s,font=ImageFont.truetype(font,size),fill=col)
text(55,28,'Local complex coupling at 256 microscopic sites',36)
text(55,81,'One-photon trials across several interference fringes; counts shown in groups of 8 sites.',22)
c=r['results'][7]
for k,ri in enumerate([0,3]):
 s=c['results'][ri];x0=85+k*700;y0=200;W=570;H=300
 observed=np.array(s['counts']).reshape(32,8).sum(axis=1);expected=s['recorded']*w.reshape(32,8).sum(axis=1);mx=max(observed.max(),expected.max())*1.12
 text(x0,y0-68,'First local threshold crossing' if k==0 else 'Competing absorption clocks',25)
 text(x0,y0-30,f'Raw counts; top of scale {mx:.0f}',17)
 d.line((x0,y0,x0,y0+H,x0+W,y0+H),fill='#8094a1',width=2)
 for j,v in enumerate(observed):
  x=x0+(j+.5)*W/32;d.rectangle((x-6,y0+H-H*v/mx,x+6,y0+H),fill='#b85642' if k==0 else '#3485a4')
 d.line([(x0+(j+.5)*W/32,y0+H-H*v/mx) for j,v in enumerate(expected)],fill='#262d58',width=4)
 text(x0,y0+H+12,'Site 0',17);text(x0+W-68,y0+H+12,'Site 255',17)
 text(x0,y0+H+50,f"{s['recorded']:,} records / 4,000 trials",22)
 text(x0,y0+H+88,'Line: Born reference scaled to recorded count',18)
text(85,660,'Parameter sweep findings',28)
text(85,706,'32 of 32 well-populated threshold settings disagree with Born weights.',24)
text(85,750,'12 of 12 clock comparisons show no significant disagreement after correction.',24)
text(85,805,'Thresholds: only the local fluctuating amplitudes determine the winner.',21)
text(85,842,'Clocks: additional exponential randomness and squared-amplitude rates are assumed.',21)
text(85,910,'Finer time steps preserve this contrast; threshold capture efficiencies are not converged.',21)
text(85,952,'Spatial vacuum correlations are QED-inspired. The full time-dependent field and selection',21)
text(85,984,'law are phenomenological; no QED-derived Born-selection mechanism is established.',21)
im.save(P/'summary.png')
rows=''
for c in r['results']:
 cells=[]
 for s in c['results']:
  if s['fit']:
   f=s['fit'];cells.append(f"{s['recorded']:,}<br>TV {f['tv']:.3f}<br>adjusted p {f['holm_p']:.4f}")
  else:cells.append(f"{s['recorded']:,}<br>insufficient records")
 rows+=f"<tr><td>{c['name'].replace('_',' ')}<br><small>a={c['electronic_strength']:g}, b={c['vacuum_strength']:g}, τ={c['correlation_time']:g}</small></td>"+''.join('<td>'+v+'</td>' for v in cells)+'</tr>'
refrows=''
for c in ref:
 s=c['results'][2];cl=c['results'][3]
 refrows+=f"<tr><td>{c['name']}</td><td>{s['recorded']:,}</td><td>{s['fit']['tv']:.3f}</td><td>{cl['recorded']:,}</td><td>{cl['fit']['p']:.3f}</td></tr>"
body=f'''<h1>Testing local complex interactions as a detector selection mechanism</h1>
<p>1 October 2026 · Autonomous continuation authorized before the user's flight</p>
<p><strong>Resolving the interaction at individual sites and adding local complex fluctuations did spread events beyond the brightest site, but none of the well-populated threshold settings reproduced Born statistics.</strong> A comparison using local absorption clocks was compatible with those statistics. That comparison assumes a squared-amplitude rate law and extra randomness, so it does not establish the proposed vacuum-only selection mechanism.</p>
<p>There were 48,000 independent incident histories in the main sweep, plus 24,000 in finer-time-step checks. Each history was assessed under four alternative recording rules; those are counterfactual comparisons, not four photons. No Born-weighted site was drawn in advance.</p>
<img src="summary.png" alt="Local threshold counts remain excessively concentrated at bright fringes; absorption clocks follow the reference pattern more closely">
<h2>What was calculated at each site</h2>
<p>A finite, one-dimensional array contains 256 individually evaluated dipole sites, spaced 0.01 wavelengths apart. An illustrative scale is a 50 nm photon wavelength, a 0.5 nm site spacing and a 127.5 nm array width. This is an idealized microscopic test strip, not a model of every atom in a commercial detector, not a calibrated semiconductor, and not the 2014 SPAD experiment. It omits electronic band structure, intersite electronic transport and material-specific energy levels.</p>
<p>Two coherent beams crossing at ±60° to the detector normal produce several fringes. Their complex projected fields are added before the modulus is taken. A Gaussian envelope gives a finite illuminated region. This directly specifies the incident interference field; a slit propagation calculation is not part of this test.</p>
<pre>ψⱼ ∝ exp[−xⱼ²/(4σ²)] {{exp(ik sinθ xⱼ) + exp(iφ) exp(−ik sinθ xⱼ)}}
σ = 0.85λ, θ = 60°, φ = 0.4;  Σⱼ |ψⱼ|² = 1
Gⱼ(t) = [1 + a exp(iϑⱼ(t)) + b ξⱼ(t)] / √C
Mⱼ(t) = ψⱼ Gⱼ(t)
C = 1 + a² + b² Kⱼⱼ</pre>
<p>The electronic phase ϑ is the relative phase of an assumed second local transition pathway, not the arbitrary overall phase of an electron. It undergoes phase diffusion, independently at each site, with uniform initial phase and correlation time τ. The complex Gaussian field ξ has the chosen spatial covariance K and Ornstein–Uhlenbeck time correlation exp(−|Δt|/τ). Both fluctuate during each photon history.</p>
<p>The scalar expressions are projected dipole couplings for a fixed polarization, rather than complete three-dimensional electronic wavefunctions. The normal-dipole mirror case assumes illumination polarized to provide a nonzero normal transition projection. Its projected incident profile is normalized to the same ψ.</p>
<p>This is <strong>an assumed modulation of a photon-assisted transition</strong>. It does not represent the full photon–atom–vacuum Hamiltonian. In particular, the common local photon phase factors out of |Mⱼ|² after its two optical paths have interfered. Merely retaining that phase in code is not proof that the model contains every physically relevant phase relation. Any additional coherent material pathways would need to be specified.</p>
<h2>How vacuum-related structure enters</h2>
<p>Four spatial choices were examined: independent sites as a control; the free-space projected electromagnetic Green tensor; and ideal perfectly conducting mirror geometries with dipoles parallel or normal to a surface at height 0.05λ. The Green tensor links field correlations, material transition rates and Casimir-related physics in <a href="https://arxiv.org/abs/0902.3586">macroscopic QED</a>.</p>
<p>At this height, the model's single-frequency local variance relative to free space is approximately 0.0773 for the parallel case and 1.9611 for the normal case. These are suppression and enhancement of the selected correlation, not a universal increase in vacuum energy. The relevant separation is the dipole site's distance from the surface, not a virtual photon's wavelength assigned from an incoming photon's distance.</p>
<p><strong>Only the spatial covariance has that QED basis.</strong> The rotating-frame Ornstein–Uhlenbeck time law, noise strength and coupling modulation are phenomenological. They are not the full ordered QED vacuum correlations and do not include the required material response, dissipation and frequency dependence. The mirror is idealized; the illustrative ultraviolet scale must not be treated as realistic optical constants for a metal.</p>
<p>The common normalization C makes the ensemble mean |Gⱼ|² equal to one at every site. This allows shape comparisons at similar mean drive. It deliberately calibrates away the change in mean total coupling and cannot predict absolute Casimir-induced changes in detection efficiency. Nearby sites share strongly correlated vacuum-inspired noise: the spatial covariance is effectively rank 13 across this small array, rather than 256 independent vacuum draws.</p>
<h2>The competing selection rules</h2>
<p>The photon pulse has normalized intensity envelope f(t) = (2/T)sin²(πt/T), 0≤t≤T, with T=10 in arbitrary slow-envelope time units. These units are not optical cycles and have no fitted conversion to seconds.</p>
<p><strong>Threshold candidates:</strong> define sⱼ(t)=|Mⱼ(t)|²sin²(πt/T). The first site crossing a common threshold η/256 wins; η takes values 0.5, 2 and 8. Crossings are interpolated between sampled times. The threshold concerns an assumed local coupling-strength excursion; it is neither deposited energy nor a derived microscopic detector record variable. It has no additional outcome draw. The single-record global stop is imposed, not derived by simulating photon depletion and detector entanglement.</p>
<p><strong>Absorption-clock comparison:</strong> set the instantaneous rate λⱼ(t)=3 f(t)|Mⱼ(t)|². Each site receives an independent unit-mean exponential random threshold for integrated rate; the first clock to fire records the event. This is a standard competing-hazard construction. It does not choose a target site beforehand, but it does insert a rate proportional to squared amplitude and extra stochastic clocks. Calling their random numbers vacuum fluctuations would not derive their physical origin.</p>
<p>Without fluctuations, the clock model gives P(record)=1−e⁻³ and P(site j | record)=|ψⱼ|² exactly. With fluctuating rates, the general prediction is E[∫λⱼ(t) exp(−∫₀ᵗΣₖλₖ(s)ds)dt]. Equal mean rates alone do not guarantee exact Born statistics: temporal correlations and survival weighting matter.</p>
<h2>Main parameter sweep and raw counts</h2>
<p>Each row uses 4,000 independent incident histories. a and b are the relative electronic and vacuum-inspired strengths, and τ is their correlation time. Results below use the final gate T=10 and the initial sampling interval Δt=0.2. Counts at gates 2 and 5, all 256 site counts, and each event's selected site and time are saved. No sparse-bin asymptotic test is used.</p>
<table><tr><th>Fluctuation case</th><th>Threshold 0.5</th><th>Threshold 2</th><th>Threshold 8</th><th>Absorption clocks</th></tr>{rows}</table>
<p>TV denotes total-variation distance from the Born site distribution: zero is identical and one is completely separated. Each formal p-value uses 10,000 multinomial reference simulations. Holm correction covers all 44 final-gate tests with at least 1,000 records. All 32 eligible threshold cases reject the Born reference; none of the 12 clock cases rejects after correction. A p-value at the simulation floor is not an exact estimate of an extremely small probability.</p>
<p>The no-fluctuation threshold model records the brightest site whenever its threshold is reachable. Electronic relative-phase fluctuations broaden this, and stronger local field fluctuations broaden it further, but the tested thresholds still concentrate records too strongly in bright fringes. The most favorable main threshold TV is about 0.461, compared with about 0.085–0.101 for the finite-sample clock comparisons. Compatibility at 4,000 trials does not establish exact equality of the noisy clock laws with Born probabilities.</p>
<p>Two high-threshold cases produced no records; the weak-field and parallel-mirror cases produced 48 and 84. These are insufficient for distribution claims and remain raw exploratory results. They are not counted as successful fits. All decisions and counts are retained; parameter settings were not dropped to improve agreement.</p>
<h2>Numerical checks and a remaining limitation</h2>
<p>All 48,000 primary histories were reconstructed from their saved event logs and matched the per-site and finite-gate counts. The no-fluctuation clock capture fractions at all three gates agree with their exact formula within 1.08 binomial standard deviations. The incident state is normalized. A global phase change leaves transition weights unchanged to 5.6×10⁻¹⁷. An exact zero of the incident field remains dark by construction; this is a property of the multiplicative ansatz, not an independent prediction about real background counts.</p>
<p>Independent 50,000-draw covariance checks reproduce the selected free-space and mirror matrices with maximum normalized errors below 0.014. Tiny negative covariance eigenvalues are numerical roundoff; retained-mode reconstruction errors are below 1.6×10⁻¹⁰.</p>
<p>Three representative cases were rerun at Δt=0.05 and 0.025 with 4,000 new histories each. These are independent statistical refinements, not matched-path error estimates. The threshold distributions remain strongly biased; all six clock comparisons have unadjusted p&gt;0.17. However, threshold-8 capture counts increase with finer sampling because more short excursions are resolved.</p>
<table><tr><th>Refinement</th><th>Threshold-8 records</th><th>Threshold-8 TV</th><th>Clock records</th><th>Clock raw p</th></tr>{refrows}</table>
<p><strong>Threshold efficiencies and precise first-crossing times are not demonstrated to be converged.</strong> The reported negative spatial result persists across these refinements, but the coarse sweep should not be read as a quantitatively converged detection-efficiency prediction. A continuous first-passage calculation or further refinement would be needed for that purpose.</p>
<h2>What the test establishes and what it leaves open</h2>
<p>The user's intuition is partly borne out: local complex fluctuations let weaker sites win sometimes. They do not automatically produce the required frequencies. The particular first-excursion threshold family tested here remains biased despite resolving every site in the finite array, retaining coherent incident phases, varying local phase dynamics and introducing spatially correlated surface-modified noise.</p>
<p>The clock comparison shows how a successful statistical pattern can arise without drawing a Born-weighted site beforehand, but its squared-amplitude rate and independent exponential thresholds already supply the missing event law. This is a benchmark, not a vacuum-driven derivation. The earlier carrier model can amplify a resulting local capture, but amplification cannot retrospectively establish its selection law.</p>
<p>No electron band structure, full entangled photon–material evolution, radiative backreaction, true many-body QED vacuum, irreversible circuit record, or real experimental fit is claimed. The bounded parameter sweep rejects these tested threshold rules; it does not reject every possible vacuum-related selection mechanism.</p>
<p><strong>The most useful next step is to derive the local event law from a specified electronic transition and reservoir coupling.</strong> It should account for dissipation and fluctuations consistently, conserve the one-photon excitation appropriately, and state whether exponential clocks follow from the derivation or remain an assumption. Further threshold tuning alone would not resolve that physical question.</p>
<h2>Saved materials</h2>
<p><a href="counts.csv">All primary raw site counts</a> · <a href="results.json">Parameters and statistics</a> · <a href="refinement.json">Time-step refinements</a> · <a href="validation.json">Verification record</a>. Eighteen compressed logs preserve the 72,000 independent histories and four alternative selections per history. Scripts are saved alongside the report; their two helper modules are included in the full research recovery archive.</p>'''
(P/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Testing local complex detector coupling</title><style>body{max-width:1080px;margin:40px auto;padding:0 24px;font:18px/1.65 system-ui;color:#18344b}h1,h2{line-height:1.2}h2{margin-top:2em}img{width:100%}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border-bottom:1px solid #cad6de;padding:10px;text-align:right}td:first-child{text-align:left}pre{white-space:pre-wrap;background:#eff4f8;padding:20px}a{color:#126a92}small{color:#526574}</style>'+body+'</html>')
