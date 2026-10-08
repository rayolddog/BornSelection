# Microscopic Complex Wave Response and the Selection of Individual Quantum Records

A phenomenological model, numerical evidence, and a proposed correlation experiment

**Working draft 0.1 — 1 October 2026**

**Provisional repository contributor order:** Codex (OpenAI); Claude (Anthropic); John Bramble. Claude is credited for prior model development and has not reviewed or approved this manuscript. The contribution statement explains the status of this byline.

## Abstract

We investigate a microscopic account of quantum detection in which complex quantum waves represent a physically existing process and local fluctuations contribute to the selection of an individual recording site. We call this individual outcome formation *Born selection*, distinguishing it from the ensemble Born distribution. A damped complex response driven by the incident amplitude, combined with a first-event recording law proportional to squared response, gives exact Born site probabilities for identical deterministic detector responses at every observation gate. Existing numerical studies support this result and identify failures of common-threshold and phase-arrival alternatives. These calculations do not derive the recording law or establish a quantum-vacuum origin of its randomness. We propose a two-channel experiment comparing site probabilities and recording times as the spatial correlations and memory of local fluctuations are varied. A newly specified quasistatic complex-gain race gives quantitative departures from a matched quantum absorption model, while a common-fluctuation limit retains Born weights. This conditional prediction makes one version of the mechanism falsifiable; it is not a universal consequence of complex-wave realism. We give the assumptions, calibration requirements, competing explanations, and remaining theoretical work needed for an experimental assessment.

## Purpose and physical postulate

The motivating postulate is that *complex quantum waves represent a physically existing process, rather than merely a mathematical tool for predicting experimental outcomes*. Here "physical" does not mean real-valued, and it does not make the real and imaginary components separately observable. An overall phase is not an observable in this model. Relative phases enter local transition amplitudes and therefore can affect detector response.

The additional dynamical hypothesis is that fluctuations acting during the interaction between an incident excitation and detector matter help determine where the first record forms. We study the sequence of coherent propagation, local coupling, reversible excitation, stochastic commitment, and amplification. The aim is a quantitatively specified microscopic mechanism that produces one record per incident excitation and the observed distribution over many trials. A phenomenological exploration can precede a derivation of every parameter from material constants. However, parameters and assumptions must be explicit if successful fitting is to become a prediction.

This manuscript reports a research model, not a completed solution of the measurement problem. In particular, a physical interpretation of amplitudes alone supplies no new experimental prediction. The event law and fluctuation dynamics must add testable content. Nor does a one-particle detector model establish a consistent ontology for general entangled many-particle states: their configuration-space structure, relativistic extension, and no-signalling requirements remain outside the demonstrated scope.

## Three stages of the proposed measurement process

The proposed mechanism organizes measurement into three stages:

1.  **Unitary propagation in free flight.** The incident photon wave packet evolves unitarily during propagation when coupling to uncontrolled environmental degrees of freedom can be neglected. Preparation and slit geometry determine the spatial amplitudes and relative phases arriving at the detector.

2.  **Formation of a local interference response at the detector surface.** The incident photon couples to the bound-electron wave packets at the detector's active surface. In the proposed mechanism, quantum-vacuum fluctuations contribute to these local complex interaction amplitudes. The resulting pattern of excitation depends on the incident amplitude and phase at each atom, the electronic transition response, and the local fluctuation field. Here "interference" refers to contributions to a joint transition amplitude, rather than direct addition of photon and electron wavefunctions. The vacuum origin of the phenomenological fluctuation field remains a hypothesis.

3.  **Evolution toward a selected site and amplification.** The complex excitations at the surface evolve, and a postulated stochastic commitment process selects an individual site at which amplification begins. Subsequent interactions in the detector convert that initial event into a persistent macroscopic record. Repeating this process with one incident excitation at a time is intended to produce the observed Born distribution of recorded sites.

These stages distinguish propagation, local response, and record formation; they need not be separated by sharp time boundaries. "Surface" denotes the active absorption region, which may extend into the detector material. The present response model represents the second stage, and its first-event law models commitment in the third. Deriving that law from the same microscopic complex dynamics, and connecting it to the amplification cascade, remain the central tasks.

## Relation to quantum measurement theory

Decoherence describes the suppression of interference in a subsystem through correlations with an environment. Its relation to definite outcomes depends on the larger interpretation or measurement framework adopted [\[1\]](https://arxiv.org/abs/quant-ph/0312059). It is therefore insufficient to compare the proposed mechanism with an unspecified "decoherence prediction." Our experimental comparator is a calibrated open quantum model, together with a specified measurement instrument that predicts site labels, recording times, and missing records. Quantum trajectory methods already describe stochastic individual detection records within standard quantum theory [\[2\]](https://arxiv.org/abs/0805.4002).

The proposed contribution is the investigation of a physical source and intermediate dynamics for stochastic commitment. It does not follow that established quantum theory lacks methods for calculating individual-record statistics. The distinction is between reproducing those statistics with an explicit model and independently deriving why that model's event law should hold.

Macroscopic QED relates material boundaries to field correlations, atomic relaxation, and level shifts [\[3\]](https://arxiv.org/abs/0902.3586). The artificial-atom mirror experiment of Hoi et al. demonstrates that a boundary can strongly modify spontaneous emission [\[4\]](https://arxiv.org/abs/1410.8840). This motivates boundary control, but does not establish that the vacuum selects a detector site, or that its correlation time is adjustable over the range required below. Reversibility experiments also motivate examining intermediate stages before macroscopic recording [\[5\]](https://arxiv.org/abs/0806.3547). Neither result supplies the missing event-selection derivation.

## Local complex response

Discretize the illuminated detector into $N$ effective channels and write the normalized incident state as 

$$
|\psi\rangle=\sum_{j=1}^{N}\psi_j|j\rangle,\qquad
 w_j=|\psi_j|^2,\qquad \sum_j w_j=1.
\tag{1}
$$

 A channel may represent a localized excitation mode or a group of atoms with a specified transition response. The exploratory numerical lattice is not a resolved semiconductor calculation.

For an electric-dipole transition, a local coupling is constructed from a contraction such as $M_j=\mathbf d_{eg,j}\cdot\boldsymbol{\mathcal E}_\psi(\mathbf r_j)$, with conjugations fixed by the chosen state and field conventions. More generally, the relevant quantity is a matrix element $\langle f|H_{\rm int}|i\rangle$. This is not a dot product of a photon wavefunction and an electron wavefunction treated as identical vectors. Contributions of indistinguishable pathways add as amplitudes within that matrix element. Taking its squared modulus then includes interference terms.

A reduced response model in a rotating frame is 

$$
\dot z_j=-(\gamma_j+\mathrm i\Delta_j)z_j+
 \eta_j\psi_j s(t)[1+\beta\xi_j(t)],\qquad z_j(0)=0.
\tag{2}
$$

 Here $s(t)$ is an imposed dimensionless temporal envelope, $z_j$ a dimensionless response amplitude, $\gamma_j$ an amplitude damping rate, $\Delta_j$ a detuning, and $\eta_j$ a coupling rate. The fluctuation strength $\beta$ and the complex process $\xi_j$ are dimensionless. Thus every term on the right has units of inverse time. The numerical units used in prior studies are arbitrary; no optical-cycle calibration is implied.

The squared modulus $|z_j|^2$ is called *excitation weight* or *squared response*. It is not automatically energy. A physical energy interpretation requires a normalization, transition spectrum, and material model. Similarly, a normalized single-photon amplitude is not a classical electric field carrying independently divisible photon energies at every site.

The fluctuation factor contains both constructive and destructive contributions: 

$$
|1+\beta\xi|^2=1+2\beta\operatorname{Re}\xi+\beta^2|\xi|^2.
\tag{3}
$$

 The cross term can have either sign; the total squared response cannot be negative. Describing that cross term as "negative power" would obscure the distinction.

For zero-mean proper complex Gaussian noise we specify 

$$
\mathbb{E}[\xi_j(t)\xi_k(t')]=0,\qquad
 \mathbb{E}[\xi_j(t)\xi_k(t')^*]=K_{jk}(t,t').
\tag{4}
$$

 $K$ must be positive semidefinite. Negative real cross-correlations are permissible within that constraint. A classical Gaussian process can serve as a phenomenological field model; it is not automatically equivalent to all ordered quantum-vacuum correlations. Establishing a QED origin requires the appropriate field operators, reservoir spectrum, ordering, and fluctuation--dissipation relation. Vacuum and material effects already included in a calibrated quantum baseline must not be counted again as additional noise.

## An explicit event law and the Born factorization

For the successful reference model, assume a conditional recording intensity 

$$
\lambda_j(t)=\alpha |z_j(t)|^2,
\tag{5}
$$

 where $\alpha$ has units of inverse time. Conditional on the field history, the channels have independent Poisson clocks. The first event is accepted and the trial ends. This gives one reported site by construction; a microscopic explanation of exhaustion of the shared excitation and inhibition of additional records is still required.

Let $\Lambda(t)=\sum_j\lambda_j(t)$. Conditional survival and the site-resolved first-event density are 

$$
S(t)=\exp\left[-\int_0^t\Lambda(u)\mathrm du\right],\qquad
 f_j(t)=\lambda_j(t)S(t).
\tag{6}
$$

 For random field histories these expressions must be averaged *after* constructing survival. Replacing each intensity by its mean generally gives a different model.

For identical deterministic responses, Eq. (2) factorizes as $z_j(t)=\psi_j F(t)$. Defining $h(t)=\alpha|F(t)|^2$ and $H(T)=\int_0^T h(t)\mathrm dt$ yields 

$$
P(J=j,t\leq T)=w_j[1-e^{-H(T)}],\qquad
 P(J=j\mid t\leq T)=w_j.
\tag{7}
$$

 This holds throughout response buildup, not only at steady state. A perfectly common multiplicative fluctuation also preserves this factorization on each history and after averaging. Unequal susceptibilities generally produce detector-response-weighted probabilities instead.

Equation (7) is an exact mathematical consequence of the model. It identifies both its success and its limit: complex linear response supplies amplitude-squared scaling, while Eq. (5) and the first-event construction supply the stochastic record. The calculation is not a derivation of the event law from vacuum dynamics. Encoding independent exponential thresholds as squared complex Gaussian variables would reproduce the same clocks but would not by itself provide that physical derivation.

## Results from the existing numerical studies

The following results are from saved synthetic calculations, not laboratory observations. Primary artifacts and their limitations are indexed in the companion review notes.

1.  **Transient complex response.** Fourteen cases of 6,000 incident histories each, totaling 84,000, tested spatial profiles, damping, phase factors, fluctuations, and detuning. With 256 channels, the recording coefficient was $\alpha=0.6$, the maximum gate was 6, and the principal integration interval was 0.05. Uniform deterministic cases obey the exact factorization above. The noisy cases without fixed response disorder were compatible with Born weights under the recorded multiple-comparison procedure; this finite-sample result does not establish exact equality. Two detuning cases departed from incident weights but agreed with independently integrated response-aware predictions. All 13 nonempty common-threshold comparisons failed the recorded distribution tests.

2.  **Phase-arrival competition.** A separate Adler-type phase race can approximate Born frequencies when its phase speed scales with incident intensity. A 200,000-trial follow-up nevertheless resolved a bias: the exact distribution differed from Born weights by total variation $0.01671$ for alignment tolerance $\pi/16$. The sample agreed with the exact race distribution. This failure is not a failure of the linear complex-response model, which is a different dynamical mechanism.

3.  **Reversible coupling.** An ideal single-excitation calculation with two coherent coupling pulses showed how a relative phase kick between them changes the final detector excitation. Reversing the entire closed evolution restored the input state. The 120,000 trial study used standard stochastic readout; it did not independently generate unique outcomes from unitary evolution.

4.  **Quantum absorption comparison.** A total of 288,000 sampled records compared coherent absorption with an undepleted response approximation. For amplitude damping $\gamma=1$ and zero detuning, the largest arrival-CDF difference was approximately $0.000915$ at coupling $g=0.05$, but $0.33124$ at $g=3$. The large discrepancy is explained by neglected depletion and coherent return. Site labels in this comparison were explicitly sampled from the common analytic Born law and are not additional evidence for its emergence.

The complex-response validation reconciled all 84,000 primary event histories. Identical-response capture probabilities agreed with analytic values to less than $4\times10^{-8}$. Paired step refinement checked selected histories, not every possible parameter. Spatial goodness-of-fit assessments used multinomial Monte Carlo reference distributions; multiple comparisons were adjusted in the saved analysis. A non-rejection is compatibility, not proof. The preliminary investigations were exploratory rather than preregistered, and should not be pooled into a single claimed significance level.

The documented Claude precursor replaces sine coupling by a linearly driven oscillator with steady response 

$$
P_i=\frac{|\kappa_i|^2}{(\omega_O-\omega_i)^2+\gamma^2}.
\tag{8}
$$

 Equal denominators and $\kappa_i$ proportional to the input amplitude give Born power ratios; separate Poisson rates select a record. The precursor explicitly leaves resonance tuning and the physical coupling open. Different-frequency cross terms vanish only under suitable averaging or channel separation. Linearity alone does not eliminate interference.

## A matched standard quantum comparator

In a coherent one-excitation benchmark, no-record amplitudes satisfy 

$$
\dot a=-\mathrm ig b,\qquad
 \dot b=-(\gamma+\mathrm i\Delta)b-\mathrm ig a.
\tag{9}
$$

 With monitored loss operator $L=\sqrt{2\gamma}|R\rangle\langle e|$, the unconditional recording density is $2\gamma|b|^2$ and the no-record norm is $|a|^2+|b|^2$. The conditional hazard divides the density by that norm. By contrast, the undepleted model $\dot z=-(\gamma+\mathrm i\Delta)z+g$ with rate $2\gamma|z|^2$ retains a prescribed drive until an event. It need not agree when depletion matters. Such disagreement is not evidence for new physics.

For the proposed two-channel comparison, a simpler exact baseline is useful. After controlled preparation of a single excitation, let $|1\rangle$ and $|2\rangle$ be two localized excitation states with initial weights $w$ and $1-w$. There is no coherent transfer between channels during the recording interval. Conditional on positive gains $G_j$, take 

$$
L_j(t)=\sqrt{h(t)G_j}|R_j\rangle\langle j|.
\tag{10}
$$

 Then standard no-record evolution gives 

$$
S_Q(t\mid G)=\sum_j w_j e^{-G_j H(t)},\qquad
 f_{Q,j}(t\mid G)=w_j h(t)G_j e^{-G_j H(t)}.
\tag{11}
$$

 If the two gain marginals are identical, ensemble site probabilities conditional on recording by any gate are $w_j$, regardless of their mutual correlation. At complete capture they are $w_j$ even for unequal strictly positive gains. Calibration losses and unequal finite-gate efficiencies must be included before making that comparison in an apparatus.

This baseline includes fluctuating detector rates within standard quantum theory. It avoids equating every observed fluctuation effect with a departure from quantum mechanics. The no-record update in Eq. (11) is crucial: unrecorded branch weights evolve, rather than remaining fixed drives.

## A conditional prediction from slow complex fluctuations

To turn the proposed mechanism into a falsifiable test, we now specify an additional, restricted variant. This is a *new analytic extension for this paper*, not a result claimed for all previous noise simulations. Suppose the field is nearly constant over one recording interval, renewed between trials, and produces 

$$
G_j=\frac{|1+\beta\xi_j|^2}{1+\beta^2},\qquad
 \mathbb E|\xi_j|^2=1,\qquad \mathbb EG_j=1.
\tag{12}
$$

 The normalization holds mean response fixed to separate mean sensitivity from fluctuation effects. It is a declared model choice, not a physical law.

The hypothesized first-event race retains rates $\lambda_j=h(t)w_jG_j$. Set $Q=\sum_jw_jG_j$. Its predictions are 

$$
\begin{aligned}
S_R(t\mid G)&=e^{-QH(t)}, \qquad (13) \\
P_R(J=j,t\leq T)&= \mathbb{E}\left[\frac{w_jG_j}{Q}(1-e^{-QH(T)})\right]. \qquad (14)
\end{aligned}
$$

 For $H(T)\rightarrow\infty$, $P_R(J=j)=\mathbb{E}[w_jG_j/Q]$. This differs structurally from Eq. (11). The distinction arises from the hypothesized competition law, not merely from assigning reality to complex waves.

For perfectly common fluctuations $G_1=G_2$, both models give Born site weights. For independent proper complex Gaussian $\xi_j$, write $q=\beta^2/(1+\beta^2)$ and $m=1-q$. The gain Laplace transform is 

$$
\mathcal L(s)=\mathbb Ee^{-sG}=\frac{\exp[-ms/(1+qs)]}{1+qs}.
\tag{15}
$$

 The race has 

$$
\begin{aligned}
S_R(t)&=\mathcal L(wH)\mathcal L((1-w)H), \qquad (16) \\
P_R(J=1)&=\int_0^\infty w[-\mathcal L'(ws)]\mathcal L((1-w)s)\mathrm ds, \qquad (17)
\end{aligned}
$$

 whereas the matched quantum model has $S_Q(t)=\mathcal L(H)$ and site probability $w$. These differences remain even when both models use the same gain distribution and mean response.

| Strength $\beta$ | Independent race $P_R(1)$ | Quantum or common-gain race |
|:----------------:|:-------------------------:|:---------------------------:|
|        0         |         0.800000          |          0.800000           |
|       0.10       |         0.798088          |          0.800000           |
|       0.25       |         0.788383          |          0.800000           |
|       0.50       |         0.760852          |          0.800000           |
|       1.00       |         0.727239          |          0.800000           |
|     $\infty$     |         0.717203          |          0.800000           |

: New conditional predictions for $w=0.8$ and complete capture. These are deterministic evaluations of Eq. (17), not observed counts or fitted experimental results. The infinite-strength row is a fluctuation-dominated mathematical limit.

In that limit the gains are unit-mean exponentials. For $a=w$, $b=1-w$, and $a\ne b$, 

$$
P_R(1)=\frac{a}{a-b}-\frac{ab}{(a-b)^2}\log\frac{a}{b}.
\tag{18}
$$

 At $w=1/2$ the continuous limit is $1/2$. Equal input weights thus conceal the spatial discriminator; unequal weights are essential. Independent 256- and 512-node quadrature evaluations agreed within $1.5\times10^{-15}$ on the tabulated parameter grid, and matched the zero-noise and exponential closed forms within $6\times10^{-16}$. These are numerical consistency checks, not uncertainty estimates for a physical detector.

Slow independent gains can therefore spoil Born frequencies even when their mean excitation weights are correct. This does not contradict finite-sample compatibility in the earlier time-dependent simulations: their correlation structure and temporal averaging differ. It also gives a useful constraint on the postulate. A viable mechanism must explain why ordinary detectors remain within observed Born statistics and specify the conditions under which a departure is expected.

The race is currently an imposed stochastic commitment law with a global first-record stop. It is not a derived many-body evolution. Its nonlinear dependence on input weights cannot simply be presumed to define a valid quantum measurement instrument. Extension to entangled inputs would require explicit tests of ensemble consistency and no-signalling. Failure of those tests could rule out this variant before an experiment.

## Proposed two-channel correlation experiment

The primary target is the joint distribution of channel label and recording time, including no-record trials. The intended platform is a controllable single-excitation system with two separately readable absorption or storage channels, such as a cavity or superconducting circuit implementation. A simple pair of ordinary photon counters would require additional modelling of absorption and amplification and may not expose the intermediate stage. The following is a proposed protocol, not a claim that a suitable apparatus has already been designed.

1.  **Prepare and verify the input.** Create a single excitation coherently distributed between the two channels. Use several weights, for example $w=0.2,0.35,0.5,0.65,0.8$, verify state preparation independently, and swap physical channel labels. Record source failures and multiple-excitation contamination.

2.  **Calibrate the quantum instrument.** Measure each channel's frequency, decay, transfer, efficiency, dark records, and timing response. Check coherent cross-coupling and common-reservoir collective decay. If present, these terms must enter the standard model; Eq. (10) alone would no longer be an adequate baseline.

3.  **Validate fluctuation control.** First use known engineered complex modulation to check the analysis and effective gains. This is an apparatus control, not evidence for vacuum selection. Measure gain marginals, spatial covariance, and memory. The slow-gain limit requires a correlation time longer than the capture window and a verified reset or adequate separation between trials.

4.  **Change the environment.** Vary a boundary condition or effective impedance to investigate a possible reservoir-linked contribution. Boundary modifications can change ordinary decay and level shifts, so repeat the relevant calibration. Demonstrate that the intended correlation change is achievable; do not infer it from mirror distance alone or assume that Casimir forces determine the noise spectrum.

5.  **Freeze and test.** Estimate the proposed mechanism's finite set of parameters on calibration settings. Before examining held-out data, freeze the model, gain normalization, temporal correlation law, gate choices, and statistical analysis. Predict both models at new weights, delays, and boundary settings, retaining all no-record and ambiguous trials.

The sharpest ideal comparison is common versus independent gains with matched marginals. For the slow race with $\beta=0.25$ and $w=0.8$, the complete-capture prediction changes from $0.800000$ to $0.788383$, while the matched quantum comparator stays at $0.8$. For finite gates, use Eq. (14); substituting the complete-capture value would bias the test. Timing survival provides a second, jointly predicted observable, rather than an independently fitted curve.

A rough planning estimate illustrates the scale. For a two-arm probability difference $\delta=0.011617$, independent trials, probability near $0.8$, two-sided 5% significance and 80% power, the normal approximation gives 

$$
n_{\rm arm}\simeq\frac{2(0.8)(0.2)(1.96+0.84)^2}{\delta^2}
 \simeq 1.9\times10^4
\tag{19}
$$

 usable records per arm. This is not a final sample-size calculation. Losses, parameter uncertainty, drift, correlations between trials, and multiple tests increase requirements. The final design should be powered by simulation of the complete calibrated likelihood. For small exploratory datasets, report raw counts and distributions rather than asserting statistical agreement.

## What the experiment could establish

Agreement with the standard baseline, with adequate precision, would constrain or reject the specified slow independent-gain race in the tested parameter range. It would not reject complex-wave realism or every possible fluctuation-based mechanism. Conversely, agreement with the race on held-out settings would support that quantitative model relative to the calibrated comparator. It would still require exclusion of ordinary response drift, technical gain noise, input miscalibration, collective decay, selection bias, and an incorrect quantum instrument.

In particular, the proposed race already predicts departures because its no-record update differs from standard quantum absorption. Observing that a simulation departs is not evidence that a real detector does. A successful experimental account must explain why the hypothesized update is physically appropriate. The experiment tests an explicit dynamical alternative, not "decoherence" in isolation.

Attributing any residual specifically to the QED vacuum would require a reservoir calculation connecting the observed boundary dependence and temporal spectrum to the fluctuation field in Eq. (2). A tunable extra noise function fitted separately at every setting would not be a discriminating prediction. If the completed mechanism predicts exactly the same observable distributions as the standard instrument, this experiment cannot distinguish their interpretations.

Objective-collapse proposals, including gravity-related ones, require their own evolution laws and parameter-dependent predictions before comparison [\[6\]](https://arxiv.org/abs/1204.4325). No collapse-rate or gravitational parameter is imported into the present model merely because it is a possible numerical timescale. Agreement between AI systems is similarly useful as an audit of calculations and assumptions, not as experimental evidence or proof of novelty.

## Remaining microscopic work

The central outstanding task is to obtain field fluctuations, local response, commitment, and amplification from a consistent joint dynamics. The present successful model retains separate stochastic recording variables. Removing those variables requires a definite field-driven event rule and a demonstration that it preserves the desired statistics across input profiles, channel numbers, observation gates, and detector responses. A single fitted histogram is insufficient.

The next theoretical checks are positivity and normalization of the proposed evolution, energy accounting, dark-event behaviour, finite propagation and detector back-action, and treatment of collective material response. The first-event stop must be replaced or justified by single-excitation dynamics. Amplification should use the detector's external energy supply, rather than assigning that energy to vacuum fluctuations. The model must also distinguish an overall electronic phase from physically relevant transition coherences.

The physical-wave postulate remains a motivation throughout these checks. Its value lies in guiding a concrete dynamics that survives them, not in changing the name of a random-number generator. The results so far identify a successful phenomenological reference model and a falsifiable extension, while retaining an explicit boundary between mathematical consistency, microscopic interpretation, and experimental support.

## Contributions and research transparency 

John Bramble proposed the physical-wave interpretation, the detector-local fluctuation mechanism, the focus on individual events, and successive questions directing the investigation. Claude contributed earlier equation development in saved conversations, including the linear-response and separate-race formulation. Codex implemented and analysed the numerical studies reported here, retrieved the precursor equations, and drafted this manuscript and its new conditional correlation calculation. The human investigator has requested the provisional repository contributor order Codex, Claude, John Bramble. Claude's review of this manuscript remains pending; no independent external AI or human peer review is claimed. Exact model identifiers should be added only where recoverable from session metadata, not reconstructed from memory.

This byline records collaboration in a working repository draft. It does not imply that AI systems accept legal or scientific accountability, or that a conventional journal will accept them as authors. The human contributor must verify and take responsibility for any submission or release. No publication, submission, or repository upload accompanies this draft.

## Data and code availability 

All reported event studies are synthetic. Saved programs, raw event logs, result summaries, and numerical validation records exist in the local research archive. The companion review notes identify the relevant directories. The new conditional predictions are supplied as a CSV with a standalone checker and validation record. These materials have not yet been assigned a public repository release or archival identifier. No raw experimental dataset was fitted in this paper. Original source recordings and transcripts were preserved without modification.

## References 

1.  M. Schlosshauer, "Decoherence, the measurement problem, and interpretations of quantum mechanics," *Rev. Mod. Phys.* **76**, 1267--1305 (2004). [arXiv:quant-ph/0312059](https://arxiv.org/abs/quant-ph/0312059).

2.  Y. Castin, J. Dalibard, and K. Mølmer, "A Wave Function approach to dissipative processes," 1992 manuscript posted in 2008. [arXiv:0805.4002](https://arxiv.org/abs/0805.4002).

3.  S. Scheel and S. Y. Buhmann, "Macroscopic QED---concepts and applications," *Acta Physica Slovaca* **58**, 675 (2008). [arXiv:0902.3586](https://arxiv.org/abs/0902.3586).

4.  I.-C. Hoi et al., "Probing the quantum vacuum with an artificial atom in front of a mirror," *Nature Physics* **11**, 1045 (2015). [arXiv:1410.8840](https://arxiv.org/abs/1410.8840); [doi:10.1038/nphys3484](https://doi.org/10.1038/nphys3484).

5.  N. Katz et al., "Reversal of the Weak Measurement of a Quantum State in a Superconducting Phase Qubit," *Phys. Rev. Lett.* **101**, 200401 (2008). Preprint titled "Uncollapsing of a quantum state in a superconducting phase qubit," [arXiv:0806.3547](https://arxiv.org/abs/0806.3547).

6.  A. Bassi, K. Lochan, S. Satin, T. P. Singh, and H. Ulbricht, "Models of wave-function collapse, underlying theories, and experimental tests," *Rev. Mod. Phys.* **85**, 471--527 (2013). [arXiv:1204.4325](https://arxiv.org/abs/1204.4325).

[^1]: Credited for earlier model development documented in the research notes. Claude has not reviewed or approved this manuscript. This is a provisional repository contributor byline, in the order requested by the human investigator, rather than a claim of compliance with journal authorship policies.

[^2]: Human investigator and proposed corresponding contributor. Identity and publication details are to be confirmed before release.
