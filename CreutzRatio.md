Gauge field smearing and controlled continuum
extrapolations
Andreas Risch𝑎,𝑏,∗
𝑎Department of Physics, University of Wuppertal, Gaussstr. 20, 42119 Wuppertal, Germany
𝑏John von Neumann-Institut für Computing NIC, Deutsches Elektronen-Synchrotron DESY,
Platanenallee 6, 15738 Zeuthen, Germany
E-mail: andreas.risch@uni-wuppertal.de
When designing lattice actions, gauge field smearing is often used in the definition of the lattice
Dirac operator. Too much smearing can result in uncontrolled continuum extrapolations as the
short distance behaviour of the theory is mutilated, which is a situation to be avoided. As a
smearing prescription we focus on the gradient flow formalism as it allows to study both smearing
and physical flow simultaneously. We investigate the effect of smearing and physical flow on the
scaling towards the continuum limit in pure gauge theory. We focus on the example of Creutz
ratios, which provide a measure of the physical forces felt by the fermions. For suitable smearing
strengths we further investigate the impact of replacing the Wilson gradient flow by stout smearing.
European network for Particle physics, Lattice field theory and Extreme computing (EuroPLEx2023)
11-15 September 2023
Berlin, Germany
∗Speaker
© Copyright owned by the author(s) under the terms of the Creative Commons
Attribution-NonCommercial-NoDerivatives 4.0 International License (CC BY-NC-ND 4.0). https://pos.sissa.it/
arXiv:2410.02794v1 [hep-lat] 18 Sep 2024
Gauge field smearing and controlled continuum extrapolations Andreas Risch
1. Introduction
A reduction of lattice artefacts is beneficial for more reliable continuum extrapolations, in
particular of short distance observables. A popular methods to alter discretisation effects is UV
filtering, which is based on the application of four-dimensional gauge field smearing. The Dirac
operator is evaluated on smeared gauge fields such that the action is altered into
𝑆[𝑈] = 𝑆g [𝑈] + Ψ 𝐷 [S [𝑈]] Ψ, (1)
where S : 𝑈 7 → S [𝑈] is a smearing transformation. Several smearing algorithms have been
developed, e.g. HYP [1], Stout [2], HEX [3] and gradient flow [4, 5] smearing. Evaluating the
Dirac operator on smeared gauge fields yields several advantages: The likelihood of finding small
eigenvalues of 𝐷 is reduced, i.e. exceptional configurations can be avoided. In [6] even at very
coarse lattice spacings the Wilson Dirac operator defined with nHYP gauge links could be shown
to exhibit a spectrum with a well-defined spectral gap. The same was shown for stout smearing
in [7]. This is particularly helpful for the simulation of mass non-degenerate quarks as the fermion
determinant is not necessarily positive in such a scenario [8]. Gauge field smearing also has an
impact on improvement coefficients and renormalisation constants. In [6] it was observed that the
improvement coefficient 𝑐SW approaches its tree-level value when gauge field smearing is applied.
The amount of renormalisation in 𝑍V is also reduced. However, the application of too much
smearing may significantly alter the UV structure of the lattice theory and therefore continuum
extrapolations based on data from insufficiently small lattice spacings may become unreliable. It is
therefore relevant to study the range of smearing strengths that still allow for controlled continuum
extrapolations. As a first step towards a smeared action setup with fermions we study smeared
observables
〈𝑂S [𝑈]〉 = 〈𝑂 [S [𝑈]]〉 (2)
in pure gauge theory. We investigate the influence of smearing on continuum extrapolations of
Creutz ratios [9], which provide a measure of the physical forces felt by the fermions caused by the
gauge field. For a previous account of this effort we refer the reader to [10, 11].
2. The gradient flow formalism, gradient flow smearing and physical gradient flow
In this work we focus on the gradient flow formalism [5] as a smearing procedure. We start
from the continuum four-dimensional Yang-Mills action 𝑆YM = − 1
2𝑔2
0
∫ d4𝑥 tr(𝐹𝜇𝜈 (𝑥)𝐹𝜇𝜈 (𝑥)).
𝐹𝜇𝜈 = 𝜕𝜇 𝐴𝜈 − 𝜕𝜈 𝐴𝜇 + [ 𝐴𝜇, 𝐴𝜈 ] denotes the field strength tensor and 𝐴𝜇 (𝑥) the corresponding
gauge field. In the gradient flow formalism a gauge field 𝐵𝜇 (𝑥, 𝑡fl) is introduced, where 𝑡fl ≥ 0
is the so called gradient flow time. At 𝑡fl = 0 the standard gauge field 𝐴𝜇 (𝑥) is used as an initial
condition for the flow time evolution, i.e. 𝐵𝜇 (𝑥, 0) = 𝐴𝜇 (𝑥). The evolution is then governed by the
gauge-covariant flow equation
𝜕
𝜕𝑡fl
𝐵𝜇 (𝑥, 𝑡fl) = − 𝛿𝑆YM [𝐵]
𝛿𝐵𝜇 (𝑥, 𝑡fl) = 𝐷𝜈 𝐺𝜈 𝜇 (𝑥, 𝑡fl), (3)
2
Gauge field smearing and controlled continuum extrapolations Andreas Risch
where 𝐺 𝜇𝜈 = 𝜕𝜇 𝐵𝜈 − 𝜕𝜈 𝐵𝜇 + [𝐵𝜇, 𝐵𝜈 ] denotes the generalised field strength tensor and 𝐷 𝜇 =
𝜕𝜇 + [𝐵𝜇, ·] the generalised covariant derivative. Performing a leading-order perturbative expansion
it was shown that the gauge field 𝐵𝜇 (𝑥, 𝑡fl) is a spherically smoothed version of 𝐴𝜇 (𝑥) with mean-
square radius 𝑟sm = √8𝑡fl [5], i.e. in the direction of positive flow time the gradient flow possesses
a smoothing property. In [12] it was shown perturbatively to all loop orders that any functional
of the flowed fields 𝐵𝜇 (𝑥, 𝑡fl) at strictly positive 𝑡fl is finite, assuming that the four-dimensional
theory has been renormalised. Consequently, no additional renormalisation has to be applied. The
Wilson gradient flow [5] is used as a lattice discretisation of the Yang-Mills gradient flow. The
flow equation is then integrated numerically using an explicit 3rd-order Runge-Kutta integration
scheme [5] with a step size Δ𝑡fl
𝑎2 never exceeding 0.01.
The gradient flow will be applied to the gauge field in two scenarios: In the first scenario,
which we refer to as gradient flow smearing, the gradient flow time and consequently the smearing
radius vanishes in the continuum limit. Hence the continuum theory is unaltered. This can be
achieved by fixing the gradient flow time in lattice units, i.e. 8𝑡fl
𝑎2 = const. The second scenario, in
which the flow time is fixed in physical units, i.e. 𝑡fl/𝑡0 = const, we refer to as a physical gradient
flow. In principle, 𝑡0 may be any physical scale of the theory. We make use of the reference flow
time introduced in [5], which we will define in section 4. In this scenario the continuum theory
is altered. This type of alteration of an observable’s continuum limit can also be understood as
a modification of the definition of the observable itself, i.e. the physical gradient flow allows to
construct new observables.
3. Combined continuum extrapolation and small flow time expansion
In the following we consider a dimensionless observable ˆ𝑂, which does not require a renormal-
isation and hence is finite in the continuum limit. We will understand this observable as a function
of the dimensionless lattice spacing parameter ˆ𝑎 ≡ 𝑎√8𝑡0
and the flow time parameter 𝜀 = 𝑡fl
𝑡0 . Due to
the finiteness of the observable the continuum limit and the zero flow time limit can be interchanged,
i.e. lim ˆ𝑎→0 lim𝜀→0 ˆ𝑂 = lim𝜀→0 lim ˆ𝑎→0 ˆ𝑂. In this case the two scenarios discussed in section 2
have a common limit where both 𝑎 = 0 and 𝑡fl = 0. Therefore, a combined Symanzik and small
flow time expansion is possible and well-defined. The double expansion of the observable reads
ˆ𝑂 =
Õ
𝑖, 𝑗 ≥0
𝑐𝑖 𝑗 ˆ𝑎𝑖 𝜀 𝑗 . (4)
We neglect logarithmic effects both in the lattice spacing [13] and in the flow time [5] as this
investigation has only intermediate precision. Evaluating this expression in the continuum ˆ𝑎 = 0, it
becomes obvious that the observable’s continuum limit ˆ𝑂 = 𝑐00 + Í𝑛
𝑗>0 𝑐0 𝑗 𝜀 𝑗 can be altered by a
physical gradient flow. 𝑐00 denotes the continuum limit at vanishing flow time. In this work, we are
primarily interested in the effect of smearing on the continuum extrapolation. To demonstrate that
eq. (4) also describes the observable’s lattice spacing dependence at fixed smearing strength 8𝑡fl
𝑎2 , we
observe that the latter is parametrised by 𝜀
ˆ𝑎2 = 8𝑡fl
𝑎2 . The expansion can therefore we rewritten as a
function of the lattice spacing and the smearing strength:
ˆ𝑂 =
Õ
𝑖, 𝑗 ≥0
𝑐𝑖 𝑗 ˆ𝑎𝑖+2 𝑗 ( 𝜀
ˆ𝑎2
) 𝑗
=
Õ
𝑖, 𝑗 ≥0
𝑐𝑖 𝑗 ˆ𝑎𝑖+2 𝑗 ( 8𝑡fl
𝑎2
) 𝑗
. (5)
3
Gauge field smearing and controlled continuum extrapolations Andreas Risch
Evaluating the smearing expansion in the continuum limit ˆ𝑎 = 0 yields ˆ𝑂 = 𝑐00, i.e. the continuum
limit is independent of the smearing strength by construction. The main advantage of this combined
Symanzik and small flow time expansion is that data measured at various small ˆ𝑎 ≡ 𝑎√8𝑡0
and 𝜀 = 𝑡fl
𝑡0
can be combined to determine the coefficients 𝑐𝑖 𝑗 , from which the lattice spacing dependence can
be reconstructed for any sufficiently small smearing strength or flow time parameter.
4. Lattice setup
ensemble 𝛽 𝑇/𝑎 𝐿/𝑎 𝑎 [fm] 𝐿 [fm] 𝑡0/𝑎2
sft1 6.0662 80 24 0.0820(5) 1.968(12) 3.990(9)
sft2 6.2556 96 32 0.0616(4) 1.971(12) 7.070(17)
sft3 6.5619 96 48 0.04031(26) 1.935(12) 16.52(6)
sft4 6.7859 192 64 0.03010(19) 1.927(12) 29.60(10)
sft5 7.1146 320 96 0.01987(13) 1.908(12) 67.94(23)
Table 1: Parameters of the SU(3) gauge ensembles [14] and computed reference flow time 𝑡0/𝑎2 in lattice
units.
This study is based on SU(3) Yang Mills theory gauge ensembles [14] using the Wilson
plaquette action, where temporal open boundary conditions [15] are imposed to alleviate topology
freezing. An overview of the gauge ensembles is given in table 1. The reference flow time 𝑡0 [5] is
used as a scale to construct dimensionless quantities. To define 𝑡0 we make use of the action density
𝐸 (𝑥, 𝑡fl) = − 1
2
Õ
𝜇,𝜈
tr (𝐺clv
𝜇𝜈 (𝑥, 𝑡fl) 𝐺clv
𝜇𝜈 (𝑥, 𝑡fl)), (6)
where 𝐺clv denotes the field strength tensor in the clover discretisation [16]. The reference flow
time 𝑡0 is then implicitly defined by [5]
𝑡2
0 〈𝐸 (𝑥, 𝑡0)〉 = 0.3. (7)
Numerical values are listed in table 1. The physical value of 𝑡0 = 0.0268(3) fm2 is obtained from
the force parameter 𝑟0 [17], where for illustration a value of 𝑟0 = 0.5 fm is used. The lattice spacing
varies between 0.08 fm and 0.02 fm and the spatial extent between 1.9 fm and 2 fm.
5. Creutz ratios and gradient flow
Creutz ratios [9] are suitable observables for a study in pure gauge theory as they possess a
finite continuum limit. The latter are constructed from planar rectangular Wilson loops 𝑊 (𝑟, 𝑡) ≡
〈tr(𝑃 exp(∮
𝛾 (𝑟 ,𝑡 ) d𝑥𝜇 𝐴𝜇 (𝑥)))〉, which are obtained from the gauge field by a path-ordered integral
along a rectangular closed path 𝛾(𝑟, 𝑡). In lattice gauge theory these objects are discretised as
𝑊 (𝑟, 𝑡) =
〈
tr
( Ö
( 𝑥,𝜇) ∈𝛾 (𝑟 ,𝑡 )
𝑈𝜇 (𝑥)
)〉
. (8)
4
Gauge field smearing and controlled continuum extrapolations Andreas Risch
0.0 0.2 0.4 0.6 0.8 1.0
ˆr
0
10
20
30
40
50
ˆχ
sft4
8tfl
a2
0
1
2
4
8
0 2 4 6 8 10
8tfl
a2
10−2
10−1
100
101
102
103
√vˆχ
ˆχ
sft4
ˆr
0.227
0.487
0.747 1.01
Figure 1: Dimensionless Creutz ratio ˆ𝜒 and relative variance
√𝑣 ˆ𝜒
ˆ𝜒 as functions of the flow time 8𝑡fl
𝑎2 and the
distance ˆ𝑟 on the ensemble sft4.
Creutz ratios are obtained from Wilson loops by 𝜒(𝑟, 𝑡) ≡ − 𝜕
𝜕𝑡
𝜕
𝜕𝑟 ln(𝑊 (𝑟, 𝑡)). To obtain 𝑂 (𝑎2)
lattice artefacts the latter definition is discretised making use of central differences [18]:
𝜒
(
𝑡 + 𝑎
2 , 𝑟 + 𝑎
2
)
≡ 1
𝑎2 ln
( 𝑊 (𝑡 + 𝑎, 𝑟) · 𝑊 (𝑡, 𝑟 + 𝑎)
𝑊 (𝑡, 𝑟) · 𝑊 (𝑡 + 𝑎, 𝑟 + 𝑎)
)
. (9)
The static quark anti-quark force can be extracted in the limit of an infinite time extent, 𝜒(𝑟, 𝑡) →
𝐹qq (𝑟) for 𝑡 → ∞ [18].
In the following discussion we will only focus on diagonal Creutz ratios 𝜒(𝑟, 𝑡) with 𝑟 = 𝑡,
which we abbreviate as 𝜒(𝑟) ≡ 𝜒(𝑟, 𝑟). We compute the latter in lattice units ( 𝜒 · 𝑎2)( 𝑟
𝑎 ) for various
half integer distances 𝑟
𝑎 = 1.5, 2.5, . . . based on gauge configurations which gradient flow smearing
was applied to. We use 𝑡0 to define dimensionless Creutz ratios, i.e. we analyse ˆ𝜒 ≡ 𝜒 · 8𝑡0 as a
function of ˆ𝑟 ≡ 𝑟√8𝑡0
. In our measurements we implement the two scenarios for scaling the flow
time via
8𝑡fl
𝑎2 =
{0, 0.25, 0.5, . . . , 2, 2.5, . . . , 3.5, 4, 5, 6, 7, 8 smearing
8𝑡0
𝑎2 × 0.0146 × 𝑗 , 𝑗 ∈ {0, 1, . . . , 4} physical flow. (10)
The computation is based on the openQCD [19] package and utilises B. Leder’s program for
measuring Wilson loops [20, 21]. For the data analysis the python3 package pyobs [22] is used,
which implements the Γ-method [23] for Monte Carlo error estimation.
As discussed in the introduction smearing is commonly used to reduce UV fluctuations in
gauge fields, which also has an impact on the variance of observables. In fig. 1 the dimensionless
diagonal Creutz ratio ˆ𝜒 and its relative variance
√𝑣 ˆ𝜒
ˆ𝜒 are displayed as functions of the distance ˆ𝑟 and
the smearing strengths 8𝑡fl
𝑎2 for the ensemble sft4. We observe that the ∼ 1
𝑟2 short distance behaviour
is smoothed by the gradient flow at distances 𝑟 / √8𝑡fl. Consequently, the path to the continuum
and hence lattice artefacts are altered in the smearing scenario. This effect becomes smaller at larger
distances where the smearing has less impact. We observe that the relative variance of the Creutz
ratio
√𝑣 ˆ𝜒
ˆ𝜒 grows with growing distances. Applying gradient flow smearing the relative variance
shrinks with growing flow time at all distances [18]. However, smearing the gauge fields does not
5




Gauge field smearing and controlled continuum extrapolations Andreas Risch
References
[1] A. Hasenfratz and F. Knechtli, Phys. Rev. D 64 (2001) 034504 [hep-lat/0103029].
[2] C. Morningstar and M. J. Peardon, Phys. Rev. D 69 (2004) 054501 [hep-lat/0311018].
[3] S. Capitani, S. Durr and C. Hoelbling, JHEP 11 (2006) 028 [hep-lat/0607006].
[4] R. Narayanan and H. Neuberger, JHEP 03 (2006) 064 [hep-th/0601210].
[5] M. Lüscher, JHEP 08 (2010) 071 [1006.4518].
[6] A. Hasenfratz, R. Hoffmann and S. Schaefer, JHEP 05 (2007) 029 [hep-lat/0702028].
[7] S. Dürr, Z. Fodor, C. Hoelbling et al., Phys. Rev. D 79 (2009) 014501 [0802.2706].
[8] D. Mohler and S. Schaefer, Phys. Rev. D 102 (2020) 074506 [2003.13359].
[9] M. Creutz, Phys. Rev. Lett. 45 (1980) 313.
[10] A. Risch, S. Schaefer and R. Sommer, PoS LATTICE2022 (2023) 384 [2212.04000].
[11] A. Risch, PoS LATTICE2023 (2024) 342 [2310.06587].
[12] M. Lüscher and P. Weisz, JHEP 02 (2011) 051 [1101.0963].
[13] N. Husung, P. Marquard and R. Sommer, Phys. Lett. B 829 (2022) 137069 [2111.02347].
[14] N. Husung, M. Koren, P. Krah et al., EPJ Web Conf. 175 (2018) 14024 [1711.01860].
[15] M. Lüscher and S. Schaefer, JHEP 07 (2011) 036 [1105.4749].
[16] B. Sheikholeslami and R. Wohlert, Nucl. Phys. B 259 (1985) 572.
[17] R. Sommer, Nucl. Phys. B 411 (1994) 839 [hep-lat/9310022].
[18] M. Okawa and A. Gonzalez-Arroyo, PoS LATTICE2014 (2014) 327 [1410.7862].
[19] http://cern.ch/luscher/openQCD/ .
[20] https://github.com/bjoern-leder/wloop .
[21] M. Donnellan, F. Knechtli, B. Leder et al., Nucl. Phys. B 849 (2011) 45 [1012.3037].
[22] https://mbruno46.github.io/pyobs/ .
[23] U. Wolff, Comput. Phys. Commun. 156 (2004) 143 [hep-lat/0306017].
[24] M. Nagatsuka, K. Sakai and S. Sasaki, Phys. Rev. D 108 (2023) 094506 [2303.09938].
10