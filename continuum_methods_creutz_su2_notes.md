# Continuum Methods, Creutz Ratios, and SU(2) Gauge-Breaking Context

Source paper: Andreas Risch, "Gauge field smearing and controlled continuum extrapolations", arXiv:2410.02794v1, DESY-24-074.

Local file: `2410.02794v1.pdf`

DESY PDF/A URL: https://bib-pubdb1.desy.de/record/640130/files/2410.02794v1.pdf?subformat=pdfa

arXiv URL: https://arxiv.org/abs/2410.02794

## Main Point

The paper studies how gauge-field smearing affects continuum extrapolations. Its central warning is that smearing reduces ultraviolet noise and lattice artifacts, but too much smearing can distort short-distance physics enough that continuum extrapolations from accessible lattice spacings become uncontrolled.

The numerical study is performed in pure SU(3) Yang-Mills theory, not SU(2). The methods are nevertheless directly relevant to SU(2) lattice-gauge analysis because they concern gauge-invariant Wilson-loop observables, Creutz ratios, and continuum scaling. For an SU(2) study with explicit gauge breaking, the main lesson is: Creutz ratios can diagnose whether gauge-breaking or smearing-like deformations modify the path to the continuum, especially at short distances.

## Gauge-Field Smearing in the Action

The paper starts from the common setup where the gauge action is left in terms of the original links, while the Dirac operator is evaluated on smeared links:

```tex
S[U] = S_g[U] + \bar{\Psi} D[\mathcal{S}[U]] \Psi .
```

Here `\mathcal{S}: U -> \mathcal{S}[U]` is a smearing transformation.

Examples mentioned:

- HYP smearing
- stout smearing
- HEX smearing
- gradient-flow smearing

Benefits of gauge-field smearing:

- Reduces ultraviolet fluctuations.
- Reduces the likelihood of small Dirac-operator eigenvalues.
- Helps avoid exceptional configurations.
- Can improve the spectral gap of Wilson-type Dirac operators.
- Can move improvement coefficients, such as `c_SW`, closer to their tree-level values.
- Can reduce the size of renormalization factors such as `Z_V`.

Risk:

- Too much smearing changes the ultraviolet structure of the lattice theory.
- Short-distance observables become especially vulnerable.
- A continuum extrapolation using not-small-enough lattice spacings can become misleading.

The paper first studies smeared observables in pure gauge theory:

```tex
<O_\mathcal{S}[U]> = <O[\mathcal{S}[U]]> .
```

This avoids the additional complication of dynamical fermions while still probing the physical gauge force through Wilson loops and Creutz ratios.

## Gradient Flow as a Continuum Method

The continuum Yang-Mills action is written as

```tex
S_YM = -1/(2 g_0^2) \int d^4x tr(F_{\mu\nu}(x) F_{\mu\nu}(x)).
```

The gradient-flow field is `B_\mu(x,t_fl)`, with initial condition

```tex
B_\mu(x,0) = A_\mu(x).
```

The flow equation is

```tex
\partial_{t_fl} B_\mu(x,t_fl)
  = - \delta S_YM[B] / \delta B_\mu(x,t_fl)
  = D_\nu G_{\nu\mu}(x,t_fl).
```

Definitions:

- `G_{\mu\nu}` is the field-strength tensor built from the flowed field `B`.
- `D_\mu = \partial_\mu + [B_\mu, .]` is the covariant derivative built from `B`.
- `t_fl` is the gradient-flow time.

At leading perturbative order, the gradient flow smooths the gauge field over the radius

```tex
r_sm = sqrt(8 t_fl).
```

Important continuum fact:

- Observables built from flowed fields at strictly positive flow time are finite once the underlying four-dimensional theory is renormalized.
- No additional renormalization is needed for such flowed observables.

The paper uses Wilson gradient flow as the lattice discretization and integrates it with a third-order Runge-Kutta scheme with step size

```tex
Delta t_fl / a^2 <= 0.01.
```

## Two Flow/Smearing Scenarios

The paper distinguishes two conceptually different uses of the gradient flow.

### 1. Gradient-Flow Smearing

The flow time is fixed in lattice units:

```tex
8 t_fl / a^2 = const.
```

Then, as `a -> 0`,

```tex
t_fl -> 0,
r_sm = sqrt(8 t_fl) -> 0.
```

So the continuum theory is unchanged. Smearing only changes the route toward the continuum limit.

This is the scenario closest to ordinary link smearing in a lattice action or measurement.

### 2. Physical Gradient Flow

The flow time is fixed in physical units:

```tex
t_fl / t_0 = const.
```

Then the smoothing radius stays finite in the continuum limit. This changes the continuum observable itself. The paper interprets this not as merely improving the same observable, but as defining a new flowed observable.

## Combined Continuum and Small-Flow-Time Expansion

The paper considers a finite dimensionless observable `\hat{O}`. It introduces

```tex
\hat{a} = a / sqrt(8 t_0),
\varepsilon = t_fl / t_0.
```

Since the observable is finite, the continuum and zero-flow-time limits can be interchanged:

```tex
lim_{\hat{a}->0} lim_{\varepsilon->0} \hat{O}
=
lim_{\varepsilon->0} lim_{\hat{a}->0} \hat{O}.
```

The common limit is `a = 0` and `t_fl = 0`. This allows a combined Symanzik and small-flow-time expansion:

```tex
\hat{O} = \sum_{i,j >= 0} c_{ij} \hat{a}^i \varepsilon^j.
```

The paper neglects logarithmic effects in both `a` and `t_fl` because the analysis has intermediate precision.

At fixed physical flow time, the continuum limit is

```tex
\hat{O}(\hat{a}=0,\varepsilon)
  = c_{00} + \sum_{j>0} c_{0j} \varepsilon^j.
```

So physical flow can change the continuum observable.

For fixed smearing strength,

```tex
\varepsilon / \hat{a}^2 = 8 t_fl / a^2.
```

The same expansion can be rewritten as

```tex
\hat{O}
  = \sum_{i,j >= 0} c_{ij}
      \hat{a}^{i+2j}
      (8 t_fl / a^2)^j.
```

At `\hat{a}=0`, this gives

```tex
\hat{O} = c_{00}.
```

Therefore, in the smearing scenario, the continuum limit is independent of the fixed smearing strength by construction, provided the expansion remains valid. Smearing changes lattice artifacts, not the target continuum value.

Practical meaning:

- Data at several lattice spacings and several small flow times can be globally fitted.
- The fitted coefficients `c_ij` reconstruct the lattice-spacing dependence at any sufficiently small smearing strength or flow time.
- If large smearing cannot be described by the truncated expansion, the continuum extrapolation is no longer controlled.

## Lattice Setup in the Paper

The study uses pure SU(3) Yang-Mills ensembles with Wilson plaquette action and temporal open boundary conditions to reduce topology freezing.

The ensemble parameters are:

| ensemble | beta | T/a | L/a | a [fm] | L [fm] | t0/a^2 |
|---|---:|---:|---:|---:|---:|---:|
| sft1 | 6.0662 | 80 | 24 | 0.0820(5) | 1.968(12) | 3.990(9) |
| sft2 | 6.2556 | 96 | 32 | 0.0616(4) | 1.971(12) | 7.070(17) |
| sft3 | 6.5619 | 96 | 48 | 0.04031(26) | 1.935(12) | 16.52(6) |
| sft4 | 6.7859 | 192 | 64 | 0.03010(19) | 1.927(12) | 29.60(10) |
| sft5 | 7.1146 | 320 | 96 | 0.01987(13) | 1.908(12) | 67.94(23) |

The reference scale `t_0` is defined from the flowed clover action density:

```tex
E(x,t_fl)
  = -1/2 \sum_{\mu,\nu}
      tr(G_{\mu\nu}^{clv}(x,t_fl) G_{\mu\nu}^{clv}(x,t_fl)).
```

The implicit definition is

```tex
t_0^2 <E(x,t_0)> = 0.3.
```

The paper uses

```tex
t_0 = 0.0268(3) fm^2
```

obtained from the Sommer scale `r_0`, using `r_0 = 0.5 fm` for illustration.

## Creutz Ratios

Creutz ratios are built from rectangular Wilson loops and are finite in the continuum limit, which makes them suitable for studying continuum extrapolations.

Continuum Wilson loop:

```tex
W(r,t) = < tr( P exp( \oint_{\gamma(r,t)} dx_\mu A_\mu(x) ) ) >.
```

Lattice Wilson loop:

```tex
W(r,t)
  = < tr( \prod_{(x,\mu) in \gamma(r,t)} U_\mu(x) ) >.
```

Continuum Creutz-ratio definition:

```tex
\chi(r,t)
  = - \partial_t \partial_r log W(r,t).
```

Central-difference lattice definition:

```tex
\chi(t+a/2, r+a/2)
  = 1/a^2
    log[
      W(t+a,r) W(t,r+a)
      /
      (W(t,r) W(t+a,r+a))
    ].
```

This central-difference form has `O(a^2)` lattice artifacts.

Relation to the static force:

```tex
\chi(r,t) -> F_{\bar{q}q}(r)    as t -> infinity.
```

The paper focuses on diagonal Creutz ratios:

```tex
\chi(r) = \chi(r,r).
```

Dimensionless form:

```tex
\hat{\chi} = 8 t_0 \chi,
\hat{r} = r / sqrt(8 t_0).
```

Distances are evaluated at half-integer lattice positions:

```tex
r/a = 1.5, 2.5, ...
```

Flow-time choices:

```tex
8 t_fl / a^2 =
  0, 0.25, 0.5, ..., 2, 2.5, ..., 3.5, 4, 5, 6, 7, 8
```

for smearing, and

```tex
8 t_fl / a^2 = (8 t_0 / a^2) * 0.0146 * j,
j = 0, 1, ..., 4
```

for physical flow.

## What Smearing Does to Creutz Ratios

For the ensemble sft4, the paper plots `\hat{\chi}` and its relative variance

```tex
sqrt(v_{\hat{\chi}}) / \hat{\chi}
```

as functions of `\hat{r}` and `8 t_fl/a^2`.

Observed behavior:

- The short-distance `~1/r^2` behavior is smoothed for distances

```tex
r \lessapprox sqrt(8 t_fl).
```

- This modifies the path to the continuum and changes lattice artifacts.
- The effect becomes smaller at larger `r`, where the smearing radius is less intrusive.
- The relative variance grows with increasing distance.
- Gradient-flow smearing reduces the relative variance at all distances.
- The variance reduction saturates; smearing does not give unlimited noise reduction.

This is the important tradeoff:

- More smearing: better noise and often smoother data.
- Too much smearing: distorted short-distance physics and non-controlled continuum extrapolation.

## Interpolation of Creutz Ratios

To extrapolate `\hat{\chi}(\hat{r})` at fixed physical distance, the value must be known at the same `\hat{r}` on all ensembles. Because measured values exist only at discrete `r/a`, the paper interpolates

```tex
\chi a^2
```

as a function of `r/a`.

The interpolation model is a source of systematic error:

- At short distances, interpolation-model differences dominate the uncertainty.
- At large distances, statistical noise dominates because Wilson-loop signals degrade.

The paper focuses on

```tex
0.3 <= \hat{r} <= 0.6,
0.14 fm <= r <= 0.28 fm.
```

This range avoids the worst short-distance lattice artifacts while retaining usable statistical and systematic errors for the ensembles considered.

## Continuum Extrapolation Fit Ansatz

At fixed `\hat{r}`, the paper performs a global continuum extrapolation that combines smearing and physical-flow data.

The truncated flow expansion is

```tex
\hat{\chi}_{tr}
  = c_{00}
  + c_{20} \hat{a}^2
  + c_{40} \hat{a}^4
  + c_{01} \varepsilon
  + c_{21} \hat{a}^2 \varepsilon
  + c_{02} \varepsilon^2.
```

Using

```tex
\varepsilon = (8 t_fl/a^2) \hat{a}^2,
```

the same ansatz can be read as a smearing expansion:

```tex
\hat{\chi}_{tr}
  = c_{00}
  + [c_{20} + c_{01} (8 t_fl/a^2)] \hat{a}^2
  + [c_{40}
     + c_{21} (8 t_fl/a^2)
     + c_{02} (8 t_fl/a^2)^2] \hat{a}^4.
```

Equivalently,

```tex
\hat{\chi}_{tr}
  = c_{00}
  + c_{20} [1 + (c_{01}/c_{20})(8 t_fl/a^2)] \hat{a}^2
  + c_{40} [1
      + (c_{21}/c_{40})(8 t_fl/a^2)
      + (c_{02}/c_{40})(8 t_fl/a^2)^2] \hat{a}^4.
```

Interpretation:

- `c_00` is the unflowed, unsmeared continuum limit.
- In the smearing scenario, all fixed-smearing extrapolations share `c_00`.
- In the physical-flow scenario, the continuum value depends on `\varepsilon`.
- The dependence on physical flow is stronger at shorter distances.

The paper reports that including too-large smearing strengths or too-large physical flow times causes poor fit `p`-values. This signals breakdown of the truncated Symanzik/small-flow-time expansion.

Distance dependence:

- At larger `\hat{r}`, larger `8 t_fl/a^2` and larger physical flow times can still be fitted.
- At shorter `\hat{r}`, less smearing is tolerable.

## Monotonicity as a Control Criterion

The paper uses monotonicity of the continuum extrapolation as a loose but practical criterion for control.

Reason:

- With limited knowledge of higher-order lattice artifacts and logarithms, a simple monotonic extrapolation is easier to trust.
- Non-monotonic behavior in `\hat{\chi}(\hat{a})` suggests that the chosen lattice spacings and smearing strengths are outside the reliable low-order expansion regime.

The paper tracks the peak location

```tex
\hat{a}_{peak}^2
```

of the fitted extrapolation as a function of smearing strength `8 t_fl/a^2` for several distances `\hat{r}`.

Rule:

```tex
\hat{a}^2 < \hat{a}_{peak}^2
```

means the considered lattice spacing lies in the monotonic part of the extrapolation.

Main numerical conclusion:

For ensembles with

```tex
a <= 0.06 fm
```

and for reliable monotonic extrapolations at distances

```tex
r >= 0.14 fm
```

or

```tex
\hat{r} >= 0.3,
```

one should choose approximately

```tex
8 t_fl / a^2 \lessapprox 1.
```

For even shorter distances:

- Use smaller lattice spacings.
- Or use less smearing.

## Stout Smearing Compared With Gradient Flow

For dynamical fermions, applying full gradient flow inside the Dirac operator is computationally expensive. Stout smearing with a small number of iterations is often used instead.

The relation used in the paper is

```tex
8 t_fl / a^2 = 8 n rho,
```

where:

- `n` is the number of stout-smearing steps.
- `rho` is the stout-smearing parameter.

In the limit

```tex
n -> infinity,
rho -> 0,
8 n rho fixed,
```

stout smearing converges to gradient-flow smearing.

The paper compares stout and gradient-flow smearing at

```tex
8 t_fl / a^2 = 8 n rho = 1.
```

Observed behavior:

- Replacing gradient flow by stout smearing increases the absolute size of lattice artifacts.
- With three stout iterations, the stout-smeared Creutz ratio almost reproduces the gradient-flow result within errors.
- The maximum of the stout-smeared extrapolation is shifted to somewhat larger `\hat{a}^2`.
- Therefore, the gradient-flow bound is conservative.
- Even one stout step with `rho = 1/8` reproduces the qualitative findings.

## Relevance for SU(2) Gauge-Breaking Studies

The paper itself is SU(3), but the continuum-method logic is useful for an SU(2) project with explicit gauge breaking.

Your thesis context appears to study whether small explicit gauge-invariance violations vanish, persist, or modify continuum extrapolations in SU(2). Creutz ratios are useful here because they are built from Wilson loops and directly probe the static force/string tension sector.

### What Transfers Directly to SU(2)

The following ingredients are group-general:

- Wilson loops as closed products of link variables.
- Creutz ratios as logarithms of neighboring Wilson-loop combinations.
- The relation of large-time Creutz ratios to the static force.
- Symanzik-style continuum extrapolation in powers of `a`.
- The idea that short-distance observables are most sensitive to ultraviolet deformations.
- The use of dimensionless variables such as `a/r_0`, `r/r_0`, `a sqrt(sigma)`, or a flow scale.
- The warning that deformation parameters can modify lattice artifacts before they vanish in the continuum.

For SU(2), Wilson loops usually use normalized trace

```tex
W(R,T) = (1/2) <Tr product_C U_\mu(x)>.
```

The normalization cancels in the standard Creutz ratio, so the core formula is unchanged.

### Gauge Breaking Versus Smearing

Smearing is gauge covariant when applied correctly; it does not intentionally break gauge symmetry. Your gauge-breaking setup is different if the action/update rule contains explicit non-gauge-invariant terms controlled by parameters such as `epsilon1` and `epsilon2`.

Still, both situations can be analyzed with a similar continuum question:

```tex
Does the deformation only change lattice artifacts,
or does it change the continuum limit?
```

For smearing at fixed `8 t_fl/a^2`, the paper expects the same continuum limit:

```tex
\hat{O}(a -> 0, fixed smearing) = c_00.
```

For explicit gauge breaking, this is not guaranteed. A useful generalized fit logic would be

```tex
\hat{O}(a,\epsilon)
  = O_0
  + A a^p
  + B \epsilon
  + C a^p \epsilon
  + D \epsilon^2
  + ...
```

or, if the breaking is expected to vanish with the cutoff,

```tex
\hat{O}(a,\epsilon)
  = O_0
  + A a^p
  + B \epsilon a^q
  + C \epsilon^2 a^{q'}
  + ...
```

The key distinction:

- If a nonzero `B epsilon` term survives as `a -> 0`, the gauge-breaking deformation changes the continuum limit.
- If breaking appears only in terms multiplied by positive powers of `a`, it behaves like a lattice artifact.
- If fits become non-monotonic or unstable at short distance, the extrapolation may be uncontrolled even if the formal continuum limit exists.

### Creutz Ratios as a Gauge-Breaking Diagnostic

Creutz ratios are especially useful because they are sensitive to the force and string tension:

```tex
\chi(R,T)
  = log[
      W(R,T) W(R+1,T+1)
      /
      (W(R+1,T) W(R,T+1))
    ]
```

depending on indexing convention. The central-difference version in the paper is equivalent up to coordinate placement.

If large Wilson loops follow an area law,

```tex
W(R,T) ~ exp[-\sigma a^2 R T + perimeter terms + ...],
```

then the Creutz ratio approaches

```tex
\chi(R,T) ~ \sigma a^2
```

at large `R,T`.

Gauge breaking can show up as:

- A shift in the inferred string tension.
- A different static force curve.
- A changed Sommer scale.
- Stronger short-distance artifacts.
- Non-monotonic continuum extrapolations.
- Different continuum intercepts for different breaking parameters.
- Failure of Wilson-loop area-law behavior or poor plateaus in effective potentials.

### SU(2) Analysis Strategy Inspired by the Paper

For each gauge-breaking parameter set, compute:

- Wilson loops `W(R,T)`.
- Creutz ratios `\chi(R,T)`.
- Diagonal Creutz ratios `\chi(R,R)`.
- Static force estimates from large `T`.
- String-tension estimates.
- Sommer scale estimates such as `r_0/a`.

Then perform continuum comparisons at fixed physical distance, not only fixed lattice `R`.

Recommended structure:

1. Choose a reference physical scale, for example `r_0`, `sqrt(sigma)`, or another stable scale from the SU(2) data.
2. Convert to dimensionless variables:

```tex
\hat{a} = a/r_0
```

or

```tex
\hat{a} = a sqrt(\sigma),
\hat{r} = r/r_0.
```

3. Interpolate `\chi a^2` or the force to common `\hat{r}` values across ensembles.
4. Fit unbroken and broken data together using an ansatz with both lattice-spacing and breaking-parameter dependence.
5. Check whether all breaking parameters extrapolate to the same continuum intercept.
6. Inspect monotonicity as a low-bar reliability check.
7. Repeat at several distances, because short distances should be most sensitive to gauge breaking.

### Practical Continuum Fit Templates for SU(2)

For standard Wilson plaquette SU(2), leading cutoff effects are often expected to be even powers of `a` for many gauge-invariant observables, but the explicit breaking term can change this expectation depending on its operator content. Start conservatively and compare alternatives.

Minimal shared-intercept fit:

```tex
\hat{\chi}(\hat{a},\epsilon)
  = c_0
  + c_a \hat{a}^2
  + c_\epsilon \epsilon \hat{a}^{p}
  + c_{a\epsilon} \epsilon \hat{a}^{p+2}.
```

This tests whether breaking vanishes as `a -> 0`.

Separate-intercept test:

```tex
\hat{\chi}(\hat{a},\epsilon)
  = c_0
  + c_\epsilon \epsilon
  + c_a \hat{a}^2
  + c_{a\epsilon} \epsilon \hat{a}^2.
```

If `c_\epsilon` is statistically nonzero, the breaking appears to survive in the continuum observable.

Quadratic breaking test:

```tex
\hat{\chi}(\hat{a},\epsilon)
  = c_0
  + c_a \hat{a}^2
  + c_\epsilon \epsilon
  + c_{\epsilon2} \epsilon^2
  + c_{a\epsilon} \hat{a}^2 \epsilon.
```

This is useful if the sign of the breaking matters or if the data are symmetric under `\epsilon -> -\epsilon`.

Monotonic-control check:

- For each fixed `epsilon`, plot `\hat{\chi}` versus `\hat{a}^2`.
- Mark the coarsest lattice spacing included.
- Repeat after removing the coarsest ensemble.
- If the extrapolated intercept changes strongly or the curve turns over, treat the fit as uncontrolled.

### How This Relates to Your Codebase

The workspace already contains Creutz-ratio analysis in `calculator.py`, `run_evaluation.py`, and `finalize_analysis.py`.

Relevant implementation ideas:

- `calculator.py` defines the adjacent-loop Creutz ratio and derived quantities.
- `finalize_analysis.py` creates Creutz summaries and diagonal `R = T` plots.
- `search_data.py` can aggregate dynamic fields such as `creutz_P`, `a_creutz`, and their errors.
- The thesis draft already frames `epsilon` as the gauge-breaking parameter.

The paper suggests adding or emphasizing:

- Continuum plots of Creutz observables against a dimensionless lattice spacing such as `(a/r_0)^2`.
- Fits at fixed physical distance, requiring interpolation across ensembles.
- Explicit fit terms mixing lattice spacing with gauge-breaking parameters.
- A monotonicity/stability criterion for deciding whether a continuum extrapolation is controlled.
- Special caution for small `R` or small physical `r`, where gauge breaking and smearing artifacts are most visible.

## Key Takeaways for the Thesis

- Creutz ratios are a good bridge between Wilson-loop data and continuum force/string-tension physics.
- The paper's main continuum-method contribution is the combined expansion in lattice spacing and flow time.
- Fixed lattice-unit smearing should not alter the continuum limit if it remains in the controlled expansion regime.
- Physical flow changes the continuum observable by construction.
- Too much smoothing gives cleaner data but can produce non-monotonic, unreliable continuum extrapolations.
- Short-distance observables tolerate less deformation.
- For the paper's SU(3) setup, `8 t_fl/a^2 <= about 1` is recommended when using `a <= 0.06 fm` and distances `r >= 0.14 fm`.
- For SU(2) gauge breaking, use the same philosophy: test whether the deformation changes only cutoff effects or also the continuum intercept.
- A shared continuum limit across breaking parameters supports restoration/irrelevance of the breaking.
- Different continuum intercepts support a surviving gauge-breaking deformation.
- Non-monotonic or unstable continuum fits indicate that the accessible lattice spacings are too coarse, the breaking/smearing is too large, or the fit ansatz is missing important terms.

## Suggested Citation Entry

```bibtex
@article{Risch2024SmearingContinuum,
  author        = {Risch, Andreas},
  title         = {Gauge field smearing and controlled continuum extrapolations},
  eprint        = {2410.02794},
  archivePrefix = {arXiv},
  primaryClass  = {hep-lat},
  reportNumber  = {DESY-24-074},
  year          = {2024}
}
```
