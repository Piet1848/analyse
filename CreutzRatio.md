PHYSICAL REVIEW 0 VOLUME 23, NUMBER 8 15 APRIL 1981
Monte Carlo stufiy of renormalixation in lattice gauge theory
Michael Creutz
Department ofPhysics, Brookhaven National Laboratory, Upton, New Fork 11973
(Received 11 December 1980)
Using Monte Carlo methods with SU(2) and U(1) lattice gauge theories, we compare physical ratios of %'ilson loops
on different length scales. An interesting renormalization-group structure is suggested for the U(1) theory. For the
SU(2) case we only see a fixed point at vanishing bare coupling. At this point asymptotic freedom is verified
numerically.
I. INTRODUCTION
Monte Carlo methods have become powerful tools
in the study of nonperturbative phenomena in gauge
theories. ' ' A discrete space-time lattice, used as
an ultraviolet regulator in a Euclidean path inte-
gral, converts a quantized gauge theory into an
equivalent classical statistical system4 which is
well suited to numerical simulation using esta-
blished techniques. Possible phase transitions
play a key role in all discussions of lattice gauge
theory. "Indeed, the appearance of a transition
in a U(l) gauge theory is essential to the existence
of a massless photon in a lattice formulation of
electrodynamics. Conversely, the absence of
transitions in four-dimensional non-Abelian gauge
models is central to our understanding of quark
confinement in the continuum limit of these theor-
ies.
In earlier papers' we studied numerically the
long-range forces between external sources with
quark quantum numbers in lattice gauge theories
based on the groups SU(2) and SU(3). Using the
long-range force to define a renormalization
scheme, we found numerical evidence that a linear
potential at long distances survives in the weak-
coupling continuum limit of the theory. One short-
coming of that work arose from the slow logarith-
mic decrease of the bare coupling with decreasing
lattice size. This only allowed modestly weak
coupling while still seeing the asymptotic linear
potential on our lattices of only 6'-10' sites.
Furthermore, there are probably "roughening"
transitions giving a nonanalytic behavior to this
force at values of coupling where short-range cor-
relation functions show no singularities. This in-
dicates an awkwardness with the renormalization
scheme of holding the asymptotic force law fixed.
In this paper we present renormalization-group
arguments based on physical quantities defined on
finite length scales. This permits study of the
renormalization of the bare coupling over a wider
range of the cutoff parameter. Monte Carlo re-
suits for SU(2) gauge theory at weak coupling veri-
fy logarithmic dependence of the bare charge on
distance scale, as predicted by asymptotic free-
dom. Passing from weak to strong coupling, we
find no evidence for a renormalization-group fixed
point, further strengthening evidence that the con-
finement inherent in Wilson's strong-coupling ex-
pansion survives in the continuum limit.
To test these techniques on a system known to
have a phase transition, "we also study the U(l)
lattice gauge theory. Here asymptotic freedom is
lost and a nontrivial renormalization-group
structure appears. We conjecture that a continu-
um limit at the critical point gives an interacting
field theory of photons and magnetic monopoles.
Away from this critical point on the weak-coupling
side, a continuum limit should yield a field theory
of free photons with the monopole mass diverging
when the cutoff is removed.
'The essence of the procedure is to define a phys-
ical correlation function at some scale and de-
mand that this quantity not change upon reducing
the lattice spacing by a factor of 2. 'This relates
the bare coupling constants for two values of the
cutoff. Iteration then drives the bare coupling to
a fixed point relevant for the continuum limit. For
a non-Abelian gauge theory this yields the asymp-
totically free behavior of a logarithmic decrease
in the bare coupling.
We only consider a single bare-coupling para-
meter; in this respect our analysis is less sophis-
ticated than that proposed by Wilson' or, as used
for spin systems, that of Swendsen. ' Additional
coupling terms could in principle reduce the fi-
nite-cutoff ambiguities at the price of substantial-
ly increasing the computing time. Our analysis
will be entirely in the context of pure gauge theory
without quarks. This is because we do not know
how to treat fermion fields with Monte Carlo
techniques.
In Sec. II we describe the general renormaliza-
tion-group procedure without specifying the pre-
cise physical quantities used. Section III is where
1815
1816 MICHAEL CREUTZ
we introduce ratios of Wilson loops with different
shapes but the same perimeters. We argue that
ultraviolet divergences should cancel from such
ratios and therefore they can serve to define a
renormalized coupling constant. In Sec. IV we
present our Monte Carlo measurements of these
ratios. We observe asymptotic freedom for the
SU(2) theory and a nontrivial structure for a U(1)
gauge group. Section V contains some concluding
remarks.
II. LENGTH SCALES AND RENORMALIZATION
In particle physics, the goal of renormalization
is to remove ultraviolet divergences from a field
theory. I'he bare-coupling parameters become
functions of an ultraviolet cutoff in such a manner
that physical quantities have a finite limit as the
cutoff is removed. A renormalization scheme be-
gins with the selection of an arbitrary set of phys-
ical measurables which is sufficiently complete to
determine the bare parameters when the cutoff is
in place. Then, as the cutoff is removed, the
bare couplings are continuously adjusted in such
a manner that these given measurables remain
fixed. For a mell-defined renormalizable theory,
this procedure should yield unique finite limits for
all physical quantities.
'The observables used in a renormalization
scheme are quite unconstrained. In quantum elec-
trodynamics one usually fixes the physical elec-
tron mass and the coefficient of the long-range
Coulomb force. In a confining theory, such as
believed for the strong interaction, the choice is
less obvious. One popular selection for nonper-
turbative studies of pure gauge theory without
fermions is, the coefficient K of the hypothesized
long-distance linear potential between external
sources with quark quantum numbers. Another
possible choice is the mass of some physical
bound state, such as the lightest glueball.
All of the quantities mentioned in the previous
paragraph are defined in terms of long-range ef-
fects. 'This is clear for the long-distance poten-
tials, but it also applies to a particle mass as
this determines how the particle propagates over
an extended range. It is, however, also conven-
ient to consider physical observables involving on-
ly finite scales. For example, in traditiona1. per-
turbative renormaliz ation-group discussions one
studies vertex functions in momentum space with
all legs off shell at some given momentum scale
Alternatively, one might be interested in some
interparticle force at a finite range r. By varying
these parameters p, or r, one studies the inter-
relationships of physics on different length scales.
From now on we will restrict our discussion to
a theory, such as quarkless gauge theory, which
has only one bare dimensionless coupling para-
meter, g, . A general physical observable P is a
function of g„as well as the cutoff scale of length
a, and the scale r on which P is to be measured,
P = P(r, a,g,(a)) . (2.1)
Here we have explicitly shown the cutoff depen-
dence of the bare coupling g,(a). The precise
form of this dependence will depend on the re-
normalization scheme used. For simplicity, as-
sume that P is dimensionless; if it were not, just
multiply by enough powers of r to make it so. For
example, from an interparticle force E(r) con-
struct P =r'E.
As a becomes small and we approach the con-
tinuum limit, P should lose cutoff dependence.
Thus we expect
P~ r, —,g, —~ ~
=P(r, c,g,(a))+O(a'),
a al) (2.2)
where we have arbitrarily compared cutoffs dif-
fering by a factor of 2. In general there are two
classes of dimensional parameters which set the
scale for the order-a' corrections in Eg. (2.2).
First, of course, is the scale r used to define P.
In addition we must consider the long-range phys-
ical parameters characterizing the continuum
theory. In particular, regardless of how large r
is, we expect corrections to Eq. (2.2) or order
a'-m' where m is a typical mass in the physical
spectrum. The important point is that it is dan-
gerous to regard the lattice theory as phenomen-
ologically useful when the spacing a is larger than
either the scale under consideration or the char-
acteristic scales of the continuum theory.
If one adopts the renormalization scheme of
holding P(r, a,go(a)) fixed at a given scale r, then
at that scale there are by definition no O(a') cor-
rections in Eq. (2.2). However, in the following
we consider physics on different scales and so at
some point these corrections must plague us. In
the Monte Carlo analysis of Sec. III, our lattice of
104 sites forces us to keep r only a few times the
lattice spacing and thus these cutoff corrections
could be substantial for any a. Nonetheless we
proceed with optimism and neglect these terms.
'The ultimate test of this lies with larger lattices,
but our results in verifying asymptotic freedom
are encouraging.
Since P is dimensionless we can scale a factor
of 2 from r and u in Eg. (2.2) to give
P~ 2~, &,8, 2 ( ~=P(~i&iso(&)) i
f a)) (2.3)
where the approximate equality represents the
neglected finite-cutoff corrections. "Zhis equation
shows the correlation between the bare coupling
MONTE CARLO STUDY OF RENORMALIZATION IN LATTICE. . . 1817
for two values of cutoff and the measured observa-
ble at two different length scales. The process
that gave Eq. (2.3) is now iterated, with x re-
placed by 2 "x, to give the pivotal formula
((slj
f (aGl gol 2-~ l= P.~
(2.12)
(2.13)
( (& l ( 2 j (2.4)
To study the renormalization of g„we use this
formula as follows. Assume that for some fixed
values of x and a we use some technique such as
Monte Carlo si:mulation to calculate the two func-
tions of bare coupling
&(g.) =P(~, s,g.),
G(g, ) =P (2r, a,g,) .
(2.5)
(2.6)
lim P(r, a,g,(a)}=P,.
g ~p (2.7)
Then on a graph of E(g,) versus g„we find g,(&)
as the value of g, where
E(g,(a)}=PO. (2.8)
From Eq. (2.3) we find the coupling at half this
cutoff:
(2.9)
Knowing go(a/2), we define P, by
~.=&l ~ I-) I.
Equation (2.4) now tells us how to find go(a/4):
'lgl- l=P ~
Iterating gives
(2.10)
(2.11)
F, G
Suppose further that we pick x so that in the phys-
ical continuum limit P has the value P,:
This entire procedure generates graphically a
"staircase" as shown in Fig. 1. This figure is
drawn for an asymptotically free theory with
go(0) = 0.
In Fig. 2 we sketch a situation where the func-
tions E and G cross each other at a nonvanishing
g, . Here the staircase asymptotically approaches
this crossing point. We have a renormalization-
group fixed point g~ defined by
&(gz) = G(gz) . (2.14)
Note that g~ can be approached either from
stronger or weaker coupling. As the bare charge
at some very small cutoff passes through g~, the
corresponding initial value P, drastically changes
as we go from a staricase on one side of g~ to the
other. This is the signal of a phase transition in
the equivalent statistical mechanical system. 'The
critical exponents of the transition are related to
the relative slopes of E and G near the critical
point. The absolute slopes depend on the value of
a/r chosen to define these functions.
The above example represents a conventional
ultraviolet attractive fixed point. One could also
imagine a fixed point where again E(g~) = G(g~)
but
d,—&(g) —~
—'G(g)
dg '
dg (2.15)
F, G
In this situation the staircase construction will al-
ways lead away from g~. A continuum limit at
such an ultraviolet repulsive fixed point is at best
possible if go is exactly equal to g~. It is also
conceivable that at some point in the construction
Eq. (2.13) may have no solution. Several authors
0P
I
P
p ——————— I
2 /' I
I
I
'a ' o
0~—,~ ~0
I
I
I
I
I
I
go(a)
FIG. 1. The 'Staircase" construction of Sec. II for an
asymptotically free theory.
Po-
I
I
I
I
I
I
I
I
I
I
I
I
I
I
I
I g
QF go(a )
FIG. 2. An example of a nontrivial fixed point.
1818 MICHAEL CREUTZ
have argued that this may be the case for four-
dimensional P' theory, which thereby may not
have a nontrivial continuum limit. '
Non-Abelian gauge theories are asymptotically
free. This means that they have an attractive
fixed point at g~ = 0. -Perturbative expansion about
this point yields g(r, a, g)= g+O(g '). (2.23)
In particular, any change in definition satisfying
The quantity g(r) represents an "effective" coup-
ling at the scale r. There is considerable free-
dom in the definition of g(r); the only require-
ment is that as go is varied with the cutoff in place
w'e have
—g (a)=~(g )=~ g +~,g '+O{g '). {2.16 g(r) -g(r)+ O(g'(r)) (2.24)
ll 15
3 8v2il i
34 1)'|-3 8,s)
Integrating Eq. (2.16) gives
(2.17)
= P, Inj, , i+ —' ln Ini, , i
+ O(g,') .g2a 0 A Ram' P (A mani
(2.18)
Here A, is an integration constant which depends
on the cutoff scheme. In Ref. 2 we used Monte
Carlo methods for pure SU(2) and SU(3) gauge
theories with Wilson's lattice cutoff to calculate
Ap in units of the square root of the string tension
iC. For SU(2) we found
Although in general P(go) depends on renormalima-
tion scheme, the first two coefficients P, and P,
are determined if g, becomes the conventionally
normalized Yang-Mills charge in the classical
limit. For SU(2) gauge theory these coefficients
are"
is a. perfectly acceptable new renormalized
charge. Indeed, we are free to define
,( ) P(r)-Po (2.25)
IH. RENORMALIZATION AND WILSON LOOPS
We use the standard Wilson formulation of
gauge fields on a hypercubical lattice. ' An ele-
ment U;& of a unitary gauge group is associated
with each nearest-neighbor pair of lattice sites i
and j. We restrict ourselves in this paper to the
groups SU(2) and U(l). The action of the theory is
a sum over all elementary squares 0 of the lat-
tice
(3.1)
With this definition the perturbation series in Eq.
(2.22) for P{r) becomes trivial. Thus, with an ap-
propriate shifting and rescaling, we can interpret
our SU(2) results in terms of any effective coup-
ling at a variable scale r.
A, = (1.3+ 0.2) x 10 'v E'. (2.19)
Equation (2.18) implies that for a factor of 2
change in cutoff
So = Tr(U;)Uq~U„, U„)
for SU(2) or
(3.2)
1 1
g( /2) 2( )+ 2p ln2+ O(g ) (2.20) S~ =g Re(U;qU)~U„, U„) (3.3)
1 11
E(g )=P —,--;ln2go2 12m (2.21)
Verification of this behavior is a primary check
on our procedure.
When the physical scale x is small enough in an
asymptotically free theory, one expects the validi-
ty of an asymptotic perturbation series
lim P(r, a, g(a)) -=P(r)
a~o
p„g'"(r)+ O(g (r)) .
n-0
(2.22)
This gives an asymptotic freedom prediction for
the functions I and 6 used in the staircase con-
struction. For small g, we expect for SU(2) gauge
theory
for U(l). Here the indices i, j, k, and I circulate
about the "plaquette" Q. 'The quantum theory is
defined via the path integral
(3.4)
where every link variable is integrated over with
the invariant Haar measure for the group. The
U(l) theory is known" to exhibit a phase transi-
tion at gp=g~= 1.0.
'The most studied "order parameter" of lattice
gauge theory is the Wilson loop. Given a closed
contour C of links in the lattice, one constructs
the expectation value of the product of link varia-
bles about that contour:
(3.5)
MONTE CARLO STUDY OF RKNORMALIZATION IN LATTICE. ..
(A)= —jl, dU;qe~'~A(U). (3.6)
( i,g)
The Wilson loops are often used as a signal for
confinement. A linearly rising potential energy
between two widely separated sources in the fun-
damental representation of the gauge group cor-
responds to large loops having an exponential fall-
off of W(C) with the minimal surface area en-
closed by C. In contrast, a nonconfining theory
such as the weakly coupled U(1) model should give
a decay of the Wilson loop no faster than exponen-
tiaQy mith the perimeter.
We mould like to use the Wilson loops to con-
struct a physical. function for the analysis of Sec.
II. Unfortunately, the bare Wilson loop by itself
cannot be used because of ultraviolet divergences.
These divergences are of a rather trivial nature,
arising from the infinitely thin contour. They
represent the self-energy of the pointlike external
sources circumnavigating the loop. The problem
already appears in the exactly solvable continuum
theory of free photons. There the expectation val-
ue of the operator
explie ~ A„Ch„ l
c (3.V)
is easily shown to vanish. Inserting an ultraviolet
cutoff, one sees that this vanishing is exponential
in the perimeter of the loop measured in units of
the cutoff. If the contour has sharp corners, as
inevitable in the lattice theory, then powers of the
cutoff also occur in W(C).
We now assume that removing these perimeter
and corner divergences as well as appropriately
renormalizing the bare charge is all that is nec-
essary to render the Wilson loop finite in the con-
tinuum limit. This immediately implies that the
ratio of two Wilson loops of equal perimeter and
number of similar corners but with different
shapes will remain finite upon cutoff removal.
Such ratios, therefore, , can serve as the physical
quantities for the analysis of Sec. II.
Denote by W(I, J) a rectangular Wilson loop of
dimension I and J in lattice units. The above dis-
cussion suggests comparing ratios such as
W(2, 2}
W(1, 3) (3.8)
Here P denotes "path ordering"; the group ele-
ments are ordered as they are encountered in a
circuit of the contour. The expectation value is
taken in the sense of the measure in Eq. (3.4); for
an arbitrary function A of the link variables we de-
fine
W(4, 4)
W(2, 6) ' (3.9}
Note, however, that in this simplest case we are
already in Eq. (3.9) asking for loops with 6 units
on a side. As our largest lattice is only 10 sites
long, we prefer to use a more compact ratio con-
sisting of four loops; thus we define
W(2, 2)W(1, 1)
[W(2, 1)]'
W(4, 4)W(2, 2)
[W(4, 2)]'
(3.10)
(3.11)
P(~, s,g,)=p, g,'+O(g, ') O+~r'
where
, (3.13)-
p, =,[8 arctan 2+ 2 arctan~ —2v —4 In(r)]
1
= 0.066 079 788. .. (3.14)
for U(l} gauge theory, and
P, =,[8arctan 2+ 2arctan —, -2w —41n(&)]= 3 1
= 0.049 559 841. .. (3.15}
for the SU(2} model.
In the strong-coupling limit we have for U(l)
'JJ(IJ) = (,) 1+ 0,(—, (3.16)
1 ~ l)~W(I, J)=,~~ 1+0 —,l
&o &
for SU(2). This translates into
&(g.) = 1 — .+ 0(Z. '),1
2go
(3.1V)
(3.18)
'The constant is added so that these quantities
vanish as go goes to zero. In Sec. IV we pres-
ent Monte Carlo measurements of the quantities in
Egs. (3.10) and (3.11}.
In Eq. (3.10}we have effectively taken r = 2a in
the more general physical ratio
(r rl fr r&
Le
For small go this quantity can be expanded pertur-
batively. We have not carried this out for arbi-
trary cutoff, but for a'/r' small a straightforward
calculation gives
with the doubled scale ~(g.) = 1 —
16 .+ o(Z, ")16g (3.19)
1820 MICHAEL CREUTZ
for U(1), or I.O
&(g.) =1 ——.+ o(g. '),
&o
G(g, )=l-—,+O(g, ")1
80
(3.20)
(3.21)
0.9
0,8 — W(1, 2) +
u(l)
for SU(2). All interesting structure should occur
in the coupling regime which interpolates between
the behaviors of Eq. (3.13) and Eqs. (3.18)-(3.21).
IV. NUMERKAL RESULTS
0.7 — W(2, 2) ~
0.6—
05—
~s
+ ++~
I.O
l6
In previous publications we have discussed the
Monte Carlo algorithms used to bring our lattices
into equilibrium. "Here we always work on a hy-
percubical lattice of 104 sites and with periodic
boundary conditions. To hasten the approach to
equilibrium, the bare coupling is first given a
damped oscillation about the values where mea-
surements are finally taken. Once satisfied with
the equilibrium of the lattice, we make a sequence
of -5-10 Monte Carlo sweeps through the lattice
and measure Wilson loops after each. Any partic-
ular shape loop is averaged over all similar loops
in the entire lattice, that is, over all translations
and rotations. To estimate the statistical errors
in some quantity, we first calculate the standard
deviation of the mean over the sequence of sweeps.
To allow for correlations between successive lat-
tices, these errors are all increased by 40%.
This factor is estimated from a measurement of
0.4— 0
0.'3—
0.2—
0. I—
h
44
O h
0 0
I/ I 6gp
I
0.5
00 I.O
Qp
2
FIG. 4. Wilson loops for U(1) gauge theory.
I.5
l.o
the successive correlations between 1 x 1 loops
over a sequence of 10000 iterations on a 24 lat-
tice at g,-= 1.75 with the group SU(2). The validi-
ty of this error estimate was also confirmed on a
measurement of E(g,) for U(1) over a sequence of
1000 iterations on a 44 lattice at go'= 0.82. For
U(1) we keep g,' at least five percent away from
0.9
0.8—
0.7—
0.6—
0.5—
0.4—
0.3—
0.2—
O. I
W(I, 2)+
W(2, 2).
W(2,4)4
W(4, 4) ~
SU (2)
0.9-
0.8—
0.7—
0.6—
C9
0.5—
0.4—
0.3
0.2—
~ F(gp)
(go)
0 0 nl
I.O
0.00 2.0
Qo
FIG. 3. Wilson loops as a function of go for SU(2)
gauge theory.
3.0
O.I
00 I 2 2
Qp
FIG. G. The quantities I' and G for the SU(2) theory.
MONTE CARLO STUDY OF RENORMALIZATION IN, LATTICE. . . 1821
1.0
0.9—
0.8—
F (gp)
0 7 G(( gp — (2 ~ In2) )
0
oe6—
C5
" 0.5—
oo4—
0.5—
0.2—
O. l— ~h
~4
0 I ~ 2
gp
FIG. 6. Testing the prediction of asymptotic freedom.
the transition value to avoid severe critical cor-
relations. Because of this we draw no conclusion
on the critical exponents of the transition.
In Figs. 3 and 4 we show as functions of g,' the
measured values of those loops needed for the
evaluation of Eqs. (3.10) and (3.11). We also plot
their strong-coupling limits. 'The absolute errors,
which are too small to show on this graph, are not
strongly dependent on loop scale. However, the
small numerical values for physically large loops
make their errors relatively more important. In-
deed, the main source of statistical error in our
analysis comes from the 4 x 4 loops, which are
barely measurable for g, in the strong-coupling
regime. Note the approach to the strong-coupling
equations (3.16) and (3.17) for g,'& 2 for SU(2) or
g & 1 for U(l).
In Fig. 5 we show for SU(2) the quantities E and
G of Eqs. (3.10) and (3.11) as functions of g,'.
Note that E(g,) apparently always lies below G(g, )
suggesting a picture like that of Fig. 1. There is
no evidence for any fixed point other than atg, = 0.
The figure also displays the strong-coupling
forms of Eqs. (3.18) and (3.19). For g,'& 2 the
numerical results already approximate this limit-
ing behavior; consequently no further structure is
expected. Also note the approach to the weak-
coupling behavior of Eqs. (3.13)-(3.15) when g,'
becomes small.
In Fig. 6 we plot E(g,) and
1 ll
( g' 12m' (4.1)
Comparing physical loop ratios on different
length scales, we have Monte Carlo evidence that
SU(2) non-Abelian lattice gauge theory does not
exhibit a renormalization-group fixed point away
from vanishing coupling. 'This provides further
evidence that confinement and asymptotic freedom
coexist in the continuum limit of the theory.
In this analysis we have been able to push the
lattice spacing to extremely small values. At the
easily reached bare coupling g,'=1, the cutoff in
I oO j F( )
G(gp)
I-I/ I6g80
I - I/2 g ~
0
O.l—
0.0l0 0.5 I
I.O l.5
go
2
FIG. 7. The quantities J' and Q for the U(1) theory.
versus g, . Note the excellent agreement with the
asymptotic freedom prediction of Eq. (2.21) for
g, & 1.8. This prediction should only apply in the
weak-coupling regime. This agreement is rather
astonishing in the light of the neglected O(a'/r')
terms.
In Fig. 7 we show the results for J" and G with
the U(l) theory. In the weak-coupling phase of
this model the loops are strongly dominated by a
perimeter behavior. As this is canceled in E and
G the resulting signal is quite small. The data
suggest that when g,' is less than g~' is less than
g~' the function E lies below G. This effect, how-
ever, is close to the error bars and must be re-
garded with some caution in light of the strong
cancellations involved in calculating E and G
from the loops. %hen g, is substantially below
g& the functions E and G are equal within errors,
indicating a scale-invariant theory. This is evi-
dence for the massless photon dominating the be-
havior of these loop ratios.
V. DISCUSSION
1822 MICHAEL CREUTZ
units of the string tension can be calculated from
Eqs. (2.18) and (2.19). This gives
a= (Gx 10 ')A ~2.
Thus me have seen asymptotic freedom at cutoff
scales much smaller than a typical hadronic size
and where direct measurement of the string ten-
sion by Monte Carlo techniques mould be intracta-
ble.
Recent conjectures based on the analogy of the
four-dimensional U(l) gauge model with the two-
dimensional XY model of statistical mechanics
suggest a continuous line of fixed points for g,
less than a critical value. ' Our results suggest
that the functions I and G may truly cross each
other at g~. If this effect is real, we conjecture
that our analysis is indicating how to take a con-
tinuum limit in which the magnetic monopoles" of
the compactified U(l} lattice theory survive. A
continuum limit should also be possible holding g,
constant below the fixed point, but in this limit the
only surviving spectrum will be free photons. The
monopole mass will go to infinity with the inverse
of the lattice spacing. In contrast, allowing gp to
go to the fixed point under the construction of Sec.
II may keep the monopole mass finite.
In this model a nontrivial continuum limit may be
possible from either side of the fixed point.
Starting on the weak-coupling side should give
monopole electrodynamics (with a possibly vanish-
ing monopole moment). Alternatively, the strong-
coupling side should give a Higgs phase where
the monopoles have condensed into a "magnetic
superconductor. *' Such a phase exhibits electric
confinement as implicit in Wilson's strong-coup-
ling expansion.
Returning now to the SU(2) case, we can calcu-
late from E(go) a renormalized charge g(r = 2a)
using Eqs. (2.25), (3.13), and (3.15). In Fig. 8 we
plot the inverse renormalized charge at 2a versus
the bare charge. For small 1/g, ' we also indicate
the strong-coupling limit 1/g'= p, (l -g, ') '
+ O(g, ~). Note that at weak coupling, i.e., large
inverse charges, the graph approaches a straight
line. The unit slope of this line demonstrates that
the a'/r' corrections to Eg. (3.15) are remarka-
bly small. 'The intercept measures the ratio of
l, 5
I 0
0.5
'0 0.5 I .0
I /go
I 5
FIG. 8. The inverse renormalized charge squared
at w = 2a versus the inverse bare charge squared for
the SU{2) theory.
the asymptotic freedom scales for g and g„
= 2P, ln —)+ O(g, ) .2A)
Here A is defined in analogy to Ap..
(5.2)
or
2e. inl —l= o.35
&2A)
(Ao ) (5.4)
A = 22Ap. (5.5)
We do not put an error on this number because of
unknown systematic effects due to finite-cutoff
corrections. The large factor in Eq. (5.5) is typi-
cal of comparisons of the lattice Ap with more
physical definitions. " Using Eg. (2.19) for Ao we
obtain
A = O.3MA. (5.6)
Verification of the number in Eg. (5.5}is in prin-
ciple possible with a one-loop perturbative calcu-
lation.
= p, lnl, , l+ —1 l lnl, , l l+ O(g'(r)) .t'l1 6, f (1)yg2 g 0 (A2~2) p g2~2
(5.3)
Estimating the intercept from Fig. 8 we obtain
M. Creutz, L. Jacobs, and C. Rebbi, Phys. Rev. Lett.
42, 1390 (1979); Phys. Rev. D 20, 1915 (1979);
M. Creutz, Phys. Rev. Lett. 43, 553 {1979);C. Rebbi,
Phys. Rev. D 21, 3350 (1980); B. Lautrup and
M. Nauenberg, Phys. Lett. 958, 63 {1980); G. Mack
and E. Pietarinen, jbjd. 948, 397 {1980);D. Petcher
and D. H. Weingarten, Phys. Bev. D 22, 2465 (1980).
2M. Creutz, Phys. Rev. D 21, 2308 (1980); Phys. Rev.
Lett. 45, 313 {1980).
3K. G. Wilson, report, 1979 (unpublished).
K. G. Wilson, Phys. Rev. D 10, 2445 (1974).
R. Balian, J. M. Drouffe, and C. Itzykson, Phys. Rev.
D 10, 3376 (1974); 11, 2098 (1975); 11, 2104 (1975);
L. P. Kadanoff, Rev. Mod. Phys. 49, 267 (1977).
MONTE CARLO STUDY OF RENORMALIZATION IN LATTICE. . . 1823
A. H. Guth, Phys. Rev. D 21, 2291 (1980).
~C. Itzykson, M. E. Peskin, and J. B. Zuber, report,
1980 (unpublished); A. Hasenfratz, E. Hasenfratz, and
P. Hasenfratz, report, 1980 (unpublished); M. Luscher,
G. Munster, and P. Weisz, report, 1980 (unpublished).
R. H. Swendsen, Phys. Rev. Lett. 42, 859 (1979).
SK. G. Wilson and J. Kogut, Phys. Rep. 12C, 75 (1974);
G. A. Baker and J. M. Kincaid, Phys. Rev. Lett. 42,
1431 (1979).
0%'. E. Caswell, Phys. Rev. Lett. 33, 244 (1974);
D. R. T. Jones, Nucl. Phys. 875, 531 (1974).
T. Banks, R. Myerson, and K. Kogut, Nucl. Phys.
8129, 493 (1977); T. A. Degrand and D. Toussaint,
Phys. Rev. D 22, 2478 (1980).
~A. Hasenfratz and P. Hasenfratz, Phys. Lett. 938, 165
(1980).