# Validation scope

Version: 0.36.0-rc.1. Runtime cache revision: 25.

Observed execution on 2026-11-23 through 2026-11-27 used macOS 27.2.1 on arm64, TeX Live 2026,
LuaTeX/LuaHBTeX 1.24.0, PGF 3.1.12, PGFPlots 1.18.3 and
ConTeXt MkIV 2026.09.10 with the LuaTeX engine. These are observed versions,
not a claim that all other engines/platforms have been tested.

The numerical regression comprises 70 Lua suites covering analytic/manufactured
fields, nonlinear material laws, permanent magnets, meshing/topology, boundaries,
profiles/integrals, static/harmonic forces, harmonic diffusion/circuits and sheet losses.
A separate rational-arithmetic oracle checks geometric predicates. Invalid TeX
inputs and operation order are checked. Cold, warm and shared frozen calculations
agree across plain LuaTeX, LuaLaTeX and ConTeXt MkIV to the tested 1e-10 tolerance.
The isolated installation check copies runtime files into a private TEXMF,
compiles consumers outside the checkout, records loaded TeX/Lua paths and fresh
PDFs, and checks failure for missing siblings even with a complete shadow copy.

The manual uses focused PGF source/result examples with explicit shared setup.
Its plain and ConTeXt loading wrappers are the two native-format exceptions;
the LaTeX wrapper executes its shared scene inside the manual. Short loading
skeletons and configuration fragments execute in declared model/axis/format
contexts. A source inventory distinguishes those contexts from rendered
examples and verifies that no displayed PGF block is omitted or repeated.
Downloadable complete documents are checked independently. The README example,
index, cross-references and rendered figures are also checked; compilation
alone is not evidence of scientific accuracy.

Independent material-loss comparisons performed on 2026-10-22 used xfemm revision
`34ac6cc586a346066984526620dfb67383c8e3c5`: 15 samples across intrinsic complex
anisotropy, ordinary/deep-skin sheets, axisymmetry and a nonuniform core. The
largest relative differences were 0.255% for potential, 2.898% for B/H, 0.280%
for energy and 1.304% for magnetic loss. These comparisons apply to those cases,
not to all magnetic physics. They were not rerun for documentation/layout changes.

Independent harmonic-load checks use three independently meshed conductor
problems with phase-shifted excitation. Mean and complex twice-frequency Lorentz
and weighted-stress forces/torques differ by at most 0.146% of the respective
reference force/torque scales in this fixture. Air-contour checks at two offsets
show a maximum 0.395% difference from the volume force on the finest mesh;
these are sensitivity results, not general force-error bounds. Analytical
factor-of-two limits, direct time reconstruction, phase covariance, contour
orientation, origin changes and axial symmetry/sign checks complement them.
First-order operating-point loads and fixed-angle harmonic P1 bands have
additional checks described below.

Independent effective-wire comparisons use the same xfemm revision
and independently generated meshes for round, untransposed-strand and Litz
windings at 400 Hz and 100 kHz. Twelve air/winding samples differ by at most
0.207% for potential, 3.678% for unsmoothed B/H, 0.137% for energy,
0.043% for total loss and 0.137% for the complex terminal voltage norm.
The field comparison is mesh sensitive; it validates the shared approximation,
not arbitrary manufactured windings. Separate tests cover the DC limit,
phase covariance, voltage/parallel ports, axisymmetric power balance and caches.

Nonlinear fundamental checks integrate the original H(B) table independently
with 20,000 angular samples, test consistent tangents and exact uniform fields,
exercise a saturated current-driven solve, and check average energy/coenergy
and reactive power. This sinusoidal-B projection is distinct from FEMM's
B(H)-based effective large-signal approximation. These fundamental-law checks
do not cover transient hysteresis, higher harmonics or nonlinear AC forces.
The separate lossy equivalent-amplitude law is checked below.

The original public harmonic interchange checks used 18 samples in each direction: native
reading of exported Lua answers, independent solves of exported FEM problems,
and Lua reading of native answers. Same-mesh differences are below 1.5e-11
relative across potential, B/H, current, energy, total loss and terminal voltage. Independent
meshes differ by at most 0.352% in potential, 11.404% in unsmoothed local B/H,
0.0208% in current, 0.220% in energy, 0.000413% in total loss and 0.0344% in
complex terminal voltage. The local B/H
maximum occurs inside a split winding and reflects mesh sensitivity; it is not
a general accuracy bound. A native circular-arc answer also recovers an analytic
uniform field, and a native nonlinear DC answer supplies an unchanged-mesh
harmonic response. Unsupported records fail explicitly.

The general regression and optional native oracle harnesses are development
infrastructure and are not required to install or rebuild the public manual.
Users can reproduce the displayed examples with the included documentation
builder, compare physical limiting cases and refine meshes/domains for their
own problem. Windows FEMM GUI acceptance, all-platform support and full FEMM
parity are not certified. General nonlinear harmonic waveforms, history-dependent
hysteresis and resolved individual-strand currents remain outside the documented
supported subset. Square-wire and CCA homogenisation and the bounded lossy
equivalent-amplitude law have separate checks below.


Signed harmonic periodicity is checked against exact one-dimensional conducting
slab solutions with periodic and antiperiodic constraints. Maximum nodal
potential errors at the tested spacing are 5.49e-6 and 2.28e-6 Wb/m respectively.
The tests check both phasor parts, including the signed conductivity mass.

Moving-band checks cover the lossless half-static mean-torque limit, doubled
phase covariance, time reconstruction, rotating complex material axes and owned
bias/response caches. First-order torque on a nonlinear rotor agrees with a
central difference of complete DC solves to 2.1e-8 relative in the fixture.
Planar and axisymmetric Lorentz, weighted and contour force variations are
compared with nonlinear central differences; harmonic response phase checks
use lossless response models. Both bias-current and bias-field cross terms are
retained. These checks do not establish finite-frequency diffusion accuracy by
static differentiation alone.

Four additional independently meshed xfemm comparisons cover harmonic torque
on lossless and conducting reluctance machines at two mesh spacings. The
maximum error among mean, real twice-frequency and imaginary twice-frequency
coefficients is 1.415% of the reference coefficient scale. Errors are not
monotonic between those two meshes; this is a comparison of two different
meshes and stress weights, not a general convergence rate or error bound.
The native revision is the same `34ac6cc586a346066984526620dfb67383c8e3c5` used above.
No FEMM 5 executable was built for these comparisons.

Periodic-sector or analytical sliding gaps, motion-induced voltage, angular
stiffness derivatives and native incremental-answer export remain outside the
validated scope. Named PGF styles, force-component scoping, signed harmonic
frontiers and DC/AC position/response workflows are exercised by the generic
format fixtures and copied-installation consumers.

### Analytical gaps and planar previous-solution interchange

The Fourier annulus tests use exact radial/logarithmic modes and the known
angular P1 projection, including analytic modal energy and Maxwell torque.
A rotating reluctance example checks torque against a central coenergy
difference: 0.0002407296377 versus 0.0002407288869 N m. This is a discrete
virtual-work check, not a bound on bulk mesh error.

The unmodified FEMM 5 CLI solver, source revision
`31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3`, independently reads the exported
DC companion mesh and solves three planar nonlinear-bias responses. Relative
maximum nodal-potential discrepancies are 1.60e-10 (incremental DC), 2.91e-11
(frozen DC) and 7.84e-11 (50 Hz incremental). Native outputs are also restored
passively. This establishes a bounded common-mesh interchange comparison,
not general compatibility with every native answer dialect.

Constrained-mesher spatial searches are checked against a separate exhaustive
candidate oracle on 32 scaled, translated and contour-reversed variants of
eight original geometries. Node coordinates, ordered triangle connectivity,
material assignments and protected segments agree exactly; affine field patches
also pass. Independent box/disk tests cover index mutation, overflow, tangency
and deterministic ordering. The corpus removes 99.48% of broad-phase candidate
visits; this is not the same as a proportional reduction in total solve time.
Three alternating CPU-time runs against the retained previous Lua runtime give
median mesh speedups of about 1.5 for nested contours and 1.7 for many slots
on the stated host, with smaller gains on simpler inputs. Retained Lua memory
was measured after collection; peak process memory was not measured.

Original topology regressions inspired by FEMM 5 test both partners at shared
periodic corners, declaration-order independence, and P1/Fourier ring references
after merging near-duplicate interface vertices. The periodic manufactured
solution's maximum nodal error decreases from 0.0796 to 0.0174 when spacing is
halved. The merged and unmerged ring solutions agree to the tested tolerance.
Strict angle requests below a constrained-corner angle remain unchanged;
requests exceeding that angle now fail before consuming the insertion budget.
This is a necessary feasibility check, not an acute-corner protection algorithm.


The shared linear core has separate algebraic regression tests for extreme
uniform scales (1e-200 and 1e200), tiny nonzero loads, exact zero sources, invalid
and nonfinite input, singular/indefinite systems, a positive-definite operator
whose zero-fill Cholesky factor breaks down, bounded recovery and exhausted
iteration budgets. Accepted solutions are checked using the original operator
and returned physical unknowns. Row-wise normwise backward error supplements
the relative residual; neither is a field-error estimate.

Eight finite-element workloads at two resolutions were compared across five
numerical policies, with three alternating repetitions (240 runs). Meshes were
identical; maximum relative nodal-potential and centroid-field differences were
5.12e-9 and 3.87e-8. Energy, loss and band mechanics meet the tested 2e-6 comparison
gate; 24 closed-circuit linkage comparisons differ by at most 2.50e-15. Packed
products give modest gains over checked row tables. Ordering and additional fill
have mixed effects, so mesh order and zero-fill factors remain the defaults.
The additional reliability checks increase solve CPU relative to the previous
core; these results do not promise a general speedup or a bound for every model.

Independent FEMM 5 exchange checks use unmodified solver revision
`31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3` and xfemm's independent answer reader.
There are 54 comparisons across exported answers, remeshed solves and reimported
answers. On the same mesh, potential/B/H differences are below 1.1e-13 relative;
current differs by at most 2.23e-9, reflecting a nearly cancelling native voltage
gradient. Independent meshes differ by at most 0.280% for potential, 8.573% for
unsmoothed local B/H, 0.0163% for current, 0.0316% for energy, 0.0000297% for
loss and 0.00567% for terminal voltage. Local derivatives remain mesh sensitive.
Three same-mesh previous-solution responses have relative potential differences
below 7.85e-11. The ordinary exports use conventional property delimiters and
carry global/local mesh sizes into background air as well as material faces.
These checks establish only the documented subset; they do not certify the
Windows FEMM GUI or enable transient/AGE exchange.


The same 54 comparisons were rerun with xfemm's pinned legacy backend after the
export corrections. Same-mesh discrepancies remain below 1.18e-11 relative.
Independent-mesh maxima are 0.424% for potential, 11.408% for local B/H,
0.0136% for current, 0.0495% for energy, 0.000077% for total loss and 0.00777%
for terminal voltage. The potential/current/energy differences are small, while
local derivatives near interfaces need their own mesh study. The two native
solvers use different meshers; these are separate comparisons against Lua,
not an assertion that their remeshed arrays should coincide.

Path-based FEM/ANS exchange rejects same-basename `.tmp` companions without
reading or modifying them. Focused tests cover destination and owned-bias
preflight, transient headers, unsupported restarts, extra nodal previous-state
columns and native AGE counts. In-memory text parsing has no sibling context.
The native comparison fixtures contain no active companions.

Numerical-policy tests exercise explicit CG/GMRES, diagonal/incomplete factors,
RCM ordering, residual acceptance, failed budgets, cache invalidation with mesh
reuse and derived-model overrides. Voltage-driven complex systems agree with
the default policy on the same mesh. Generic-format consumers exercise PGF
styles, appended styles, grouping, FEM overrides and operating/sliding models.
The checked inventory contains 662 public/helper declarations. Solver tolerances
and iteration controls are algebraic controls, not discretisation-error bounds.


The default material collection is pinned to FEMM 5 revision
`31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3`: 282 records, with an explicit MIT
notice. All records are checked, alongside the optional 246-record legacy input.
Nonlinear directional DC sheets are tested against constituent interface
conditions, finite-difference tangents/energy, linear and unit-fill limits,
explicitly meshed iron/air layers, and planar/axisymmetric frozen/incremental
responses. Axis-aligned DC FEM/ANS records preserve the curve and fill; oblique
sheet and nonlinear-sheet response exports are rejected.

Four independently meshed sheet-inclusion solves against the unmodified FEMM 5
CLI use two sheet directions and spacings 0.12 and 0.06 m on a 2 m domain.
Twelve local B samples differ by at most 5.47%; at the finer mesh the central
field differences are 0.139% and 0.088%. These are finite-mesh comparisons,
not error bounds. The same files are read and solved by xfemm, whose older
mixing law differs by up to 12.1% at these samples. File compatibility therefore
does not assert equality of FEMM 4 and FEMM 5 nonlinear sheet physics.


Electrical network checks compare DC and harmonic ports with independent
current-driven impedance solutions, including mutual coupling between two FE
windings, parallel internal branches, analytic RLC limits, KCL and absorbed
active/reactive power balance. Nonlinear DC and linear harmonic linkage targets
recover independently prescribed currents; axisymmetric linkage and frozen
circuit-state restoration are included. Incompatible ideal-source constraints
and disconnected electrical nodes fail explicitly.

Reduced-sector periodic/antiperiodic gaps are checked against annular Fourier
modes and the coenergy derivative of torque. Harmonic peak/mean factors and
refinement sensitivity are tested. A separate native FEMM 5 answer containing
10,502 nodes, AGE quad records and raw rings is preserved numerically and read by
xfemm, both with and without the optional raw-ring trailer. This is a file-
preservation check; it does not validate a native AGE field operator in Lua.
The separate operator and eccentric-position validations below cover those
additional capabilities.


## Native air-gap operators and eccentric positions

The runtime independently factors the native cubic/quintic stiffness through
Hermite endpoint data and a four-corner rectangular energy. On 120 combinations
of fractional shift and aspect ratio per basis, comparison with native expanded
C++ arithmetic gives maximum scaled entry errors of 4.67e-14 (cubic) and
6.67e-14 (quintic). The pinned FEMM 5 source is
31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3. The separate xfemm postprocessor
reconstructs the exported cubic gap field and potential to 2.61e-14 relative
error at the tested angles. An unmodified FEMM 5 same-mesh solve of the exported
linear reluctance problem agrees to 2.89e-12 relative maximum nodal error.
The native Windows GUIs were not executed.

Native static and harmonic snapshots, native stiffness energy, optional raw
rings, explicit interpolation selection and frozen caches are checked. The
native formulation is a thin-gap approximation; midgap stress torque and the
coenergy derivative differ by about 2.2 percent in the tested relatively thick
annulus. This is documented rather than hidden by claiming an exact discrete
virtual-work identity.

Eccentric P1 positions retain bulk connectivity and pass affine-field tests at
several rotor angles, including near a complete revolution. Collision and
polygonal nesting failures are rejected. For a nonlinear current-fed reluctance
rotor at offsets (0.35,0.2) mm, the weighted band force component Fx agrees with
a central coenergy derivative within 1.93e-4 relative error on the coarse mesh and
2.94e-5 on the finer mesh and smaller step. The tested mesh spacings are
0.8 and 0.4 mm, with displacement steps 0.02 and 0.01 mm. These are sensitivity
checks for a specified fixture, not universal error bounds.

Reduced-sector views transform bulk P1 values, rotate vectors and apply signed
periodic continuation before clipping. Regression checks exercise a rotor
crossing the radial seam, static and harmonic sampling, positive display
triangles and unchanged calculation signatures. Polar gap rendering is separate
from the computation operator. The generic consumers exercise native answer
export/import, eccentric positions and nonzero-angle sector fields in plain
LuaTeX, LuaLaTeX and ConTeXt MkIV, including shared frozen caches.

A separate ordinary FEMM 5 remesh of the exported problem also passes. The
sampled gap-field difference is 2.34 percent for these different meshes. Native
quintic cell-centre field expressions agree to 9.66e-14 relative error.
Native FEMM 5 periodic and antiperiodic quarter-sector answers at 13 degrees
(162 nodes each) also import successfully; signed continuation of the gap field
and positive native stiffness energy are checked. Regional volume averages are
not claimed to survive AGE interchange: native answers omit gap volume triangles.


## Transient fields and material-coordinate motion

Backward Euler is checked against a planar analytical diffusion mode, exact
linear discrete energy balance and first-order time refinement. On the tested
fixed mesh, errors 0.01235384 and 0.00592465 against a finer-time reference give
ratio 2.085 when the step halves. Axisymmetric linear/nonlinear cases dissipate
energy; isolated planar conductors satisfy their integrated current constraint.
Affine boundary scaling persists across steps, and a nonzero initial clock
survives state capture/restoration. A stationary initial field stays stationary.

Restart continuation agrees with uninterrupted calculation. Tests reject changed
state data, altered mesh coordinates, invalid time steps and incompatible frozen
caches. Full-circle P1 tests preserve conducting bulk connectivity and previous
material-node potentials; lossless imposed motion agrees with an independent
stationary position. Current-fed winding emf equals the linkage difference over
the step, and the terminal relation includes resistance. The partitioned torque
update is checked directly. With viscous damping and zero magnetic torque,
angular-speed errors 0.01596462, 0.00815389 and 0.00412146 show first-order
convergence towards the exponential solution as the step halves.

Independent, unmodified FEMM 5 solver revision
`31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3` was run on the same linear planar mesh
at steps 0.00005, 0.0002 and 0.001 seconds. Relative nodal potential differences
are 2.0132e-12, 2.3558e-13 and 1.7483e-12. Imported induced-current samples also
agree. This checks the bounded real planar time-answer contract, including
separate previous-potential records and an explicit absolute import clock.
It does not validate transient AGE, axisymmetric native time exchange or the
Windows GUI. Native solver controls are private comparison inputs; runtime
imports neither execute nor consume those sidecars.

Three source/result examples show a switched conductor, reluctance-rotor motion
and native time-answer interchange. They were compiled and inspected. The generic
consumers exercise state capture, stepping, PGF expressions/styles, time plots,
passive restart, dynamic motion and native time import in all three formats.
Time-domain electrical networks, arbitrary conducting remeshing, irreversible
hysteresis and moving native/Fourier gaps remain unsupported.


The 46 Lua suites, catalogue checks, exact-predicate oracle, existing TeX
regressions and invalid-input checks pass. The extended time fixtures also pass
cold, warm and shared frozen caches across all three formats to 1e-10, with
no solver execution on a cache hit. Copied isolated installations pass with
TeX/Lua provenance in plain LuaTeX, LuaLaTeX and ConTeXt MkIV.

A separate prescribed-translation check uses a nonmagnetic conductor with
conductivity 1000 S/m in an imposed nearly uniform field of -1 T along y.
Translation at 0.01 m/s along x gives induced current density -9.99999950 A/m2,
consistent with the weak-reaction limit `sigma (v cross B)_z = -10 A/m2`.
No zero-total-current constraint is imposed in this electric-field reference.

The static circular multipole exterior is explicitly rejected for time states:
its existing zero-net-current audit includes impressed current only. An induced
return-current/exterior compatibility audit is required before extending this
approximation to time diffusion. Finite exterior domains and explicit boundary
data remain available. The diagnostic has a retained regression.

The time-domain documentation snapshot had 267 pages, 286 executed PGF blocks, 57 worked cases and
71 executable reference fragments. References and index are stable; strict
checks find no overfull boxes or missing glyphs. The new figures, key reference,
Lua API descriptions and index were inspected. All 769 public declarations have
reviewed manual references.


## Analytical coplanar electric charges (2026-11-23)

The new electric provider is checked against the three-dimensional Coulomb law,
including signed superposition, inverse-distance potential, inverse-square
field, homogeneous permittivity and coordinate-unit scaling. Independent centred
finite differences check E = -grad V; spherical surface quadrature checks
Gauss's law. These are genuine point charges in C, not planar line charges.

An exact unit-radius equipotential provides a rendering-convergence test:
maximum radius errors at horizontal resolutions 20, 40 and 80 are
0.0416326531, 0.0110609104 and 0.00277595704 coordinate units. Contour refinement
does not change the analytical samples. Tracing checks radial lines, decreasing
potential along directed dipole lines, orientation after backward integration,
null points, source/window termination and normalized vector lengths.

Owned snapshots, source cancellation, raw singular evaluations, masked profiles,
and gaps narrower than the profile spacing are checked. Twelve invalid PGF cases
exercise unsupported values, missing data, declaration order, mixed physics and
singularity diagnostics. Transformed nodes, normal anchor queries, local styles,
curved profile/seed paths, all view kinds and PGFPlots consumers execute in plain
LuaTeX, LuaLaTeX and ConTeXt MkIV. Copied private installation and TeX/Lua
provenance checks pass with the electric modules included.

The complete 47-suite regression passes, together with catalogue checks,
independent geometric predicates, TeX input/scoping checks and cold/warm/shared
frozen-cache comparisons in all three formats to 1e-10. Electric analytical
calculations do not use a disk solution cache. No electrostatic FEM comparison,
induced-conductor calculation or native electric interchange is claimed.


The rebuilt manual has 280 pages, 300 inline PGF example blocks and four
inventoried source-file listings, 58 worked cases and 71 executed reference
fragments/loading skeletons. The README example and named format exceptions
also pass. References and index are stable; no overfull boxes or missing glyphs
remain. The electric tutorial, plots, systematic reference and index were
inspected in the completed PDF. All 855 interface declarations have reviewed
manual coverage, and the final build fingerprints match the source tree.


## Classical electrostatic symmetry models (2026-11-23)

The seven homogeneous-medium configurations are checked against Gauss flux,
independent finite differences of potential, centre/interior limits, potential
continuity and the displacement jump at prescribed charged surfaces. Capacitor
energy is independently integrated from epsilon E^2/2 with the appropriate
volume, per-length or per-area measure; the relative gate is 1e-7. Signed and
zero voltages, charge/capacitance/energy normalization, finite cylinder potential
references and owned snapshots are checked. Equivalent centimetre/millimetre
models and permittivity scaling at fixed charge versus fixed voltage agree to
1e-13 in the tested cases.

Profiles preserve both field limits and a continuous potential even when a gap
lies between two ordinary samples. Reversed paths, tangencies, endpoints and
zero-source absence of artificial jumps are covered. Direct thin-gap field and
potential checks extend to 1e-10 relative electrode spacing. Exact radial/normal
lines preserve physical orientation and the shared null-field threshold.
These are bounded floating-point tests, not a guarantee for arbitrarily close
or extremely large coordinates.

All seven model choices, model-derived geometry, common views, electric plots,
owned snapshots and formatted integral quantities execute in plain LuaTeX,
LuaLaTeX and ConTeXt MkIV. The copied installation passes in all three formats
with runtime provenance and incomplete-installation guards. Nineteen public
PGF input errors are rejected. The two new runnable examples compile; the
focused source/result tutorial and systematic reference were inspected visually.
The interface inventory now covers 960 declarations.

These analytical providers do not solve finite ends, fringe fields, dielectric
interfaces, floating conductors or arbitrary conductor geometries. No native
electrostatic FEMM solver comparison or electric file interchange is claimed.

The 49 Lua suites, catalogue/predicate checks, direct TikZ and invalid-input
checks pass. After separating a paragraph of figures in the plain-format test
fixture, the complete cold/warm/shared-frozen matrix passes again in all three
formats to 1e-10. These magnetic cache comparisons guard integration regressions;
the analytical electric models themselves have no disk solution cache.


The complete manual now has 309 pages, with 333 inline PGF blocks, four
inventoried source-file listings, 63 worked cases and 71 executed reference
fragments/loading skeletons. Its references and index are stable, with no
overfull boxes or missing glyphs. The charged-body/capacitor tutorial, reference,
profiles and new index entries were visually reviewed in the final PDF.
Manual input and isolated-installation runtime fingerprints match the checkout.


## Planar electrostatic finite elements

Scalar planar Poisson tests independently check exact two-dielectric slabs,
internal sheet-charge jumps, electrode reactions, nonnegative field energy,
outward displacement signs, physical units and extrusion depth. Compound
holes, independent snapshots and contradictory boundary inputs are covered.
A manufactured quadratic potential gives RMS errors .024879, .00605683 and
.0011951 for meshes .4, .2 and .1; corresponding raw-gradient errors are
.214882, .105723 and .0530497. These establish convergence for that fixture.

A coaxial capacitor refines both its polygonal boundaries and triangular mesh.
Relative capacitance errors against 2*pi*epsilon/log(b/a) are .575971%,
.173290% and .047374%. The source-free solution also satisfies 2W=QU. Every
extracted P1 contour vertex has its requested potential; field lines stop at
conductors and vectors do not bridge dielectric or charge-layer jumps.
Finite-strip tests distinguish mesh refinement from exterior enlargement,
including their different influence on local fringe fields and electrode charge.
No universal pointwise accuracy or native FEMM electrostatic comparison is claimed.

The three-format fixture executes regions, electrodes, charge layers, compound
fills, curves, raw views/profiles, quantities and owned snapshots. It checks
rotation/scaling and native physical pics using TikZ's transform-shape option.
The isolated installation passes with loaded-file provenance and matching runtime
fingerprints. The full regression has 49 Lua suites; 29 electric invalid-input
cases fail explicitly. The interface inventory covers 960 declarations. The
finite-plate, dielectric/sheet and charged-cylinder examples compile, and the
focused tutorial/reference figures were inspected alongside their printed code.

The planar scalar provider uses isotropic linear media and a finite rectangular
exterior. Its axisymmetric extension is validated below. Floating conductors, capacitance matrices,
electric forces and native electric interchange are not included. Raw P1 gradients
remain discontinuous; sharp electrode-end maxima and artificial-domain truncation
require separate sensitivity checks.

The complete source/result manual builds with stable references and index,
including the planar electrostatic tutorial and systematic reference. It has
309 pages, 333 inline PGF blocks, four executed source-file listings,
63 worked cases and 71 checked reference fragments/loading skeletons. The
finite-plate comparison uses explicit seeds within its cropped wide-domain
view. The final tutorial, field/profile figures, key reference and index were
visually inspected. Runtime and documentation input fingerprints match the
copied installation report and successful manual build.

## Axisymmetric electrostatics and dielectric data

The scalar axisymmetric solver integrates stiffness with 2*pi*r and integrates
constant volume/sheet sources with exact weighted nodal moments. Independent
checks cover the uniform axial field, a charged disk, prescribed outward
normal displacement, m/cm scaling, full electrode charge, W=QU/2 and Gauss
balance. The axis retains a varying potential; magnetic A_phi axis constraints
are not used. Negative radial domains, explicit depth and incompatible axis
conditions fail explicitly.

With h=.4, .2 and .1, the axial quadratic-potential fixture gives RMS potential
errors .0236931, .00459360 and .00158386. The radial quadratic-potential fixture
gives weighted RMS field errors .126839, .0591895 and .0296416. The independently
refined spherical capacitor has 274, 934 and 3256 triangles and relative
capacitance errors 2.34650%, .688499% and .181192%. Volume- and surface-charged
sphere exterior-field RMS errors decrease from .0786075 to .0111109 and from
.0784675 to .0110734, respectively. A concentric grounded shell separates these
spherical comparisons from rectangular-box truncation error. These figures
are fixture-specific convergence evidence, not universal error bounds or a
same-mesh native FEMM equivalence claim.

All 95 dielectric records are independently compared with the pinned original
statlib.dat for names and numeric ex/ey/qv fields. Generation is byte-reproducible.
The library has 94 isotropic records and one anisotropic record; the latter
retains both original permittivities. The tensor FE provider accepts all records;
analytical providers continue to reject the anisotropic catalogue entry.
The MIT notice and upstream material-library terms are retained with the data.
The scalar adapter and solver are original LPPL code.

A copied isolated installation passes actual LuaLaTeX, plain LuaTeX and ConTeXt
MkIV consumers with TeX/Lua provenance and matching runtime SHA-256 checks.
Those consumers exercise PTFE/glass styles, grouping, symmetry ownership,
radial/axial plot aliases, discontinuous profiles and full-body charge units.
The spherical-capacitor, dielectric-disk and charged-body examples compile;
their source/result tutorial, plots and complete catalogue table were visually
inspected. Raw P1 gradients remain unsmoothed, including near the axis.
Floating conductors, dielectric tensors and electric forces are covered by
the extended validation below. Native electric interchange remains unsupported.

The 0.26.0-dev regression passed 51 Lua suites, 34 electric invalid-input cases,
generic TikZ consumers and cold/warm/shared-frozen cache comparisons in all
three formats. The interface review covers 974 declarations and fingerprints
the generated material-choice sources. The original magnetic catalogue also
passes its byte-for-byte regeneration check after the dielectric addition.

The 0.26.0-dev 323-page manual executed 347 inline PGF blocks and four source-file
listings, covering 66 worked cases and 71 checked reference/loading fragments.
Its references and index are stable, with no overfull boxes or missing glyphs.
The spherical field/profile figures, axisymmetric reference, corrected catalogue
labels, complete data table and index were visually inspected. Recorded manual
input fingerprints and copied-installation runtime hashes match the final sources.

## Extended electrostatic models

The 0.27.0-dev regression includes 52 Lua suites and 41 invalid electric PGF
cases. An isolated copied runtime passes LuaLaTeX, plain LuaTeX and ConTeXt MkIV
with matching SHA-256 provenance, including the advanced electric fixture.
Cold, warm and shared frozen magnetic results continue to agree across the
three formats to 1e-10. The interface review covers 1033 declarations.

Independent identities verify neutral/charged equipotentials, separated pieces
sharing one total charge, planar depth and full-revolution normalization,
Maxwell-matrix reciprocity and energy, independent-source offsets, owned matrix
states and Robin boundary storage. A rotated dielectric tensor is checked
against a manufactured uniform field with independently prescribed displacement:
E and D are not assumed parallel. Periodic tests synchronize unequal boundary
subdivisions, remove the constant mode by antiperiodicity, reject incompatible
corner voltages and include signed periodic reactions in charge balance.

Raw Maxwell-contour forces are checked against zero force in a uniform field
and Q times the external uniform bias for a symmetric charged sheet. Relative
errors for mesh sizes .2, .1, .05, .025 and .0125 are 4.346%, 2.400%, 2.722%,
2.477% and .132%; intermediate convergence is not monotonic. At h=.05, three
contour offsets give 3.077%, 2.909% and 2.474% errors. A charged sphere in a
uniform axial bias has 1.062% relative full-revolution force error at h=.07.
Contour orientation, contact rejection, planar torque, axis closure and
anisotropic contour rejection are also checked. These are fixture-specific
results, not universal force-accuracy bounds or claims of smoothed fields.

FEMM 4 and FEMM 5 source/manual inspection confirms that native electrostatic
point-property charge is in C/m: a transverse line in planar geometry or a ring
weighted by 2*pi*r in axisymmetry. A zero native point density selects prescribed
potential. The analytical Coulomb node has different semantics and remains
rejected by FE. This is a source audit; no native electrostatic solver execution
or numerical comparison is claimed. Native electric file exchange, singular
nodal adapters, general rotational periodicity, nonlinear/lossy dielectrics and
arbitrary 3D FE remain outside this implementation.

Additional regression cases check rotated/scaled axisymmetric force contours in
all three formats, automatic electrode-name collisions, complete tensor-interface
comparisons at mesh vertices, and absence of artificial gaps when a region repeats
the background tensor. The physical boundary tolerance preserves regular-axis
closures after PGF fixed-point coordinate transformations.

The complete 0.27.0-dev manual was rebuilt on 2026-11-24: 332 pages, 353
inline PGF examples and four executed source listings, covering 67 worked
cases and 71 checked reference/loading fragments. References and the index
are stable, with no overfull boxes, missing glyphs or unresolved references.
All 222 recorded manual input fingerprints and 95 copied-installation runtime
fingerprints match the final sources. The floating-screen tutorial, matrix
plot, tensor and mixed-boundary examples, Maxwell-force figure and index
were visually inspected, including source/result pagination.

## Native electrostatic interchange

Version 0.28.0-dev adds bounded electrostatic format-1 problem/answer exchange.
Passive readers reject unsupported records and unrepresentable physical data;
no imported text is executed. Tests cover nonrectangular labelled cells, holes,
circular arcs, floating conductors, diagonal tensors, independent volume/sheet
sources, translated periodic problem input and owned raw-solution restoration.
Unknown references, duplicate answer elements, conflicting physical overrides,
unsafe source replacement and unsupported native exports fail explicitly.

Native line-source refinement in a grounded polygonal disk gives off-source
potential errors of 0.56264%, 0.09947% and 0.05856% at mesh sizes 0.25, 0.12
and 0.06 m. The source is at the centre and the probe at (0.5, 0) m. Full-ring
charge, zero radial source weight on the axis and prescribed nodal voltage
are checked separately. These source idealizations do not give a finite
mesh-independent potential or self-energy at the source.

Independent xfemm and FEMM 5 CLI solves were executed on planar, axisymmetric,
internal charged-sheet and Robin fixtures. The xfemm postprocessor reads both
LuaFEMM answers and native answers without smoothing; LuaFEMM also reads the
native answers. Maximum absolute differences across potential/field components
in this corpus are 2.774e-8 for xfemm and 3.215e-9 for FEMM 5. These numbers
are fixture-specific, not a general interoperability or accuracy certificate.
The FEMM 5 electrostatic CLI was built and run on macOS; the GUI was not tested.

The compared FEMM 5 revision is `31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3`;
the FEMM 4 source audit uses `7d9e8ed90772ddb939a3256f976e231bcb5cca97`.
The native comparison found distinct internal sheet-charge assembly conventions:
FEMM 4/xfemm assigns the source to one adjacent element, this FEMM 5 loader to
both. The explicit edition setting preserves the total displacement jump.
FEMM 5 is the default; format 1 alone does not identify an edition.

The generic PGF fixture exports both file types, restores the answer, recomputes
the problem, captures a potential profile and renders it with PGFPlots. A final
copied installation passes plain LuaTeX, LuaLaTeX and actual ConTeXt MkIV;
all 99 runtime file hashes match the tested sources. The invalid electric PGF
corpus contains 48 cases and the interface review covers 1049 declarations.
The new five-page reference was compiled and its figures visually inspected.
Periodic native answer import/export, rotational interfaces, off-diagonal native
permittivity, nonlinear/lossy laws, exterior transforms and arbitrary 3D FE are
outside this bounded exchange implementation.

The final 336-page manual executes 357 inline PGF blocks and four source-file
listings (361 markers), covering 68 worked cases and 71 separately checked
reference/loading fragments. The complete exchange chapter and final index were
visually inspected. References/index are stable, without overfull boxes, missing
glyphs or unresolved references; all 228 manual input fingerprints match the
final sources. The 53-suite regression and all three generic cache modes passed.


## Image-charge fields and induced surface charge — 2026-11-26

The exact grounded-half-space construction uses a real Coulomb source and its
opposite image, with the analytical view restricted to the physical half-space.
The executed manual example captures the normal displacement along the plane
and compares the induced surface density with the independent image formula.
The corresponding generic consumers check the zero plane potential and the
normal displacement using the actual captured source coordinates, accounting
separately for TeX coordinate rounding after picture transformations.

Validation passes in plain LuaTeX, LuaLaTeX and LuaTeX-based ConTeXt MkIV on
macOS arm64 with TeX Live 2026 (LuaTeX/LuaHBTeX 1.24.0). The copied installation
checks all 99 runtime files and their provenance. The complete regression run
passes 54 Lua suites, the 48 electric invalid-input cases, and three-format
cold/warm/shared-frozen cache agreement to 1e-10. The public inventory remains
1049 declarations; the image construction uses existing keys and commands.

The rebuilt 339-page manual has stable references/index and executes 359 inline
PGF blocks plus four source listings, covering 69 worked cases and 71 separately
checked fragments. Its 230 recorded input hashes match the final sources.
The new field diagram, induced-charge plot and displayed source were visually
reviewed on pages 219--221 of that build. The subsequent source-correction
provider is validated separately below.


## Singular electric sources with a regular correction — 2026-11-26

The pure-Lua source-correction provider supports planar line charges and axial
point charges, scalar dielectric contrasts and prescribed-voltage conductors.
Independent image solutions test both source kinds beside grounded planes and
dielectric interfaces at three mesh scales. Potential-error norms decrease on
every refinement. At the finest spacing, relative vector-field errors at the
three fixed probes are 1.089%, 1.679%, 0.565% and 0.496%, respectively. These are
case-specific raw-field observations, not global error bounds; point-sampled
raw-gradient errors can be nonmonotonic on intermediate meshes.

Off-centre sources inside a grounded circular cylinder or spherical shell are
compared with independent image formulae. Refining geometry and mesh from 16
to 64 boundary subdivisions reduces the three-probe potential-error norm from
.02513 to .001469 V (planar) and .001474 to .00009856 V (axial), with normalized
source coefficients. A neutral-dipole experiment varies the outer domain at
fixed mesh spacing. A dielectric-rectangle case is compared with an independent
ordinary-Poisson finite charged-disk regularization. The remaining differences
include discretization and finite source-geometry effects.

Other checks cover SI scaling, logarithmic-reference invariance, exact kernels,
finite-difference gradients, source cancellation, solid/thin electrodes,
homogeneous guards, nonlinear contour roots, field-line orientation and
termination, profile exclusions, owned snapshots, and adaptive quadrature
refinement/exhaustion. Additional response checks below validate induced
electrode charge, matrices, interaction energy and forces separately.
Divergent total self-energy and native hybrid export remain rejected. There is no arbitrary 3D
correction, off-axis ring source, or hybrid disk-restart claim.

## Singular-source boundary and mechanical responses — 2026-11-26

Off-centre sources in a grounded circular cylinder and sphere validate induced
charge, finite source interaction energy and total Maxwell stress against
independent image formulae. At offset 0.3 of the radius, normalized planar
forces converge from .33735911 to .33202250 to .32999424 for 16/32/64 boundary
subdivisions; the exact value is .32967033. The corresponding axial sequence is
.35804894, .36180248, .36253009, versus .36227509. Mesh size decreases with
boundary refinement. Charge and interaction-energy errors also decrease.

Force reversal/mask invariance is tested to 2e-8 relative. Moving the contour
agrees within 3%; fixed-mesh position differences of the independently defined
interaction energy agree with stress within 6%. These are case-specific
finite-element tolerances, not universal accuracy guarantees. Applied-voltage
charging work has factor one; induced response has factor one half. Changing
the planar logarithmic reference produces the independently derived constant
shift of self-subtracted energy without changing total potential.

Floating thin and solid conductors satisfy total-charge constraints and total
nodal equipotential relations. Maxwell matrices preserve reciprocity and
reconstruct operating charges from source offsets. A finite charged disk in
ordinary Poisson elements beside a rotated tensor gives decreasing probe
norms .00110864 and .000613403 relative to a hybrid reference; finite radius
and mesh errors remain. Tensor displacement components satisfy the constitutive
law. Homogeneous isotropic source guards remain required.

Independent same-sign image data test insulated and Robin boundaries with
potential norms .00105963 and .000468463 at the finest spacing. Independent
Fourier Green series test periodic/antiperiodic grounded strips: three-mesh
potential norms decrease from .01674489 to .00447222 to .000909800, and from
.01326755 to .00207968 to .000248778. Axial translated pairs separately satisfy
full-revolution balance and total nodal voltage constraints. Prescribed normal
flux, zero-coefficient Robin limits and positive-reference conditions are tested.

The copied installation passes in plain LuaTeX, LuaLaTeX and LuaTeX-based
ConTeXt MkIV, with all 106 runtime hashes matching the checkout and actual load
provenance outside it. The generic fixture covers floating/tensor response,
source-offset matrices, interaction-energy quantities, nonlinear contour force
and mixed/periodic boundaries, in addition to typed nodes and views. The public
interface inventory covers 1070 declarations; 67 invalid electric PGF inputs
have explicit diagnostics. Total self-energy, native hybrid answer export,
additional volume/sheet loads and arbitrary 3D correction remain unavailable.

All 56 Lua suites and generic integration checks pass. Cold, warm and shared
frozen calculations agree in the three formats to the tested 1e-10 tolerance.
The final expanded response suite, including Fourier strip checks, also passes
separately. The new complete source-response and source-force examples compile;
their focused PGF chapter has been visually reviewed without overfull boxes.

The final manual has **354 pages** and executes **370 inline PGF blocks** plus
four source-file listings (**374 markers**), with **73 worked cases** and
**71 separate reference/loading fragments**. All **242 source fingerprints**
match the final sources. References and index are stable, with no overfull
boxes, missing glyphs or unresolved references. The source-correction chapter
starts on page 237. New source/result diagrams, capacitance chart and image
comparison on pages 246--247 and 250--251, the key reference and final index
were visually inspected.


## Native persistence and extended magnetic models

Native schema 1 is checked with analytical Coulomb/classical symmetry providers,
planar and axisymmetric electric FE/source correction, floating electrodes,
rotated dielectric tensors, Maxwell matrices, periodic traces and ordinary
static/harmonic magnetic fields. A voltage-driven winding is saved before its
first explicit solve, then restored with matching coefficients and terminal
power. Unknown schemas/capabilities/bases, damaged data, bad mesh indices,
callbacks and inconsistent fixed traces are rejected. Restoring answers is
checked with solver entry points disabled. These checks do not establish a
restart or general 3D-field capability.

The nonlinear equivalent-amplitude harmonic law agrees with linear complex
permeability in its linear limit and with the lossless fundamental law when the
loss angle vanishes. Its nonlinear fixed-point defect is checked independently
of the inner linear residual. Waveform history, higher harmonics and exact
cycle-averaged nonlinear stored energy are outside this approximation.
Square and CCA wire tests include DC limits, passivity, and terminal real/reactive
power identities. Independent FEMM 5 solves at 100 kHz give CCA vector-field
differences of 0.32--4.47% on unrelated meshes. Square-wire differences of
0.35--6.24% are recorded separately: native preprocessor and postprocessor use
different diffusion factors, whereas luafemm uses one consistent convention.
The native solver source revision was
`31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3`, without source patches.

Axisymmetric Kelvin mapping is checked against an exact mapped dipole and an
independently solved native coilgun example. The analytic field-error measure
decreases from 0.31673 to 0.06087 on refinement. The coilgun comparison uses
11513 Lua nodes and 43293 native nodes: three air probes differ by 3.23--4.68%
in B and 2.80--3.10% in flux. Scaled flux constraints, ordered ANS restoration
and opening the exported snapshot in an independent xfemm reader are checked.
Different interpolation and meshes limit the comparison; Windows GUI acceptance
has not been tested. Mapped coordinates are not physical exterior positions.

Quantitative magnetic maps have affine-patch, material-separated recovery,
vector-orientation and elliptical-phasor RMS checks. Invalid static RMS views,
reversed ranges and excessive vector sample requests fail. Generic controls,
styles, scoping and square/CCA native answer restoration execute in all three
formats, including the copied isolated installation. Nine shared usage scenes
exercise native storage, maps, refinement and mapped exterior drawing.

Planar field-driven refinement reduces a manufactured sine-field error from
0.59536 to 0.25211, with 25 then 216 nodes. A separate double-sine-series Poisson
reference gives harmonic-loss errors of 13.55% then 0.505%, and flux errors of
1.024% then 0.248%, after four refinements. The weighted conductor force remains
within 1.2e-6 relative of its exact J-cross-B value. Force errors need not decrease
monotonically. Narrow-interface material area, conforming connectivity and
transactional budget failures are checked separately.

For one localized Gaussian field, adaptive and uniform meshes give comparable
L2 field errors (0.35484 and 0.35509) with respectively 868 and 5184 nodes.
Observed CPU times were 0.321 and 0.441 s; retained Lua heaps after collection
were 2258 and 11255 KiB. These are local workload measurements, not portable
performance guarantees or peak process memory. The residual indicator is
heuristic and does not certify flux, loss or force accuracy.


### Prescribed-current three-dimensional fields

The Biot--Savart provider is checked against the circular-loop axial formula,
finite straight segment and square-loop formulae, superposed coaxial loops,
and Ampere's infinite-wire field inside/outside a uniform circular section.
A long finite convolved wire approaches that independent interior/exterior
limit. Off-axis helical fields agree with separately coded 12,000-panel
composite Simpson integration to the tested 2e-8 relative tolerance. Tests
also cover current/length scaling, rigid rotation, explicit feed closure,
finite fields on thick-wire centre lines, exterior divergence/curl differences,
polygonal circle refinement and convergence of discrete toroidal turns toward
the continuous winding limit.

True 3D integration around an infinite wire retains its axial coordinate and
circular radius to the tested 1e-6 relative tolerance. The camera is a subsequent
orthographic projection. Native prescribed-current answers retain source data
and the registered evaluation law; restoration is tested with field evaluation
disabled. Generic declarations, arbitrary planar paths, profile/plot options,
source freezing and native round trips execute in all three formats. Sixteen
invalid controls are explicitly rejected.

The thick curved-wire law is the documented divergence-preserving radial
convolution, not a solved metal conduction problem or a claim of uniform density
on every normal cross-section. General 3D material response, skin/proximity
physics, magnetic forces and imported external 3D fields are not validated by
these tests. Quadrature error estimates and visually smooth traces are not
certified field-error bounds.

### Reusable spatial scenes

Orthographic projection is checked against the existing Biot camera and an
orthonormal-basis identity. Scene ownership, plane/grid budgets, stable layer
order, edge-on axes and explicit failures are covered. Generic fixtures draw
multiple views with field evaluation disabled, retaining scoped styles and
TeX math labels. Twenty-eight invalid controls fail explicitly. Viewport clipping is
two-dimensional and does not change physical trace bounds.

The optional vector renderer is checked against analytical line/plane crossing
intervals, independent sloping-plane point classifications and independent
vertical-ray tests for two intersecting planes viewed from both sides. Tests
cover reversed cameras, collinear depth-order changes, crossing gaps, coplanar
and edge-on cases, non-occluding annotations, bounded failure and original-curve
mark positions. Generic fixtures disable the field evaluator while redrawing
both layered and vector views. Fragment decoration suppression is exercised in
all three formats. These tests concern centre-line visibility and opaque planar
surfaces; they do not validate solid tubes, finite stroke footprints, general
translucent composition or exact geometric predicates.


## Faceted conductor views

Tube tests independently check the oriented volume of a regular-prism wire,
its projected occlusion interval, 2037 independent oblique half-space rays,
camera reversal, and unchanged physical field
values when lead/return geometry is hidden. Generic scenes also test foreground
alpha ordering (100 pointwise painter-order probes) and rejection of a camera which places a translucent plane behind
an opaque object. Face budgets and invalid radius/size/side choices fail explicitly.
The manual includes executed tube and translucent-plane examples. These tests
establish the bounded display contract, not exact swept-solid geometry or
arbitrary transparent-volume rendering.


## Curved field marks and electric spatial capture

Arrow tests check local support discrepancy on a sampled circle, original
arclength, signed tangent orientation, size reduction, hidden-interval clearance,
multiple subpaths and work budgets. Native PGF bending is exercised on circles,
S-curves and both cameras of a helical toroid. Generic tests include short
finite-element contour segments, safe decimal TeX dimensions, scoped native tip
selection and ordinary picture transforms. The sampled local fit is a display
criterion, not a certified continuous geometric bound or new field accuracy.

Electric spatial checks compare all three E components with finite differences
of V, positive/negative source signs, interior/exterior charged-sphere laws and
infinite-cylinder translation invariance. Exact planar and axial finite-element
solutions check translation and the cylindrical-to-Cartesian basis. Integrated
paths retain field orientation. Captured arrays and physical samples are
unchanged under camera and visibility changes; equipotential roles remain
unarrowed. Executed examples compare a dipole section with spatial lines, a
charged sphere's placed section with an offset potential plane, and a revolved
finite-element cylinder. These checks do not establish arbitrary electric
isosurfaces, opaque source solids or a general three-dimensional FE solver.
