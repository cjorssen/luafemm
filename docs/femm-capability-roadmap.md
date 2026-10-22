<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# FEMM capability audit and development roadmap

Date: 2026-10-22. Baseline: luafemm 0.9.0-dev,
`93620f0d04569877491a0c1079ea27bd32ab79f4`.
This is a capability assessment and proposed sequence, not an implementation
promise, release schedule or reference for currently available keys.
For the current release status, use the [capability matrix](capabilities.md);
the tables below retain the historical audit and are not a current backlog.
The [Tangle audit](tangle-audit.md) is retained separately for later work.

Implementation update (0.10.0-dev): the planar postprocessor now provides named
material/path selections; area, volume, current and B/H averages; nonlinear
energy/coenergy; weighted and contour Maxwell forces/torques; conductor Lorentz
forces/torques; force glyphs, formatted results and PGFPlots snapshot series.
See `doc/sections/postprocessing.tex` and `examples/forces.tex`. The baseline
assessment below is retained historically.

Static implementation update (0.10.0-dev): custom native domains/path frontiers;
strong periodic/antiperiodic constraints with synchronised traces; circular
multipole open boundaries; rotated linear tensors and linear sheet mixtures;
passive angular magnetisation; DC named series/parallel windings with
current/voltage drive; resistance, power, closed-winding linkage and secant
inductance are implemented. See `doc/sections/static.tex`. These cover the first
static increment, not every extension in the historical table: nonlinear
anisotropy, arbitrary imported expressions, periodic arc import, full FEMM
parallel-circuit import, point sources, incremental inductance, sliding gaps,
axisymmetry and harmonic physics remain open. Loss models, scalar colour maps
and field-vector plots are still future work. Tangle work remains deferred.

## Scope and evidence

The target is equivalent supported physical models, numerical observables and
problem/solution interchange expressed through TikZ/PGF/PGFPlots. Mouse tools,
GUI windows, printing dialogs and legacy Lua 4 commands do not need literal
copies. TikZ paths, nodes, pics, transformations and styles remain the modelling
language; PGFPlots handles quantitative plots. Runtime stays pure Lua and generic
TeX, with plain LuaTeX, LuaLaTeX and ConTeXt MkIV supported together.

Local evidence: FEMM revision `7d9e8ed90772ddb939a3256f976e231bcb5cca97`,
particularly `manual/magnlua.tex`, `manual/eleclua.tex`, `manual/heatlua.tex`,
`manual/cflowlua.tex`, `fkn/`, `belasolv/`, `hsolv/` and `csolv/`.
Our baseline is checked against the runtime modules, [interchange limits](interoperability.md),
[meshing contract](meshing.md), README and existing validation records.
The [official documentation index](https://www.femm.info/wiki/Documentation/)
also links the FEMM 5 manual dated 2026-10-01 and its format specifications.
That manual confirms planar/axisymmetric magnetic, electrostatic, thermal and
current-flow scope, including transient heat flow. A complete FEMM 5 parity claim
will need a separately pinned feature/format inventory; the local source revision
and the online manual must not be assumed identical.

Physical scope is bounded: full 3D, electromagnetic wave propagation, arbitrary
magnetic hysteresis history, general transient magnetic field-circuit-motion
simulation and fluid flow are separate projects, not automatic requirements of
the established FEMM magnetostatic/time-harmonic formulations. Frozen and
incremental permeability supply useful linearizations, not a complete transient
magnetic integrator. FEMM's harmonic loss models do not reconstruct arbitrary
hysteresis loops or irreversible magnet histories.

## Current magnetic baseline

- Planar DC, P1 triangular finite elements; isotropic nonlinear B-H laws and
  prescribed permanent-magnet direction; optional FEMM-like B-H interpolation.
- Prescribed current densities, component ampere-turns and winding weights;
  imported series circuits converted to the supported excitation model.
- Native TikZ contours, compound fills, nested regions, explicit topology for
  imports, holes, constrained meshing and bounded quality refinement.
- Dirichlet, Neumann and nonnegative Robin data, native rectangular side keys,
  supported imported boundary tags, compatibility checks and potential gauges.
- Oriented magnetic field contours, mesh/geometry rendering, B/H profiles,
  circulation and section flux, angular scans and spatial harmonic analysis.
- U/E/toroidal devices, rotor/stator slot geometry and an additional ideal-flux
  model; the latter is a luafemm extension, not FEMM interoperability parity.
- Computation caches, a planar DC subset of FEM import/export, ANS output and
  independent xfemm checks. ANS input is not implemented.

The material catalogue contains 246 records, but stored properties do not imply
all associated physics is active. Conductivity, loss and winding metadata are
not a harmonic solver. Some isotropic static mixtures are supported; anisotropic
lamination laws and the broader import/material set remain incomplete.

## Remaining magnetic work

Effort labels below are relative engineering assessments, not calendar estimates.
Every row includes Lua implementation, generic keys, documentation and independent
tests; geometry or a user-visible key alone is not completion.

| Area | Remaining work | Effort/dependencies |
| --- | --- | --- |
| Native boundaries and sources | Attach named conditions to arbitrary native curves and internal graph constraints; point potential/current properties; shared interfaces; paired periodic/antiperiodic boundaries. Preserve IDs through flattening and splitting. | Medium/high; topology, mesher and DOF constraints together. |
| Open domains | Validated asymptotic/mixed boundaries, shell-based open-boundary constructions and axisymmetric exterior transformations. Demonstrate domain-size convergence. | Medium/high; axisymmetric transformations depend on that formulation. |
| Materials | Linear anisotropy, supported lamination directions/fill laws, complete wire/strand models where applicable, spatial magnetization definitions and supported nonlinear combinations. Keep unsupported cases explicit. | High; real tensors first, complex constitutive laws later. |
| Axisymmetry | Correct magnetic unknown and radial factors, regularity at r=0, r/z components, flux interpretation, 2*pi*r volume measures and FEM/ANS conventions. | High; not a coordinate-label change. |
| Circuits | Named electrical circuits with signed turns, series/parallel semantics, conductor unknowns where required, circuit current/voltage/flux linkage and resistance; consistent native/interchange representation. | High; DC first, harmonic coupling later. |
| Harmonic magnetics | Complex A, B, H and sources, eddy-current coupling, skin/proximity effects, conductive-region current constraints, harmonic material approximations and losses. Audit phasor sign and peak/RMS conventions. | Very high; new complex algebra and suitable linear solvers. |
| Operating-point modes | Read previous solutions; frozen and incremental permeability about a validated nonlinear state, including incremental force/flux-linkage observables. | High; material tangents, mesh identity and solution provenance. |
| Rotating machines | Periodic sectors, synchronized boundaries, sliding air-gap formulation, angle sequences and torque/flux-linkage vs. angle without remeshing noise. | Very high; boundary, circuit and postprocessing foundations. |
| Integral postprocessing | Named region/group integrals: area, volume, current, field averages, magnetic energy/coenergy, flux linkage, inductances and losses. Respect nonlinear and permanent-magnet conventions. | Medium/high; observables depend on available physics. |
| Forces and torques | Lorentz contributions, Maxwell contour stress, weighted stress tensor and coenergy checks; appropriate air contours and bodies; DC and harmonic components. | High; weighted stress needs an auxiliary problem, not a single B sample. |
| Visualization | Scalar colour maps with colour bars, vector plots, named selections and optional documented field recovery/smoothing; complex magnitude/phase and snapshots. Keep raw element values available. | Medium; sampling and units must remain consistent. |
| Interchange | Extend FEM/ANS coverage with each physical feature; load native FEMM solutions, preserve identifiers/encoding policy, add controlled editing of imported models and Windows GUI acceptance checks. | High across milestones; preserve the polygonal TikZ export contract. |
| Meshing and scale | Deferred Tangle improvements, narrow-feature tests, better performance and memory use, convergence studies and optional adaptive refinement. | Medium/high; performance is not field-error certification. |

Inductance needs particular care: in a nonlinear biased problem, secant and
incremental inductances differ. Energy/coenergy also require a declared material
convention. Force evaluation is especially sensitive to discretization and must
be checked by independent methods and mesh refinement.

An AC extension cannot reuse conjugate gradients merely by replacing real
numbers with complex tables. Matrix symmetry, definiteness, circuit constraints
and preconditioning must be revisited. Similarly, a spatial winding Fourier
analysis already available in luafemm is not a time-harmonic field solve.

## Other FEMM modules

These modules share geometric and finite-element infrastructure, but have their
own unknowns, constitutive laws, boundary conditions, conservation tests and file
formats. Source locations below are relative to the pinned FEMM checkout.

### Electrostatics

Evidence: `belasolv/prob1big.cpp`, `belasolv/mesh.h`, and
`manual/eleclua.tex` (materials, conductors, point values and integrals).

Solve for electric potential V with E = -grad(V), D = epsilon E and
div(D) = rho. Start with linear isotropic dielectrics, then the supported
anisotropic permittivity, planar/axisymmetric geometry, volume/surface/point
charges and prescribed potentials. Add floating equipotential conductors with
specified total charge, mixed/periodic boundaries and open-domain treatment.

Postprocessing needs E/D components and magnitudes, potential contours, field
lines, charges, energy, conductor properties, capacitance coefficients and
electrostatic forces/torques. The native files are `.fee` and `.res`.
Suggested validation: parallel plates, coaxial geometry, dielectric interfaces,
a floating conductor and reciprocal capacitance coefficients. A planar linear
prototype is relatively small; complete boundary/conductor/force parity is a
substantial separate stage.

### Current flow

Evidence: `csolv/MESH.H`, `csolv/PROB1BIG.CPP`, and
`manual/cflowlua.tex` (frequency, anisotropic conductivity/permittivity,
conductor currents and point values).

Start with DC conduction: E = -grad(V), J = sigma E and charge conservation.
Add terminals specified by voltage or total current, insulating/flux/mixed
boundaries, anisotropic conductivity, floating/equipotential conductors and
planar/axisymmetric geometry. The harmonic extension uses complex admittivity
sigma + j*omega*epsilon for the chosen phasor convention, with the dielectric
loss model accounted for consistently. Report conduction and displacement
currents separately, impedance, dissipated power and terminal properties.

This electric-potential formulation does not replace magnetic induction and
eddy-current physics in the magnetic solver. The native files are `.fec` and
`.anc`. Validate a resistor, interface current balance, parallel conductors and
an analytical harmonic RC configuration. DC is relatively accessible;
harmonic parity benefits from the common complex algebra work but still needs
its own formulation and tests.

### Heat flow

Evidence: `hsolv/MESH.H`, `hsolv/prob1big.cpp:179` (time-step contribution),
`manual/heatlua.tex:150` (previous solution/time step), `:320` (materials) and
`:334` (boundary properties).

Solve for temperature T with heat flux q = -k grad(T). Steady conduction needs
volumetric sources, anisotropic conductivity, temperature-dependent k(T), fixed
temperature, prescribed flux, convection and nonlinear radiation. Isothermal
conductors may have prescribed temperature or total heat flow. Preserve
Kelvin-based radiation evaluation even if display units are Celsius.

FEMM also supports transient heat flow: implement volumetric heat capacity,
initial/previous temperature fields, time-step assembly and consistent restart
semantics. Check time-discretization and mass-lumping conventions against the
reference solver rather than assuming a stationary solve with a new source is
equivalent. Planar and axisymmetric cases are both needed.

Postprocessing needs temperature maps/isotherms, heat-flux vectors and streamlines,
profiles, heat rates, volume averages and thermal balances. Native files are
`.feh` and `.anh`. Validate a slab, radial conduction, a convection boundary,
radiative equilibrium and an analytical transient diffusion mode. Linear steady
conduction is a modest first stage; nonlinear boundaries and transient parity
make the complete module a larger undertaking.

### Coupling beyond separate modules

Magnetic/current-flow losses can drive heat sources; temperature can feed back
into conductivity and magnetic properties. This is an additional orchestration
layer, not automatic functionality obtained by implementing the four solvers.
It requires conservative source transfer between meshes, dependency-aware cache
invalidation, coupled convergence checks and explicit time/frequency conventions.
Fluid convection simulation and mechanical deformation are separate physics;
prescribing a heat-transfer coefficient or rotor angle does not implement them.

## Shared architecture and rendering

Extract reusable geometry, topology, meshing, sparse assembly, boundary tagging,
unit handling and passive caches progressively. Keep physics-specific residuals,
material derivatives, source terms and observables separate. A common scalar
diffusion assembler is useful for electrostatics, conduction and heat; it must
not obscure the distinct axisymmetric magnetic operator or harmonic circuit DOFs.

Namespaced key families should select a problem's physics, geometry and analysis
type with `.is choice`; material, region, boundary, terminal and output families
should compose through normal PGF `.style` definitions. No public key spelling
is fixed here. Preserve the existing magnetic interface and reject incompatible
material/physics combinations. Store physics, frequency, time, operating point,
geometry, material law and algorithm revisions in computation fingerprints.

One important rendering distinction: magnetic field lines in our planar model
are contours of Az. Electrostatic and thermal potential contours are instead
equipotentials/isotherms; E and heat flux follow negative gradients. A generic
streamline tracer or suitable flux-potential construction is needed for these
vector fields. With anisotropic coefficients, flux and gradient need not even
be parallel. Do not relabel scalar contour plots as field lines.

A physical catalogue must also be selected by physics. The magnetic catalogue's
nonmagnetic entries are not a thermal or dielectric database. Audit the other
FEMM catalogues, units, temperature/frequency data and redistribution notices
before importing them. Original code remains LPPL; upstream data or translated
code retain applicable notices as explained in `LICENSES.md`.

## Recommended milestones and completion criteria

1. **Consolidate planar DC utility:** named integrals/selections, energy and
   coenergy, circuit observables, validated force/torque calculations, scalar
   and vector rendering. Start a machine/actuator quantitative tutorial.
2. **Broaden static models:** arbitrary native boundary tagging, periodicity,
   anisotropy, spatial magnetization, open boundaries and fuller static circuits.
   Use the deferred Tangle work where measured failures or cost justify it.
3. **Add axisymmetric magnetics:** validate axis regularity and all quantities,
   file exchange and plots with analytic/reference cases.
4. **Add harmonic magnetics:** complex algebra, conducting regions/circuits,
   materials and losses, with explicit phasor conventions and skin-depth tests.
5. **Complete advanced machine workflows:** previous-field modes, sliding gaps,
   angle sweeps and incremental observables; extend exchange throughout.
6. **Add a scalar-physics branch:** electrostatics first to test the shared
   architecture, followed by DC current flow and stationary heat. This branch
   can begin after the static infrastructure is stable without waiting for all
   harmonic machine features.
7. **Complete the additional modules:** harmonic current flow, nonlinear and
   transient thermal models, axisymmetric coverage, full postprocessing and
   problem/solution exchange; consider coupling only after independent validation.

For each claimed feature require: a documented weak formulation and units;
analytic/manufactured and independent FEMM/xfemm comparisons as available;
mesh/domain/time convergence; conservation and sign checks; invalid-input tests;
cache and file round trips; identical supported behaviour in all three TeX
formats; and complete English key documentation with executed source/result
examples. Native tools remain optional development references. No numerical
solver tests were rerun for this documentation-only assessment.

The existing comparative checks are useful evidence for their specific DC
cases. They are not a parity certificate for all FEMM physics. A future parity
matrix should track each material, boundary, observable and interchange feature
as supported, partial, rejected or unverified, against pinned references.
