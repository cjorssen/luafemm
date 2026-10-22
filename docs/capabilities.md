<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Capability matrix — 0.10.0-dev

This matrix describes the supported **planar, static** release, not all physics
represented by the FEMM material catalogue. “Supported” means implemented and
covered by the named checks within the stated scope. “Partial” means an explicit
restriction or approximation. “Deferred” means unavailable, not an alternative
interpretation of an existing key. No row is a general accuracy guarantee.

The package uses pure Lua at runtime and generic TikZ/PGF/PGFPlots interfaces.
Independent xfemm executables are optional development oracles only.

## Geometry, physics and output

| Capability | Status | Scope and limits | Evidence |
| --- | --- | --- | --- |
| plain LuaTeX, LuaLaTeX, ConTeXt MkIV | Supported | Same PGF styles and numerical results; MkIV explicitly uses LuaTeX. | `scripts/build.py formats`; extracted TDS installation |
| Arbitrary TikZ paths and transformations | Supported | Curves are flattened to straight constrained segments with a declared tolerance; no exact curved finite elements. | `tikz-interface.tex`, `test_api.lua`, TikZ integration tests |
| Nested/compound regions, holes and exclusions | Supported | Fill rules and paint order; computational exclusions differ from material holes. | `test_domains.lua`, `test_interchange.lua` |
| Parametric U/E/toroid/rotor/stator shapes | Supported | Planar constant-depth models, slots, conductor layers, signed winding weights; no 3D conical model. | Component, shape, machine tests; Kamil tutorial |
| Constrained triangulation and quality refinement | Partial | Pure Lua; bounded refinement can fail explicitly on unsuitable geometry. Periodic synchronisation can reduce angles after refinement. No adaptive field-error estimator. | `test_mesh*.lua`, `test_quality.lua`, predicate oracle |
| Planar isotropic nonlinear magnetostatics | Supported | Single-valued B–H curves; optional FEMM interpolation. No magnetic history/hysteresis integrator. | `test_solver.lua`, `test_interchange.lua`, xfemm comparisons |
| Catalogue | Partial | 246 magnetic records; stored frequency/loss/wire metadata do not activate absent physics. | `test_materials.lua`, deterministic regeneration |
| Linear anisotropy and sheet mixtures | Supported | Rotated positive tensors; linear parallel/perpendicular sheet laws. | `test_static.lua`, diagonal-tensor xfemm case |
| Nonlinear anisotropy, detailed strand/wire physics | Deferred | Isotropic nonlinear curves remain supported. | Explicit rejection tests for nonlinear anisotropy |
| Permanent magnets | Supported | Uniform, radial, tangential or angular coercive direction; three-point triangle quadrature. No irreversible demagnetisation history. | Magnet tests, `test_static.lua`, angular xfemm case |
| Dirichlet/Neumann/Robin boundaries | Supported | Rectangle sides or named paths; internal Dirichlet constraints. Robin coefficient must be nonnegative. | `test_boundary.lua`, `test_static.lua` |
| Periodic/antiperiodic boundaries | Supported | Rigid paired traces, synchronised nodes, exact signed DOF elimination; repeated-field rendering. | `test_static.lua`, `static_convergence.lua`, format/cache matrix |
| Open exterior | Partial | Circular multipole Robin approximation; adjacent air and zero net longitudinal current. No arbitrary exact infinite domain or shell transform. | Analytic dipole and off-centre radius/mesh study |
| Named DC circuits | Supported | Current/voltage drive; series/parallel resistive winding branches; prescribed uniform leg currents. No induced-current conductor unknowns. | `test_static.lua`, generic formats, xfemm circuit case |
| Circuit observables | Partial | Resistance, power and gauge-independent closed-winding linkage; secant inductance where documented. No incremental or general multiport inductance matrix. | Ohmic and energy/linkage identities |
| Ideal confined circuits | Supported | Flux conservation, junction balance, mean contours and Ampere profiles; distinct from finite-permeability open-field models. | `test_ideal.lua` |
| Field lines and quantitative profiles | Supported | Oriented Az contours, B/H components, polar scans, spatial Fourier analysis, PGFPlots curves. | Profile, integral, machine and generic-format tests |
| Regional integrals | Supported | Area, volume, current, B/H averages, nonlinear energy/coenergy under the stated magnet reference. | `test_post.lua`, xfemm postprocessor |
| Force and torque | Partial | Planar Maxwell contour/weighted stress and conductor Lorentz methods; assumptions on selected bodies/air and mesh sensitivity apply. No harmonic components. | `test_post.lua`, `post-interop.py`, three-mesh armature study |
| Scalar colour maps and vector arrows | Deferred | Force arrows exist; they are not a general B/H vector renderer. | No implementation claim |
| Persistent native caches | Supported | Mesh/solution reuse and cold/warm/frozen modes; passive format 7. Recompute older snapshots. | `test_cache.lua`, three-format frozen comparisons |
| Axisymmetry | Deferred | Requires its own operator, axis regularity, measures and exchange conventions. | Roadmap |
| Time-harmonic magnetics and eddy currents | Deferred | No complex field solve, skin/proximity effects or harmonic loss prediction. Spatial winding harmonics are not time-harmonic physics. | Roadmap |
| Frozen/incremental operating-point solves | Deferred | No prior-solution linearisation or incremental inductance. | Roadmap |
| Sliding air gaps and motion | Deferred | Static periodic sectors exist; no sliding-gap coupling without remeshing. | Roadmap |
| Electrostatic/current-flow/thermal solvers | Deferred | Separate planned modules, including their own material catalogues and formats. | Roadmap |

## FEMM interchange

Native TikZ curves remain available. Export approximates them by segments;
imported circular arcs are flattened with the declared tolerance. An unsupported
physical feature is rejected instead of being silently replaced. See
[interoperability](interoperability.md) for the precise contract.

| Feature | FEM import | FEM export | ANS output |
| --- | --- | --- | --- |
| Planar DC geometry/materials | Supported subset | Supported subset; polygonal curves | Solved mesh and nodal potential |
| Linear anisotropy | Diagonal tensors and supported sheets | Diagonal model-frame tensors; rotated cross terms rejected | Same representability restriction |
| Spatial magnetisation | Restricted passive angle grammar | Canonical safe angular expression | Elementwise source snapshot preserves discrete H as well as A/B |
| Periodic boundaries | Straight paired chains; periodic arcs deferred | Paired equal-length straight primitives | Signed nodal pair trailer |
| Circuits | Supported series records; parallel import deferred | Resolved operating currents in series winding records | Operating-point currents; spatial-source snapshots can lose circuit association |
| Native ideal model | Unsupported | Unsupported | Unsupported |
| Native cached solution | Separate native format | Separate native format | Export is not a lossless replacement for the native cache |
| Loading an existing FEMM ANS solution | Deferred | Not applicable | Writing only |

Export does not preserve a native voltage controller or parallel network as an
editable FEMM circuit network. The original TikZ source and native cache retain
those semantics. Windows GUI acceptance remains **unverified**; native xfemm
reader and independent solver comparisons validate only the documented fixtures.

## Evidence and release boundary

Run `make test`, `make static-convergence` and `make dist`. The last rebuilds
all 30 source/output examples and the indexed manual, regenerates the catalogue,
and tests an extracted TDS installation in all three formats. Optional independent
checks are `scripts/interop.py`, `scripts/static-interop.py` and
`scripts/post-interop.py`; they require the development xfemm build.

[Validation](validation.md) records measured errors and sensitivity results.
Agreement with xfemm on a shared mesh tests representation/postprocessing;
it does not measure discretisation error. Mesh refinement, domain enlargement
and geometric tolerance are separate controls. Cancellation between their errors
can make one apparently accurate coarse result misleading.

The [capability roadmap](femm-capability-roadmap.md) preserves the original audit
and future milestones. The [Tangle audit](tangle-audit.md) remains a deferred
engineering recommendation; this release does not include a Tangle port.
