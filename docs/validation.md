<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Development-version validation

Version **0.8.0-dev**, verified on 2026-09-27. Local environment: TeX Live 2026,
LuaTeX/LuaHBTeX 1.24.0, PGF 3.1.12, PGFPlots 1.18.3 (`compat=1.18`),
ConTeXt MkIV 2026.09.10 10:02. Test commands and detailed logs are produced by
`make test`, `make examples`, `make manual` and `make dist` under `build/`.
A remote CI workflow is supplied; local verification does not imply a completed
GitHub Actions run.

## Ideal-circuit evidence

`tests/test_ideal.lua` compares a confined linear annulus against
`Phi = mu NI depth log(rout/rin)/(2 pi)` and checks its Ampere contour.
Across different full sections, flux agrees to a relative tolerance of 1e-9.
The same conservation checks cover a nonlinear unequal-section U, iron-to-gap
passage, a nonlinear asymmetric E (`Phi_center = Phi_left + Phi_right`) and a
tapered toroid with an alnico arc. The annulus accuracy assertions allow 1.2%
on flux and 2.5% on circulation on the selected polygonal mesh; these are
regression tolerances, not general error bounds.

The generic-format matrix additionally compiles an asymmetric ideal E,
its mean contour and integrated H profile with cold, warm and shared frozen
caches. Mesh, flux, profile and field calculations remain pure Lua. The new
shapes represent planar constant-depth geometry. They do not implement a 3D
conical-pole model, infinite permeability or irreversible demagnetisation.

## Automated evidence

- API: independent constructor/material/region inputs, finite sampling inputs,
  escaped numeric-key validation and no globals from Lua module imports.
- Solver: affine exact solution Az=x+2y and B=(2,-1), manufactured Poisson
  convergence, constitutive secant/tangent/extrapolation, zero excitation,
  invalid input rejection, nonlinear U electromagnet, current balance and
  reversal, contour continuity, constant potential and orientation.
- Permanent magnets: exact piecewise-affine magnet/air solution, coercive source
  sign, H components, disk reference, nonlinear alnico U, magnetization reversal
  and zero-coercivity limit.
- Boundaries: affine Dirichlet, Neumann, Robin and mixed-side solutions on grid
  and Delaunay meshes, two model-unit scales, nonzero data and nonlinear permanent
  magnets. A manufactured Robin problem has first-order L2 field convergence.
  A nonlinear symmetric half-U agrees with the full-device gap field within 1%.
  Tests reject incompatible Neumann loads (including tiny ones), corner conflicts,
  negative Robin coefficients and late declarations, and check gauge restoration.
  TikZ tests check picture-key ordering, parameterised and appended styles,
  repeated problem keys, scope persistence, inherited defaults and exact
  equivalence with command declarations. The format matrix also uses the key.
- Mesh: constraint recovery, area coverage, topology, positive triangles,
  material areas, oblique/concave boundaries, intersections, shared edges,
  holes and narrow gaps; 49 deterministic stress configurations.
- Compound domains: nested and intersecting contour areas on both meshers,
  nonzero/parity rules, reversed contours, multiple islands, current integrals,
  preservation of earlier inner materials, local refinement outside cavities,
  copied input data and cache invalidation for contour/fill-rule changes.
  A circular shell in a transverse field agrees with the analytic cavity field
  within 3% (about 0.3% on the tested mesh), using exact exterior Dirichlet data.
  TeX checks native fill-rule inheritance, transforms, curves and profile lifecycle;
  the three-format cache matrix includes a compound passive annulus.
- Quality: requested angle/area bounds recomputed from returned elements,
  conserved material areas and compound cavities, scale/translation changes,
  thin gaps, affine Robin boundary patch after exterior splitting, cache reuse
  and invalidation, invalid options and bounded failure on acute corners.
  The generic-format cache matrix exercises the quality controls through styles.
- Predicates: 4,530 orientation/incircle/diametral-disk signs compared with independent exact
  Python rational determinants, including uncertain cases that trigger the
  Lua expansion fallback.
- Materials: all 246 catalogue records instantiated, source and coercivity
  conventions, ten fill-factor curve transformations, 69 nonzero material
  currents, explicit region cancellation and independent metadata copies.
- Profiles: affine exact fields and potential, every component/abscissa,
  traversal direction, left-normal offsets, out-of-domain handling, curved
  path capture, caching and ownership of the associated model.
- TikZ: named/transformed coordinates, rotations, rounded corners, nodes on
  physical paths, inherited material settings, PGF current expressions,
  profile plotting, legend/cycle state and deferred solving.
- Components: core/armature/winding areas, physical anchors, invalid parameter
  rejection and preservation of signed ampere-turns under transformations.
  TeX tests check text/font/padding independence, rotated/scaled placement,
  positioning-library use, multiple sizes and anchor queries after model changes.
  The generic-format matrix also solves a component electromagnet, samples both
  gaps and reuses its solution with changed node text and display scale.
- Persistent cache: grid and Delaunay round trips preserve nodal potentials,
  B/H components, profiles and contours exactly. Tests cover physical changes
  that reuse the mesh, geometry/mesh changes that invalidate it, disabled and
  refresh modes, missing/stale frozen results, mesh-only records, corrupt and
  incompatible files, malformed geometry, failed nonlinear solves, interrupted
  replacement, unwritable destinations and unsupported closure callbacks.
- Formats: identical nonlinear alnico/iron geometry and a curved profile in
  plain LuaTeX, LuaLaTeX and ConTeXt MkIV. A mixed-boundary affine patch exercises
  SI expressions, independent declaration defaults and all three boundary types.
  Bx, By, node-count, component-field and boundary-field snapshots
  agree within 1e-10. The tests verify loader idempotence and catcode restoration.
  Each format runs cold, warm with changed display scale/styles/sample counts,
  and frozen from the same cache originally written by plain LuaTeX. Warm and
  frozen runs emit no Newton iteration and frozen reads leave the file intact.

Plain and LaTeX documents are compiled with shell escape disabled. The stock
MkIV launcher enables shell escape itself; luafemm does not invoke external
commands or depend on that setting.

The example and manual builds reject errors, missing characters and overfull
boxes. The manual executes its worked examples and builds its key index.
The distribution check verifies generated-data reproducibility and replays the
format tests against an extracted TDS installation outside the runtime tree.

## Manual examples

Every standalone document under `examples/` is included in the manual as a
complete source/output pair: 14 documents, including the three format wrappers.
Their shared `scene.tex` has its own listing. The coverage check rejects omitted
examples and stale names. PGF's native `codeexample` machinery reads the actual
files, and the manual includes every page of the resulting PDFs. pdfcrop with
Ghostscript removes empty paper margins from the manual's copies, preserving
vector graphics; the original example PDFs remain unchanged. Short complete
snippets also execute in place; partial declarations remain code-only fragments.

The opening tutorial is compiled from the same two standalone sources supplied
under `examples/tutorial-coil.tex` and `examples/tutorial-alnico.tex`.
`make manual` builds all example PDFs once before the manual's multiple index
passes; the resulting field and profile figures are embedded directly. This
also runs the examples through the same log checks as `make examples`.
The tutorial's own small TikZ construction is executed inside the manual.

Both tutorial documents now enable automatic persistent caching. Their cache
files follow the TeX output directory; rebuilding the manual reuses the solved
models when only the presentation or manual text changes.

The draft cases use 5 mm global spacing, with 1 mm local spacing in the coil's
3 mm gaps and 0.5 mm in the permanent magnet's 1 mm gaps. They contain 4,852
and 4,970 triangles respectively. These counts are observations for the current
sources, not mesh-size guarantees. The figures show the draft solutions openly,
and the tutorial separates mesh refinement from contour-level counts and profile
sample counts. The two finer settings in its table are also compiled locally;
this checks that the suggested examples run, not that either mesh has a certified
accuracy. Geometry, labels, arrow orientation, profiles and manual page layout
are inspected in rendered PDFs.

## Previous sensitivity measurements

These retained prototype measurements document sensitivity; they are not newly
measured release guarantees. Cocircular diagonal choices and Newton counts can
change when the mesher changes. A direct numerical comparison against FEMM or
xfemm remains to be performed.

For the alnico U, By is sampled at (-32.5,30.5) and (32.5,30.5) mm. All cases
converged in six Newton iterations.

| Global / gap spacing (mm) | Domain (mm²) | Nodes | Triangles | Left By (T) | Right By (T) | Relative residual | Minimum angle |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 / .5 | 200 × 185 | 11376 | 22364 | .782876 | −.782907 | 2.66e-12 | 1.74° |
| 1.5 / .25 | 200 × 185 | 20409 | 40300 | .777974 | −.778597 | 2.11e-11 | .874° |
| 2 / .5 | 280 × 265 | 22293 | 44038 | .775140 | −.775140 | 1.02e-11 | 2.00° |

Refinement changes the gap field by roughly 0.6%; enlarging the domain changes
it by roughly 1%. Refinement lowers the smallest angle: local density is not
quality refinement with an angle target. These historical runs did not enable
the new bounded quality-refinement option.

For the nonlinear electromagnet, the corresponding Delaunay prototype measured:

| Spacing (mm) | Domain (mm²) | Nodes | Triangles | Left By (T) | Right By (T) | Newton |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 3 | 220 × 210 | 6363 | 12436 | 1.115694 | −.958717 | 11 |
| 2 | 220 × 210 | 13849 | 27264 | 1.103039 | −.955858 | 10 |
| 1.5 | 220 × 210 | 24314 | 48052 | 1.097133 | −.950098 | 12 |
| 2 | 320 × 310 | 29307 | 57982 | 1.104500 | −.939216 | 11 |

Neither a tiny algebraic residual nor these relative changes bound the total
physical error. Domain truncation, geometry flattening, material accuracy and
finite-element resolution must be checked for each application.
