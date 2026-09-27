<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Changelog

## 0.9.0-dev — 2026-09-27

- Add rotor and circular/square stator nodes with compound iron domains, shafts,
  casings, detailed slot placement and shape controls, and conductor layers.
- Add balanced coil records, per-conductor materials and excitation, and uniform
  or sinusoidally weighted distributed windings with conserved ampere-turns.
- Add local air-gap refinement, exact circular scans, polar field components,
  mechanical/electrical angle plots, and spatial Fourier spectra and statistics.
- Add Kamil's four-stage machine tutorial and a complete node-construction
  example, with executed source/output pairs and an extensive key reference.
- Validate machine geometry and numerical fields, generic-format cache reuse,
  and optional independent xfemm interchange comparisons.

## 0.8.0-dev — 2026-09-27

- Add confined planar magnetic domains with one floating potential per window.
  Enforce zero wall flux, exact discrete flux conservation and junction balance;
  retain nonlinear materials and permanent-magnet sources.
- Add ideal-domain membership and queued window-excitation keys, exact section
  flux output, geometric/half-flux mean contours and cumulative H circulation
  in PGFPlots. Preserve native TikZ paths and generic-format support.
- Add E core/electromagnet and toroidal/tapered-toroidal node shapes, unequal
  U/E section widths, independent E gaps and optional toroidal magnet arcs.
- Expand the English manual with all new keys, conventions, numerical limits
  and four complete examples paired with their compiled output.
- Upgrade private caches to format 6; test ideal cold/warm/frozen snapshots
  across plain LuaTeX, LuaLaTeX and ConTeXt MkIV.
- Reject ideal-mode FEMM interchange explicitly: floating-window constraints
  are outside the existing FEM/ANS subset. Ordinary interchange is unchanged.

## 0.7.0-dev — 2026-09-27

- Remove the prototype command compatibility layer and its manual appendix.
  Use `femm/problem`, native `femm/region` paths and `pic {femm field}`.
  Remove old Lua-to-TeX forwarding aliases and the model-region `points` alias;
  use the TeX bridge directly and `region.contours[1]` respectively.
  Rewrite the minimal plain-LuaTeX example using the current key interface.

- Add a passive magnetic FEMM 4.0 reader/writer and planar DC import with units,
  circular arcs, labeled faces, series-current excitation and true No Mesh holes.
- Preserve native TikZ curves; export their computational polygonal interfaces.
- Generalize imported domains and boundary assembly, including disconnected
  components with independent Neumann compatibility checks and gauges.
- Write static non-incremental ANS files from the actual converged mesh and
  potentials; reject stale solutions and incompatible nonlinear material laws.
- Add natural cubic FEMM B-H interpolation, bounded smoothing and endpoint-slope
  extrapolation as an explicit policy; retain linear interpolation by default.
- Add generic import/export keys, named delayed exports and geometry styles;
  test interchange under plain LuaTeX, LuaLaTeX and ConTeXt MkIV.
- Bump the private cache to format 5 and validate carved-domain coverage.
- Add independent xfemm disk-reader and solver checks, plus three complete
  interchange tutorials with their rendered output. Windows GUI validation
  remains outstanding; no native library is required by the package.

## 0.6.0-dev — 2026-09-27

Development series following the first public release candidate.

- Optional `minimum angle`, `maximum area` and `max Steiner points` problem
  keys for bounded, pure-Lua Delaunay quality refinement. Encroached constraints
  are bisected before circumcentre insertion; local flips and a final audit
  preserve interfaces, exterior tags, topology and requested quality bounds.
- Exact diametral-disk predicate, deterministic triangle priority, insertion
  diagnostics, and explicit failure at budgets or geometric precision limits.
  This is not a full Triangle port or an acute-corner protection algorithm.
- Before/after U-electromagnet mesh example with full source/output in the
  manual, quality reference chapter, independent quality/Robin/cache tests,
  and 4,530 rational predicate comparisons. Three-format cold/warm/frozen tests
  include quality refinement; cache format 4 records its controls and statistics.

- Complete source and actual output for every distributed example in the manual,
  using PGF's native code highlighter and including the three format wrappers.
  Manual builds compile all examples and reject missing source/output pairs;
  short standalone lifecycle and cache examples execute beside their code.

- Compound physical paths use native TikZ `even odd rule` and `nonzero rule`.
  Nested holes preserve existing media or background air; all subcontours are
  constrained on both meshers. Local Delaunay refinement follows the filled area.
- Shared winding classification, a compound-region Lua API, complete contour
  cache descriptors (revision 3), and preserved single-curve profile semantics.
- Explicit lifecycle documentation covering styles, picture hooks, physical
  declarations, mesh freezing, deferred solving, profile registration and caches.
- Nested alnico/air/iron example and tests for material areas, current integrals,
  fill rules, local refinement, cache changes and analytic circular-shell shielding.

- Native `femm/boundary` TikZ key, usable in picture options, reusable styles
  and setup scopes. Picture declarations are deferred until model creation,
  independently of their order relative to `femm/problem`; body declarations
  share the command implementation. Scoped defaults do not leak across models.
- Repeated `femm/problem` keys in composed styles initialise one model with
  the final options and apply its boundary declarations once.

- Key-based exterior Dirichlet, Neumann and nonnegative Robin conditions,
  independently selectable on all four sides, with affine SI data and PGF
  expressions. Default zero potential remains unchanged.
- Consistent edge assembly, conflicting-corner diagnostics, all-Neumann current
  compatibility and a deterministic potential reference. Boundary changes
  reuse geometry but invalidate the solution (cache revision 2).
- Boundary reference chapter, symmetric half-U example and analytic tests on
  both meshers, including nonlinear magnets, convergence and three-format caches.

- Genuine PGF component nodes: U core and complete U electromagnet, with fixed
  dimensions independent of text, named physical anchors, configurable parts
  and per-part graphic styles. Components use existing simple region paths.
- Winding ampere-turns normalized by captured section area, including local
  scaling; component-local magnetization directions follow transformations.
- Pure Lua component geometry shared by paths and anchors, single registration
  during drawing, and clear errors for invalid geometry or delayed placement.
- Complete node-based field/profile example, generic-format cache coverage and
  regression tests for anchors, text independence, transforms and positioning.

- Opt-in persistent cache with `femm/cache={file=...,mode=auto}` and PGF
  choice modes `auto`, `refresh`, `off` and `frozen`.
- Separate exact descriptors for meshing and physical/solver inputs: changing
  excitations or material laws reuses the mesh; presentation and measurement
  changes reuse the converged solution. Current picture frames are retained.
- Passive, versioned `.lfc` records with round-trip numeric precision, integrity
  checks, structural validation and replacement through a neighbouring temporary
  file. Missing or stale frozen results are errors. Cache files are not FEMM
  `.fem` or `.ans` files; exchange with FEMM remains future work.
- Relative TikZ cache paths respect TeX's output directory. `\femmcachestatus`
  and Lua `model.cache_info` report reuse; numerical statistics retain the
  original calculation's iteration counts and timings.
- Idempotent mesh/solve requests, a shared mesh-preparation path and English
  cache reference documentation. Tutorial sources enable automatic reuse.
- Regression coverage for both meshers, physical/geometric invalidation,
  damaged records, interrupted writes, callback rejection, and cold/warm/frozen
  compilations across plain LuaTeX, LuaLaTeX and ConTeXt MkIV.

## 0.5.0-rc.1 — 2026-09-25

First release candidate of the generic TikZ and PGFPlots libraries.

- Opening tutorial with two complete U-magnet examples, field lines and
  contour profiles, following the progressive presentation of the PGF tutorials.
  Coarse drafting meshes, local gap refinement and final convergence checks
  are explained separately from graphical detail and profile sampling.

- Original work licensed under LPPL 1.3c, with Christophe Jorssen as author and
  Current Maintainer; upstream FEMM data keep their separate terms.
- Author contact, GitHub development/issue links and an explicit Codex/GPT-6-Astra
  "vibe coding" notice with guidance on independent numerical validation.

- Installable TEXMF layout, with plain LuaTeX, LuaLaTeX and ConTeXt MkIV loaders.
- Explicit `\usepgfplotslibrary{femm}` entry point; PGFPlots stays optional.
- Lua numerical modules separated from their TeX output boundary, documented
  public API, independent input tables, finite numeric-key validation, consistent
  formatting and namespaced private TeX macros.
- English reference manual with indexed keys, defaults, units, executed examples,
  numerical limitations, Lua API and the full 246-record material catalogue.
- English examples and technical notes, automated format comparisons and local
  source/TDS packaging with provenance and installation verification.

### Migration from the prototypes

Install the `tex/` tree instead of copying files from the repository root.
The documented TikZ path interface and `\femmplot` syntax are preserved.
The original polygon commands were retained in 0.5 and 0.6, then removed
in 0.7. Use `femm/region` on native TikZ paths.

Private macros are now under `\luafemm@...`; old undocumented implementation
macros are not compatibility interfaces. Numeric literal keys reject non-finite
values or executable Lua text; keys explicitly documented as PGF expressions
continue to accept expressions. Lua model constructors and declarations copy
input tables, so changing a caller's table no longer changes an existing model.

## Prototype history

- **0.4:** FEMM material catalogue and named straight/curved PGFPlots profiles.
- **0.3:** alnico/soft-iron permanent magnet, local mesh density, filtered exact
  orientation/incircle predicates and mesh stress tests.
- **0.2:** native TikZ paths, transformations, curved geometry and constrained
  Delaunay meshing written in Lua.
- **0.1:** nonlinear planar P1 solver, oriented potential contours and a U-core
  electromagnet through an initial generic command interface.
