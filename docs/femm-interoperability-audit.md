<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# FEMM problem and solution interchange: audit and proposed work plan

Date: 2026-09-27. Baseline: luafemm `3783c6934889d057e315570d134a3a1dbfc96d53`
(`v0.6.0-dev`). This is a source audit and a design proposal, not an implemented
interchange API or a claim of experimentally verified file compatibility.
The proposed next development series is `0.7.0-dev`.

Implementation update: the seven work packages have now been implemented for
an explicit planar DC subset. See [the implementation and validation record](interoperability.md)
and the manual's interchange chapter for the actual API and limits. The design
below is retained as the pre-implementation audit; proposed module names and
features should not be mistaken for a promise of broader current support.

## Scope and source evidence

The requested operations are:

1. Load a magnetic `.fem` problem created in the FEMM preprocessor.
2. Write a TikZ-defined magnetic problem as a `.fem` that FEMM can mesh and solve.
3. Write luafemm's converged finite-element solution as an `.ans` that FEMM can
   open for post-processing, without substituting a new FEMM solve.

Runtime code remains generic TeX and pure Lua. FEMM, xfemm and Triangle are
algorithm/format references and independent development-time test tools.
The existing `.lfc` cache remains a separate, validated computation cache.

**Design requirement:** FEMM interchange must preserve TikZ's geometric freedom.
Do not restrict native material paths to straight segments and circular arcs.
Keep the supported evaluated PGF path vocabulary, transforms, compound contours,
fill rules and node/pic composition. The exchange format is an output adapter,
not the definition of the native modelling language. This requirement does not
claim that every TikZ drawing operation already defines a valid physical region:
material contours must still be closed and satisfy the geometry checks.

The sibling repositories were clean when audited:

| Repository | Revision |
| --- | --- |
| `../FEMM` | `7d9e8ed90772ddb939a3256f976e231bcb5cca97` |
| `../xfemm` | `34ac6cc586a346066984526620dfb67383c8e3c5` |
| `../Triangle` | `d3d0ccc94789e7e760f71de568d9d605127bb954` |

The following locations are primary evidence; paths are relative to those
repositories, and line numbers refer to the audited revisions.

| Source | What it establishes |
| --- | --- |
| FEMM `femm/FemmeDoc.cpp:2430`, `OnSaveDocument` | Magnetic format 4.0, problem/property sections, nodes, segments, arcs, holes and block labels |
| FEMM `fkn/prob1big.cpp:755`, static solution writer | `.fem` preamble, solution nodes/elements, per-label current data, periodic and air-gap trailers |
| FEMM `femm/FemmviewDoc.cpp:978` | What the magnetic postprocessor actually reads; element-to-label mapping and optional solution columns |
| FEMM `fkn/prob1big.cpp:354` and `:619` | Mixed-boundary assembly and prescribed-potential coordinate convention |
| FEMM `fkn/femmedoccore.cpp:1070` | Expansion of series circuits into label-specific excitation using signed turns |
| FEMM `fkn/matprop.cpp:37` and `:242` | B-H slope construction, cubic interpolation and endpoint-slope extrapolation |
| xfemm `cfemm/libfemm/FemmReader.cpp` and `FemmProblem.cpp` | Independent problem reader/writer and optional-field handling |
| xfemm `cfemm/libfemm/CBlockLabel.cpp:110` | One-based property references, diameter-like mesh parameter, default/exterior flag bits |
| xfemm `cfemm/fsolver/static2d.cpp:1042` | An `.ans` writer with four-column static element records |
| xfemm `cfemm/fpproc/fpproc.cpp:1187` | Independent postprocessor reader for compatibility testing |
| Triangle `src/Triangle/triangle.h:65` and `triangle.c:6853` | Tagged segments, hole seeds, region attributes and flood-based domain carving |

A header inventory of the 28 `.fem` fixtures in the sibling source trees found
28 planar DC problems, 17 files with arcs and 3 with unmeshed holes. Units were
metres (21), inches (3), centimetres (3) and millimetres (1). Many files use
periodic or air-gap conditions. These are useful positive and rejection cases,
but do not cover the entire format; axisymmetric/AC rejection and all six length
units need additional authored fixtures. No native compatibility tests were run
as part of this audit, and no upstream fixture has been copied into this project.

## Findings: geometry and topology

FEMM stores a planar graph, not a sequence of filled polygons. Nodes, straight
segments and circular arcs delimit faces; point-like block labels assign material,
mesh settings, circuits, groups and magnetisation to those faces. Hole seeds mark
faces that must not be meshed. A default block label can supply otherwise unlabeled
regions. Visible/hidden segment flags are drawing metadata, not permission to
delete a physical constraint.

Luafemm currently stores ordered material contours with nonzero/parity fill rules.
The last filled region wins, and the rectangular exterior is filled with air by
default. `luafemm-tex.lua:145` captures flattened paths and rounds coordinates;
the high-level TikZ node/pic or original Bezier construction is not retained.
The mesher constructs its rectangle at `luafemm-mesh.lua:118`; boundary assembly
recognises its four sides; mesh and cache audits assume rectangular coverage and
Euler characteristic one. A material cavity is currently still part of the mesh.

Consequently:

- A circular FEMM exterior must not silently become a rectangular air box.
- A FEMM `No Mesh` hole must not silently become an air region.
- Exporting one block label per original TikZ path is incorrect for overlaps,
  compound contours and disconnected filled pieces. Labels belong to final faces.
- A polygon centroid is not a reliable interior label for a concave region or
  a region with holes. Use a verified interior point with boundary clearance.
- Importing segments/arcs and guessing that each loop is a standalone material
  polygon loses shared interfaces, region connectivity and default-label semantics.
- Missing, duplicate or conflicting block labels require explicit diagnostics;
  apply default-label semantics deliberately rather than guessing a material.

Introduce a common planar topology layer. Preserve the semantic input document,
then compile it into split, tagged constraints, oriented face boundaries and
explicit material/excitation assignments. Resolve intersections, shared sides,
T junctions, nested cycles and TikZ painting precedence before meshing. Preserve
source IDs through splitting and refinement. A face ID must remain available on
every final triangle, separately from its material name.

For imports, retain circular arcs in the source document and flatten them only
for computation, with a documented chord-error bound and their FEMM angular mesh
constraint. Native TikZ curves are exported through the common polygonal
computational geometry described below. An unchanged imported document can retain
its original arcs in a source-preserving FEM export, but this is distinct from
exporting the exact polygonal problem solved by luafemm.
Load imported coordinates directly into the Lua document and topology layers.
Drawing them with TikZ and recapturing the resulting paths would introduce the
frontend's rounding and snapping into an otherwise lossless import. Rendering
must remain a separate view of the imported geometry.

Generalised domain support must update meshing, sampling, boundary assembly,
refinement and cache validation together. For a mesh with C connected components
and H unmeshed holes, the relevant audit is `V-E+T=C-H`. Boundary tags must survive
all subdivisions. Pure-Neumann compatibility and potential gauges must be checked
per connected component, not only for the global current sum.

### Preserving TikZ paths through computation and interchange

The current solver already works on polygonal approximations, not analytic TikZ
curves. `luafemm-path.lua` reads evaluated PGF lines and cubic Beziers and uses
adaptive de Casteljau subdivision before meshing. The existing problem/region
`curve tolerance` keys control that approximation. Consequently a noncircular
curve need not be rejected merely because FEMM cannot store it analytically.

Separate three representations:

1. **Source geometry:** evaluated PGF primitives in model coordinates, with
   cubic control points, contour closure, fill rule and property provenance;
   or imported FEMM segments/arcs and labels. Retain this data before flattening
   so geometry can be regenerated at a new tolerance. This is a planned extension:
   0.6 currently retains the flattened contours, not their cubic control points.
   It preserves evaluated geometry, not arbitrary TeX macros or their source text.
2. **Computational geometry:** one canonical, tagged polygonal arrangement,
   compiled from the source using explicit tolerances. It contains final faces,
   interfaces, exclusions and boundary assignments and feeds both meshing and
   the default FEM export. Do not independently approximate curves in each writer.
3. **Solution mesh:** the refined triangles and physical nodal A, tied to the
   computational geometry by stable IDs and a revision/fingerprint. Preserve this
   mesh exactly for ANS output; its internal triangle edges are not FEM geometry.

The export contracts are:

| Operation | Geometric contract |
| --- | --- |
| FEM import | Preserve supported source segments/arcs and labels; compile them into the same computational representation as TikZ paths. |
| Default FEM export | Write the canonical computational interfaces as straight segments, including polygonal approximations of Beziers, ellipses and other evaluated curves. This is exact relative to that polygonal model, within numeric serialization precision, and approximate relative to the source curve. |
| Source-preserving FEM re-export | Retain imported arcs only where the source remains representable and consistent after edits. Treat this as a separate, explicit policy; never silently switch geometry during a paired FEM/ANS export. Analytic export of arbitrary TikZ curves is not promised. |
| ANS export | Embed the polygonal problem corresponding to the solved snapshot, then its existing mesh and A. Use the same problem snapshot for a paired default FEM export. A curved TikZ source requires no extra geometric approximation at this stage. |

Do not fit circular arcs automatically to arbitrary Beziers: this adds another
approximation and complicates topology and tolerance guarantees without being
necessary for interoperability. A future analytic-preservation option may retain
proven circular primitives; an affine transform can turn a circle into an ellipse,
and evaluated PGF cubic data alone do not establish exact circular provenance.

Keep geometric and finite-element resolution independent. Reuse the existing
problem/region `curve tolerance` keys and their model-unit convention; `mesh size`
and refinement keys control elements inside the resulting domain. Finer triangles
cannot recover curvature discarded earlier. Changing geometry tolerance after a
solve must invalidate geometry-dependent mesh/solution caches and require a new
solve before ANS output. It must never silently produce a different FEM preamble
for an existing solution. Extend cache snapshots with the source/compiled geometry
and approximation policy needed for reproducible export after restoration.

Audit rounding, snapping, intersection splitting and subdivision limits together
with the chord tolerance. Do not advertise the chord tolerance alone as a bound
on the complete geometry pipeline. Preserve shared interfaces once, retain source
tags during subdivision, and detect collapsed gaps, ambiguous contacts or changed
face assignments. Report unresolved geometry instead of silently deleting small
regions. In particular, test narrow curved air gaps at several geometry tolerances.

The manual must distinguish the smooth TikZ drawing from the polygonal domain
actually solved. Provide a style-controlled diagnostic overlay of computational
interfaces and optionally the mesh, so users can inspect the approximation.
Start tutorials with modest geometric and mesh resolution, then tighten each
separately and compare the measured field as the model approaches its final form.

## Findings: numerical and format semantics

| Item | Audit result and required handling |
| --- | --- |
| Indexing | Geometry endpoint IDs are zero-based. `.fem` property references are one-based with zero meaning absent. `.ans` triangle label indices are zero-based. Convert only at the codec boundary. |
| Length units | Six units exist. Normalize deliberately; preserve the declared file unit and convert `x`, `y`, depth, arc geometry and area constraints consistently. The TikZ frontend uses model millimetres; solver nodes use metres. |
| Coordinates flag | `[Coordinates]=polar` does not mean stored node pairs should blindly be converted from polar to Cartesian. The static boundary polynomial uses this flag to select Cartesian or polar variables. |
| Current density | FEMM material `J_re` is in MA/m²; luafemm uses A/m². The conversion is a factor of one million. Real/imaginary fields must not be silently conflated. |
| Block mesh setting | The fourth block-label field is converted by xfemm to `area=pi*d*d/4`. It is not directly luafemm's point spacing or an area in SI units. Preserve the original control and distinguish its numerical interpretation. |
| Mesher/solver controls | FEMM `MinAngle`, smart-mesh settings and `Precision` are not promises of identical luafemm algorithms or stopping criteria. Retain file metadata and report how it is applied; do not overwrite it with library defaults just because their keys were not explicitly set. |
| Depth | Luafemm currently has no model depth. It does not affect the planar A/B solve, but FEMM uses it in integrated quantities. Preserve imported depth and require an explicitly specified, meaningful depth for native exports. |
| Materials | Imported custom definitions are authoritative; do not replace them by a catalogue record with the same name. Keep raw properties alongside the effective supported solver law to avoid applying fill factors twice. |
| Region excitation | TikZ paths can override current density per region. FEMM stores prescribed density in material properties. Native export may need distinct effective material records for equal magnetic laws with different currents. |
| Circuits | Store original circuit identity, current and signed turns. Resolve series-current constraints by labeled-region area, with tests for total current and return sign. Do not blindly add `NI/S` to a material source: FEMM's circuit calculation accounts for the prescribed material contribution. General parallel/conductive circuits need separate support. |
| Node components | Native U components currently provide NI and flatten into regions. Actual I and N cannot be recovered independently from NI. An equivalent current-density export must be identified as such; imported circuit metadata can be retained separately. |
| Boundary placement | FEMM attaches properties to individual segments/arcs and sometimes points. The current four-side interface is insufficient. Introduce tagged-edge boundary data; keep side keys as a convenient rectangular adapter. |
| Boundary defaults | An unassigned FEMM boundary must retain its natural condition; it must not acquire luafemm's default outer `Az=0`. Missing/default properties and pure-Neumann gauges need explicit tests. |
| Dirichlet data | FEMM evaluates `A0+A1*x+A2*y` in file coordinates, including phase projection in DC. SI-based luafemm gradient coefficients therefore need a length conversion; polar-varying data need additional support or rejection. |
| Mixed data | Constant real FEMM type-2 data are candidates for `c0=c`, `c1=g` in `Ht=c*Az+g`, consistent with the audited weak-form signs. Prove this with analytic edge tests, including coercive material touching the boundary. Negative coefficients exceed the current solver contract. |
| Varying Neumann/Robin data | Luafemm accepts affine g; FEMM's mixed c1 is constant per property. General affine data along an edge cannot be represented exactly by assigning one constant to each subsegment. Strict export should reject this case; a later explicit approximation mode would need a tolerance and tests. |
| Point sources/constraints | They are real FEMM features, not decoration. Preserve in the document and reject solving until the corresponding nodal semantics exist. |
| Magnetisation expressions | FEMM may store Lua expressions such as `theta+180`. Parse and retain them as data, never execute imported text. Initially accept constant directions for solving; expression support would need a deliberately restricted evaluator. |
| Native Lua callbacks | Arbitrary source or boundary functions cannot generally be represented by FEMM's finite property records. Exact export must reject unsupported callbacks rather than silently sample them into different physics. |
| Previous solutions | Preserve previous-solution references and incremental/frozen-permeability settings as metadata. Reject solving these modes until their dependence on the prior field is implemented; the first ANS writer targets ordinary non-incremental DC solutions. |

### Nonlinear constitutive laws are a compatibility issue

`luafemm.constitutive` uses piecewise-linear H(B) and a vacuum slope beyond the
last B-H point. FEMM uses cubic interpolation with computed slopes and extrapolates
using the last slope. Identical point tables therefore do not guarantee identical
H(B), differential permeability, energy or a nonlinear recomputed solution.

This matters particularly for `.ans`: the postprocessor reconstructs B from nodal
A, then evaluates H and integrated quantities using its own material law. A file
can be syntactically valid and preserve A/B while producing different H or energy.
Make constitutive-law compatibility a named workstream, not an unexplained test
tolerance. Preserve the existing law for existing documents; design an explicit
FEMM-compatible interpolation/extrapolation option and validate its derivatives,
monotonicity, saturation and permanent-magnet shift before claiming nonlinear
post-processing equivalence. If unavailable at an intermediate milestone, state
the narrower A/B guarantee explicitly and do not certify H/energy equivalence.

### The static `.ans` layout

The writer should use the same computational problem snapshot as the default `.fem` writer,
then append `[Solution]` and these records:

1. Mesh-node count; coordinates in declared file units, real Az in T m and the
   appropriate boundary marker. Do not reuse FEMM's internal centimetre scaling
   or its solver's scaled unknown: luafemm already stores physical Az.
2. Triangle count; three zero-based node indices and a zero-based block-label ID.
   Modern FEMM writes additional edge tags and a previous-current field; the
   audited normal static readers consume the leading four fields. Select and
   document a canonical non-incremental dialect, verified in both readers.
3. A current/circuit record for every block label, not merely for every named
   circuit. For an equivalent-material export without circuits, the upstream
   dummy record `1 0` means zero additional current; avoid double counting J.
4. Periodic-pair count and air-gap-element count, both explicitly zero for the
   supported nonperiodic subset. Do not omit these trailers and rely on EOF.

The exact original converged mesh and A values must be written. Remeshing during
ANS export would require projecting the solution and would no longer preserve
the solved field. The writer must also work after restoring a `.lfc` solution;
therefore stable face IDs, effective material data and exchange metadata need
cache descriptors, storage, validation and a new cache-format revision.

## Proposed architecture and public interface

Keep numerical Lua modules independent of TeX:

- `luafemm-document.lua`: owned native/imported source geometry, exchange metadata,
  stable IDs, explicit support/capability reports and validation; no global
  current-model dependency and no FEMM-only restriction on native curves.
- `luafemm-path.lua`: separate evaluated PGF primitive capture from adaptive
  flattening, preserving the existing public path and fill-rule behaviour.
- `luafemm-fem.lua`: passive parser/writer for magnetic format 4.0. Separate file
  syntax acceptance from permission to solve the represented physics.
- `luafemm-topology.lua`: shared planar arrangement and labeled-face compilation;
  reuse robust predicates and extract the current segment preprocessing here.
- `luafemm-ans.lua`: writer from an immutable problem/solution snapshot. A small
  independent test reader can inspect output without promising public ANS import.
- Existing solver, boundary, refinement and cache modules: consume the same
  compiled topology and labels; retain compatibility with existing region calls.
- Generic TeX bridge: only key parsing, lifecycle hooks, string escaping and
  drawing. Plain, LaTeX and MkIV loaders remain thin.

Proposed syntax, deliberately not executable in 0.6:

```tex
\begin{tikzpicture}[
  femm/import={file=u-magnet.fem},
  femm/export={format=fem,file=u-magnet-copy.fem},
  femm/export={format=ans,file=u-magnet.ans}]
  \pic {femm geometry};
  \pic {femm field};
\end{tikzpicture}
```

`format` should be a PGF `.is choice` (`fem`, `ans`); drawing styles should use
ordinary `.style` keys with separate region, interface, group and label styles.
Avoid forcing imported material names into TeX control-sequence names. A later
document-only inspection/drawing mode can preserve unsupported physics without
pretending that it is solvable.

Define lifecycle semantics before implementing these keys:

1. Collect picture options. Read the import once, at model creation, then apply
   only explicitly provided problem overrides. Defaults must not erase file data.
2. Register imported topology and any allowed additional physical declarations.
   Render the imported geometry separately; drawing alone must not register it twice.
3. Compile and freeze the computational geometry and its tolerance policy before
   meshing. Existing declaration-order rules still apply; measurement paths do not
   become physical interfaces. A changed source or geometry tolerance creates a
   new revision and invalidates the old solution for export.
4. Repeated `femm/export` options append output requests bound to this model.
   Process them at the end of the picture after all declarations. FEM output
   compiles geometry but requires no magnetic solve. ANS output requires a
   converged solution, whether freshly solved or restored; otherwise report an
   error instead of starting an unexpected computation or writing empty results.
5. Provide an explicit `\femmexport[...]` using the same backend for delayed
   output, including after a profile outside the source picture triggered solving.
   Specify model selection so starting another picture cannot export the wrong one.

The parser must support CRLF/LF, optional fields, quoted names and the observed
case-insensitive keywords, with line-numbered errors and bounded record counts.
Choose an explicit encoding policy for Windows-origin files and test accented
names and TeX-special characters. Unknown fields can be retained in the document
with a capability diagnostic; unknown physics must not silently reach a solve.
Never use `load`, `dofile` or imported TeX execution to parse a data file.
Resolve input paths through the document's normal search mechanism and output
paths consistently with `--output-directory` and the existing cache. Write a
validated temporary file and replace the destination only after success; do not
accidentally overwrite the input file. Keep numeric serialization locale-neutral
and round-trip-safe, and isolate Lua strings with `\luaescapestring` at the bridge.

## Work packages and acceptance gates

| Stage | Deliverable | Gate before proceeding |
| --- | --- | --- |
| 1. Specification and corpus | Field/index/unit tables, supported-feature matrix, three-level geometry schema, authored minimal FEM/ANS fixtures and malformed-input cases | Native TikZ freedom and polygonal export contracts are explicit; readers agree with fixtures; unsupported AC, axisymmetry, periodic/air-gap, expressions and unknown physics receive explicit diagnostics |
| 2. Path capture, topology and domains | Retained evaluated PGF primitives; reproducible flattening; shared tagged graph/faces; nonrectangular exterior; true exclusions; robust seed placement; generalized boundary and component audits | Areas, connectivity, labels, currents and boundary tags are preserved for overlapping/nested TikZ paths and imported segmented/arced models; transformed Beziers and ellipses remain supported; old 0.6 tests remain valid |
| 3. FEM writer | Export of the canonical computational geometry, equivalent effective materials/excitations, explicit depth and deterministic labels | FEMM GUI and xfemm open and mesh polygonal exports from curved TikZ models; uniform-field and coil cases solve correctly; repeated export does not mesh or solve the luafemm model |
| 4. FEM import and round trips | Import of the declared planar DC subset, faithful metadata retention, arcs, units, labels and series excitations | GUI-created examples draw and solve in luafemm; import-write-reload preserves canonical geometry and physics; nonrectangular and No Mesh cases do not acquire artificial air |
| 5. Nonlinear compatibility | Explicit B-H interpolation/extrapolation policy and FEMM-compatible option | H(B), derivatives and saturation comparisons pass, then soft-iron and alnico-device comparisons converge under independent mesh/domain refinement |
| 6. ANS writer | Embedded computational FEM problem, original mesh/A, label/current tables and zero trailers | Both postprocessors load the file; counts, connectivity, A and unsmoothed B match; tests cover curved TikZ sources, paired FEM/ANS consistency, changed-tolerance rejection, circuits, nonlinear H policy, depth-dependent integrals and cold/warm/frozen caches |
| 7. Generic interface and documentation | Key families, `.is choice`, reusable styles, lifecycle hooks, complete tutorials and compatibility reference | The same import/export/ANS examples compile under plain LuaTeX, LuaLaTeX and MkIV; manual source/output coverage, format tests and extracted TDS installation pass |

The TeX interface should be exercised in small vertical examples throughout
stages 2–6, not added for the first time in stage 7. The first useful vertical
slice is a linear air rectangle with affine prescribed Az: create in TikZ,
export FEM, open/solve in xfemm, import back, solve in Lua and export ANS.
It exposes units, signs and numbering without nonlinear or mesher ambiguity.
Immediately follow that first slice with a noncircular cubic boundary and a
transformed ellipse, before expanding to an oblique interface, a circular exterior,
a real No Mesh hole, nested materials, a series coil and the two U-magnet tutorials.
Include a narrow curved gap and compound transformed node/pic paths. Vary curve
tolerance at fixed fine mesh, then mesh resolution at fixed geometry, to distinguish
geometric error from finite-element error. Compare default FEM geometry with the
ANS problem preamble and solved boundary to detect inconsistent approximations.

The intended first release profile is planar, real magnetostatics with supported
isotropic materials, prescribed density/validated series-current excitation,
constant magnetisation directions and supported tagged Dirichlet/natural/mixed
conditions. Include ordinary circular exteriors and No Mesh holes in the target;
without them, label the importer explicitly as a restricted preview. Preserve but
reject solving anisotropy, unsupported lamination physics, general parallel
conductive circuits, point sources, functional magnetisation, AC, axisymmetry,
periodic/antiperiodic and sliding air-gap models until independently implemented.
Do not make these extra solver features implicit prerequisites for reading or
inspecting a document.

## Validation and documentation strategy

Use three independent checks, since a round trip through the same new reader
and writer can hide reciprocal errors:

- **Structural:** authored fixtures with known indices/units, raw metadata round
  trips, labeled-face membership, domain area, total current, deterministic output,
  malformed/truncated inputs, stale/mismatched solution snapshots and cache reloads.
- **Numerical:** uniform-field affine patches, mixed-boundary signs, per-component
  Neumann compatibility, current-normalized coils, remanence and nonlinear material
  samples. Compare independent solves using convergence criteria rather than
  requiring different meshers to produce the same mesh or nodal ordering.
- **External:** xfemm's file reader, mesher, solver and FPProc, plus actual opening
  and post-processing in Windows FEMM. Load the emitted `.ans` from disk rather
  than only using xfemm's in-memory session shortcut. Compare unsmoothed B away
  from interfaces; smoothing and interface-side selection need explicit settings.

Native tools are optional local development dependencies and an explicit CI
interoperability job, never subprocesses launched by the installed TeX package.
Record their revisions, commands, input hashes and tolerances. Include actual
Windows GUI evidence before claiming GUI compatibility; source inspection and
xfemm acceptance alone are insufficient. Check fixture redistribution terms and
prefer newly authored small fixtures over bundling upstream test collections.

The English manual should explain the three workflows separately, list every
new key/default/unit/scope/error, document import precision, supported physics,
material-law differences, provenance and result-file structure, and state exactly
when files are read or written. Each tutorial must provide runnable TeX, its FEM
input where applicable, and the computed drawing/profile beside or after the code.
Include both directions of the U-coil and alnico examples, the geometry-only
inspection case and a documented rejection example. Clearly separate physical
interchange from reproducing TikZ drawing styles, labels, anchors and source code.
Add a noncircular TikZ-curve export tutorial with the smooth drawing, polygonal
interface overlay and field profile. Document the existing curve/mesh controls
separately, the three geometry representations, cache invalidation and the limits
of geometric fidelity. Retain generic TeX, pure Lua, ordinary PGF keys/styles and
thin, escaped TeX/Lua bridges throughout; the export adapter must not introduce
native runtime dependencies or a second restricted native drawing language.

This audit changes no runtime code, version number or public syntax. Start with
the document specification and the linear vertical slice; the principal workload
is topology/boundary generalisation and constitutive-law fidelity, not text I/O.
