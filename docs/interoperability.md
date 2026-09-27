<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# FEMM interchange implementation and validation

The 0.7.0-dev implementation follows the seven work packages in
[femm-interoperability-audit.md](femm-interoperability-audit.md), preserving
native TikZ paths and a pure-Lua runtime. The English manual is the public key
reference; this document records implementation choices and reproducible evidence.

## Work packages

1. **Specification and fixtures.** `luafemm-fem.lua` owns the passive magnetic
   4.0 document schema and codecs. The authored `tests/fixtures/uniform.fem`
   supplies an independent affine-field reference. Parser tests cover metadata
   round trips, indexing, truncation, units and rejected physics. Unknown scalar
   fields survive passive parsing but unknown physical fields are rejected by
   the solver adapter.
2. **Paths, topology and domains.** Native evaluated PGF soft-path data are retained
   with their frame and approximation tolerance. `luafemm-topology.lua` splits
   intersections/shared edges and constructs oriented bounded faces. Imported
   arcs are subdivided by angular and chord bounds. The existing constrained
   triangulator meshes a rectangular staging domain; carving then removes
   excluded/outside triangles, compacts vertices and audits active area, exposed
   edges and Euler characteristic. Consequently staging refinement may do extra
   work outside an imported domain. Shared child boundaries count as a union,
   not independent holes. Disconnected pure-Neumann components have separate
   compatibility checks and gauges.
3. **FEM output.** Native paths compile into a polygonal arrangement without
   generating triangles. Export retains computational interfaces and creates
   effective materials and verified interior labels per final face. Region
   current density is included in each effective material. The writer does not
   attempt to recover separate winding turns/current from a native NI input.
4. **FEM input.** The adapter supports isotropic planar DC materials, constant
   magnetisation, series excitation, default labels, circular arcs, true holes
   and supported tagged boundaries. Source metadata are copied and retained;
   drawing is independent of physical registration. Imported physical topology
   is immutable in this release, with explicit errors for region/side overrides.
   The passive serializer can separately retain the original arcs/properties;
   the computational exporter emits the domain used by Lua.
5. **Nonlinear fidelity.** `interpolation=femm` selects an independently written
   natural cubic solver, the monotonicity repair used by the audited static
   material implementation, and endpoint-slope extrapolation. Input B-H data
   remain owned and unmodified. The existing piecewise-linear law stays the
   native default. Derivative tests and independent postprocessor H samples
   exercise this distinction. Nonlinear ANS output requires the compatible law.
6. **ANS output.** `luafemm-ans.lua` writes physical nodal A and the actual mesh,
   with labels, current records and empty periodic/air-gap trailers. It rejects
   missing/stale solutions and never remeshes. Cache format 5 includes the new
   topology/physics descriptors and audits carved domains before restoration.
7. **Generic interface and documentation.** Import/export families, format and
   interpolation choices, named delayed exports and ordinary geometry/group
   styles are tested under plain LuaTeX, LuaLaTeX and ConTeXt MkIV. The three
   interchange examples include source/output pairs in the manual: noncircular
   curves, a U electromagnet node and a U magnet with an alnico section. Each
   exports and imports its computational problem.

The initially proposed `luafemm-document.lua` is folded into the passive FEM
codec rather than becoming a separate module with duplicated record ownership.
The public API is documented in `doc/sections/interchange.tex`.

## Reproducing the independent checks

Run from the repository root:

```sh
make test
python3 scripts/interop.py --build --xfemm ../xfemm
make dist
```

The optional native runner copies the checkout into `build/interop/xfemm-source`,
configures CMake and compiles the magnetic reader/mesher/solver/postprocessor.
CMake may download the checkout's Tangle dependency. Native tools never become
installed package requirements and the original sibling checkout is untouched.
On macOS the disposable copy needs three compatibility adaptations: use
`__sincos`, replace obsolete `ptr_fun` trimming expressions with lambdas, and
replace `malloc.h` includes with `stdlib.h`. The runner applies these explicitly.

`tests/xfemm_probe.cpp` loads the actual emitted answer file through
`FPProc::OpenDocument`. The audited xfemm Lua `open()` command rejects ANS files,
so using that command would not test the required disk-reader path. The probe
turns field smoothing off and returns A, Bx, By, Hx and Hy at specified points.

Development evidence on 2026-09-27:

- xfemm revision `34ac6cc586a346066984526620dfb67383c8e3c5`;
- Apple Clang 21.0.0, macOS arm64;
- 72 nonlinear material samples from Lua-written U-device ANS files agree with
  the independent disk postprocessor to a scaled tolerance of `1e-7` for A/B/H;
- four independent xfemm solves of exported problems meet a 2% relative By
  threshold at the selected gap sample; two Lua mesh resolutions are included;
- uniform-field rectangle, circular annulus and restored annulus cache pass
  independent affine-field checks;
- constant real mixed boundary data are checked by both ANS reading and a fresh
  xfemm solve;
- cold, warm and frozen computations, imports and exports compile in all three
  TeX formats, with numerical agreement at `1e-10` in the format matrix.

Selected independent gap samples at (-32.5, 31.5) mm:

| Device | Lua spacing (mm) | Lua By (T) | xfemm By (T) |
| --- | ---: | ---: | ---: |
| U coil | 6 | 0.6532007402 | 0.6534732801 |
| U coil | 3 | 0.6516734450 | 0.6534732801 |
| U alnico | 6 | 0.4327833082 | 0.4286297327 |
| U alnico | 3 | 0.4270283392 | 0.4286297327 |

These are specific draft-mesh comparisons, not a universal numerical error bound
or a claim of monotonic convergence. FEMM mesh settings differ from Lua spacing.
Reports, input hashes, compiler/revision metadata and logs are generated under
`build/interop/`. All fixture geometry and the probe source are original project
work; upstream fixture collections and native binaries are not redistributed.

## Deliberate limits

Windows FEMM GUI acceptance remains unverified on the macOS host. Compatibility
claims are limited to the audited source conventions and actual xfemm tests.
AC, axisymmetry, previous-field modes, periodic/sliding-gap boundaries, point
properties, anisotropy, unsupported laminations, parallel circuits and executable
magnetisation expressions remain outside the solver subset. Open/isolated graph
constraints and ambiguous geometry are diagnosed. Imported filenames and names
are passive byte strings; no implicit Windows code-page conversion is performed.

The `.ans` file is a solution exchange, not a complete luafemm cache or a TeX
source round trip. Analytic Bezier/ellipse export is unnecessary: the polygonal
interfaces used by the solver are the exchange contract. The manual explains
geometry tolerance separately from element refinement and shows an overlay of
the actual computational interfaces.
