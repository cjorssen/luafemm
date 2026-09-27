<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# luafemm

**0.8.0-dev** — a generic TikZ and PGFPlots library for planar magnetostatics,
written entirely in Lua and TeX. Declare materials on ordinary TikZ paths,
solve a finite-element model inside LuaTeX, and draw oriented field lines or
sample the field along a path. Works with **plain LuaTeX, LuaLaTeX and
ConTeXt MkIV**. No native solver, external mesher or shell escape is required.

This development version includes nonlinear B–H laws, permanent magnets,
246 FEMM material records, constrained triangular meshing, local mesh spacing,
and PGFPlots profiles of B, H and their components. It is a planar static
solver: it does not implement hysteresis history, eddy currents, axisymmetry
or all of FEMM's material physics. Optional angle/area quality targets are
enforced by bounded Lua refinement, with explicit failure when unattainable.
This is not a complete Triangle port. See the manual for the model and limits.

## Quick start

Install the distribution's `tex/` tree in a TEXMF tree, or use the development
commands below from this checkout.

```tex
\documentclass{article}
\usepackage{tikz,pgfplots}
\usetikzlibrary{femm}
\usepgfplotslibrary{femm}
\pgfplotsset{compat=1.18}
\begin{document}
\begin{tikzpicture}[femm/problem={
  xmin=-15,xmax=15,ymin=-15,ymax=15,mesh size=2}]
  \filldraw[femm/region={material=cast-alnico-5-lng37,
    magnetization angle=90},fill=red!20]
    (-3,-5) rectangle (3,5);
  \pic[femm/lines=13] {femm field};
  \draw[femm/profile={name=probe,samples=101},dashed]
    (-10,8)--(10,8);
\end{tikzpicture}
\begin{tikzpicture}
  \begin{axis}[xlabel={$x$ (mm)},ylabel={$B_y$ (T)}]
    \femmplot[profile=probe,component=By,abscissa=x]
  \end{axis}
\end{tikzpicture}
\end{document}
```

Compile with `lualatex --no-shell-escape document.tex`.
`\usepackage{luafemm}` is a convenience wrapper for the TikZ library;
PGFPlots remains optional.

For plain LuaTeX, load `tikz.tex`, optionally `pgfplots.tex`, then use the same
library commands and `\tikzpicture` / `\endtikzpicture`, `\axis` / `\endaxis`.
For ConTeXt MkIV, use `\usemodule[luafemm]`, optionally
`\usemodule[pgfplots]`, and `\starttikzpicture` / `\stoptikzpicture`.
Run ConTeXt explicitly with `--luatex`. Complete format examples are in
[`examples/formats/`](examples/formats/).

## Ideal magnetic circuits

Choose `field model=ideal` to confine flux to iron, magnets and explicitly
included air gaps. Flux is conserved within every branch and balances at
junctions. This retains finite, potentially nonlinear material laws.

```tex
\begin{tikzpicture}[femm/problem={field model=ideal,depth=10,mesh size=3}]
  \node[femm/e electromagnet={ampere turns=900,coil width=3,
    left gap=1,center gap=2,right gap=4}] (M) {};
  \pic[femm/lines=19] {femm field};
  \pic[femm/mean={component=M,cycle=left,profile=left}] {femm mean path};
  \pic[femm/mean={component=M,cycle=right}] {femm mean path};
\end{tikzpicture}
```

The U shapes accept unequal leg/yoke widths. New shapes include `femm/e core`,
`femm/e electromagnet`, `femm/toroid` and `femm/tapered toroid`, with an optional
permanent-magnet arc. Mean contours can follow section centres or an isolated
half-flux field line. Use `\femmplot[profile=left,component=circulation]` to plot
the accumulated integral of H along a mean contour, or `\femmflux[profile=cut]`
to print a named section's exact discrete flux in webers.

Native paths use `femm/region={material=...,ideal domain=true}`. Excitation seeds
can be declared with `femm/excitation={x=0,y=0,ampere turns=1000}`. The manual
explains declaration order, signs, reference contours and current counting.
See the complete `examples/ideal-*.tex` documents and their rendered results.
Tapered toroids are **planar, constant-depth sections**, not a 3D or axisymmetric
calculation of conical bodies. Ideal models use the Lua Delaunay mesher and
private caches; their floating constraints are not exported to FEMM files.

## FEMM interchange

Version 0.7 adds planar DC `.fem` import, polygonal `.fem` export and solved
`.ans` export in pure Lua. TikZ curves remain available: FEM output uses the
computational polygonal interfaces. Specify `depth` before exporting and choose
`interpolation=femm` before solving a nonlinear problem intended for ANS output.

```tex
\begin{tikzpicture}[femm/import={file=problem.fem},
  femm/problem={mesh size=2},
  femm/export={format=ans,file=problem-lua.ans}]
  \pic {femm geometry};
  \pic {femm field};
\end{tikzpicture}
```

Import supports isotropic DC materials, constant magnetisation, series circuits,
closed regions, circular arcs, true holes and supported tagged boundaries.
Unsupported physics is rejected. Imported physical topology is immutable in
this release; TikZ drawing overlays and field profiles remain available.
See the manual's interchange chapter and `examples/interchange-*.tex`.

Optional independent checks use `python3 scripts/interop.py --build --xfemm ../xfemm`.
The native tools are for development only. Windows FEMM GUI acceptance has not
been tested on the macOS development host.

## Documentation and examples

The [manual source](doc/luafemm-manual.tex) uses PGF's documentation environments,
following tikz-ext's reference-and-example presentation. It includes every
public key and default, all material identifiers, executed examples, units,
scoping, numerical methods, the Lua API and a key index. An opening tutorial
follows a U-core electromagnet and an alnico/iron magnet from drawing to field
lines and straight/curved profiles, with coarse meshes during drafting and
refinement checks before final figures.
Every distributed example appears with its complete source and compiled result,
including the plain LuaTeX, LuaLaTeX and ConTeXt MkIV wrappers. Short worked
examples use PGF's side-by-side layout; longer documents place the output below
the source. `make manual` rebuilds all example PDFs and checks this coverage.

```sh
make manual       # build/doc/luafemm-manual.pdf
make examples     # build/examples/*.pdf
make test         # numerical tests and the three-format integration matrix
make dist         # source and installable TDS archives in build/dist/
```

Development needs Python 3 (standard library) and a TeX Live installation
with LuaTeX, TikZ, PGFPlots, ConTeXt MkIV, the LaTeX packages used by the manual,
MakeIndex and pdfcrop; the manual's margin cropping also needs Ghostscript.
Runtime documents only need their TeX format, LuaTeX and the
selected graphics libraries. The tested versions and regression evidence are
in [validation](docs/validation.md). The build helper handles the explicit
MkIV script launcher needed by some installations with a LuaMetaTeX launcher.

| Example | Contents |
| --- | --- |
| [u-electromagnet.tex](examples/u-electromagnet.tex) | Saturating U core, winding and two air gaps |
| [curved-core.tex](examples/curved-core.tex) | Rounded paths, rotated armature and mesh overlay |
| [mesh-quality.tex](examples/mesh-quality.tex) | U electromagnet before/after angle-based quality refinement |
| [alnico-u.tex](examples/alnico-u.tex) | Alnico / soft-iron magnet with local gap refinement |
| [materials-profiles.tex](examples/materials-profiles.tex) | Catalogue materials, straight and curved field profiles |
| [boundary-conditions.tex](examples/boundary-conditions.tex) | Symmetric half-U, Neumann symmetry plane and gap profile |
| [nested-domains.tex](examples/nested-domains.tex) | Alnico magnet inside an air cavity and a nonlinear iron shell |
| [tutorial-coil.tex](examples/tutorial-coil.tex) | Tutorial: winding, field lines and a profile across both gaps |
| [tutorial-alnico.tex](examples/tutorial-alnico.tex) | Tutorial: permanent magnet, straight and curved profiles |
| [plain.tex](examples/plain.tex) | Native TikZ interface in plain LuaTeX |

## Declaration order and nested domains

Put `femm/problem` and cache settings in the picture options, then declare
**all materials, physical paths, component nodes and boundary conditions before
the first mesh or solve request**. Even `\pic {femm mesh};` freezes the physical
model. Registering a `femm/profile` curve does not solve; plotting it does.
The manual's “The order of operations” subsection lists every trigger, including
cache restoration and begin-picture hooks.

Separate paths paint materials in declaration order: draw the outer region
first, then the inner one. For an explicit material cavity, use a compound path
and native TikZ fill rules, which can be included in a reusable style:

```tex
\tikzset{iron shell/.style={even odd rule,
  femm/region={material=pure-iron},draw,fill=gray!25}}
% Inside a magnetic picture, before any mesh/field request:
\path[iron shell]
  (-20,-15) rectangle (20,15)
  (-10,-6) rectangle (10,6);
```

The hole preserves any previously declared inner material, or air by default.
It remains meshed and can be sampled. This is a material cavity, not an unmeshed
exclusion. `nonzero rule` is also supported, with PGF's contour orientation
semantics. See [the nested-domain example](examples/nested-domains.tex).

## Parametric component nodes

The U core and complete electromagnet are available as genuine TikZ nodes:

~~~tex
\node[femm/u electromagnet={width=80,height=60,leg thickness=15,
  gap=3,ampere turns=6000},anchor=core center] (M) at (0,0) {};
\draw[femm/profile={name=gap},dashed]
  (M.left pole)--(M.left armature);
~~~

Use these inside a picture with `femm/problem`, before solving. A component
owns fixed dimensions and physical anchors; text does not resize its parts.
`femm/u core` creates only the iron core. The electromagnet adds an armature,
two winding sections and optional gap refinement. Winding current density is
computed from the actual section area, preserving the specified ampere-turns.
See [the complete node example](examples/u-electromagnet-node.tex) and the
manual's “Parametric components as nodes” chapter for all keys and anchors.

## Boundary conditions

Set exterior conditions through TikZ options and reusable styles:

```tex
\tikzset{half model/.style={
  femm/boundary={side=left,type=neumann,value=0}}}
\begin{tikzpicture}[half model,
  femm/problem={xmin=0,xmax=85,ymin=-75,ymax=95,mesh size=4}]
  % Declare the half-device, then draw its field.
\end{tikzpicture}
```

The boundary key may precede or follow `femm/problem`; it is applied once
the new model exists. In the picture body, `\tikzset{femm/boundary={...}}`
and `\femmboundary[...]` apply the same condition immediately, before meshing.
Conditions set inside a graphical scope persist in the physical model.
The other sides keep their default zero potential. Choices are `dirichlet`,
`neumann` and `robin`; `side` accepts `all`, `left`, `right`, `bottom` and `top`.
Dirichlet prescribes Az; Neumann prescribes H along the counterclockwise
boundary tangent; Robin prescribes `Ht = coefficient * Az + value`.
Optional `x coefficient` and `y coefficient` add affine spatial terms.
All numerical boundary data use SI units, including coordinates in these terms.
For example, `\femmboundary[y coefficient=.2]` imposes Az = 0.2 y on all sides,
giving Bx = 0.2 T in an empty domain. The manual explains signs, symmetry,
corner compatibility and the potential reference for all-Neumann problems.
See [the half-U example](examples/boundary-conditions.tex).

## Reusing computations

Add a cache option to the magnetic picture:

~~~tex
\begin{tikzpicture}[
  femm/problem={xmin=-20,xmax=20,ymin=-20,ymax=20,mesh size=2},
  femm/cache={file=magnet-field,mode=auto}]
  % Declare the regions, then request field lines or a profile.
\end{tikzpicture}
~~~

`auto` reuses a compatible solution, or just the mesh when physical inputs
change. `refresh` recomputes, `off` leaves files untouched, and `frozen` requires
a compatible saved result. Drawing styles and measurement profiles do not
invalidate the solution. `\femmcachestatus` reports the reuse outcome.

Files use the private `.lfc` format (version 6). FEMM `.fem`/`.ans` interchange
is available separately for ordinary open-air models.
Relative names follow TeX's output directory, if set. Parent directories must
exist, and distinct problems need distinct filenames. The tutorial examples
enable caching. The manual's “Saving and reusing computations” chapter documents
every option, invalidation rule and limitation.

## Repository layout

- `tex/generic/luafemm/`: generic TeX libraries and Lua modules.
- `tex/latex/luafemm/`, `tex/context/third/luafemm/`: thin format wrappers.
- `doc/`: complete English reference manual; `docs/`: supporting technical notes.
- `examples/`, `tests/`, `scripts/`: documents, regression cases and development tools.
- `data/femm/`: unchanged upstream material data and notices.
- `build/`: ignored generated PDFs, logs, caches and local release archives.

See [CONTRIBUTING.md](CONTRIBUTING.md) for coding conventions and release checks,
and [CHANGELOG.md](CHANGELOG.md) for version history and breaking changes.

## Author, maintenance and distribution

Copyright (C) 2026 **Christophe Jorssen**, author and Current Maintainer.  
Contact: [christophe.jorssen@gmail.com](mailto:christophe.jorssen@gmail.com).  
LPPL maintenance status: **maintained**.

Publication is planned on CTAN. Development takes place at
[cjorssen/luafemm](https://github.com/cjorssen/luafemm); please report bugs at the
[GitHub issue tracker](https://github.com/cjorssen/luafemm/issues).

## Development notice

This package was **"vibe coded" with Codex/GPT-6-Astra**, under Christophe
Jorssen's direction. AI-assisted code can contain subtle programming and
mathematical errors even when examples compile and regression tests pass.
Independently check computed fields, meshes and material interpretations against
analytical results or a validated solver, and study mesh and domain convergence.
Do not rely on these results for safety-critical engineering decisions without
independent validation. The current tests do not certify equivalence with FEMM
or Triangle.

## Licensing and provenance

The original luafemm code and documentation are distributed under the
**LaTeX Project Public License, version 1.3c (LPPL-1.3c)**; see [LICENSE](LICENSE),
[LICENSES.md](LICENSES.md) and the file inventory in [MANIFEST.md](MANIFEST.md).

The work is based on studying the source code of FEMM/xfemm and Triangle and
reimplementing selected formulations and algorithms in Lua/LuaTeX. No native
FEMM, xfemm or Triangle source or binary is embedded. The manual identifies
the algorithm references and the differences from those programs.

The imported FEMM catalogue and its generated derivatives retain **separate
upstream terms**, including the Aladdin Free Public License; see
[data/femm/README.md](data/femm/README.md). They are not relicensed as LPPL data.
The current archives therefore contain mixed licensing and must not be described
as entirely LPPL when preparing a CTAN submission. The data's redistribution
status must be settled or explicitly declared before publication.
