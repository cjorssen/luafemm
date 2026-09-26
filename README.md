<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# luafemm

**0.5.0-rc.1** — a generic TikZ and PGFPlots library for planar magnetostatics,
written entirely in Lua and TeX. Declare materials on ordinary TikZ paths,
solve a finite-element model inside LuaTeX, and draw oriented field lines or
sample the field along a path. Works with **plain LuaTeX, LuaLaTeX and
ConTeXt MkIV**. No native solver, external mesher or shell escape is required.

The first release candidate includes nonlinear B–H laws, permanent magnets,
246 FEMM material records, constrained triangular meshing, local mesh spacing,
and PGFPlots profiles of B, H and their components. It is a planar static
solver: it does not implement hysteresis history, eddy currents, axisymmetry
or all of FEMM's material physics. Mesh refinement does not guarantee a
minimum triangle angle. See the manual for the precise model and limits.

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

## Documentation and examples

The [manual source](doc/luafemm-manual.tex) uses PGF's documentation environments,
following tikz-ext's reference-and-example presentation. It includes every
public key and default, all material identifiers, executed examples, units,
scoping, numerical methods, the Lua API and a key index. An opening tutorial
follows a U-core electromagnet and an alnico/iron magnet from drawing to field
lines and straight/curved profiles, with coarse meshes during drafting and
refinement checks before final figures.

```sh
make manual       # build/doc/luafemm-manual.pdf
make examples     # build/examples/*.pdf
make test         # numerical tests and the three-format integration matrix
make dist         # source and installable TDS archives in build/dist/
```

Development needs Python 3 (standard library) and a TeX Live installation
with LuaTeX, TikZ, PGFPlots, ConTeXt MkIV, the LaTeX packages used by the manual,
and MakeIndex. Runtime documents only need their TeX format, LuaTeX and the
selected graphics libraries. The tested versions and regression evidence are
in [validation](docs/validation.md). The build helper handles the explicit
MkIV script launcher needed by some installations with a LuaMetaTeX launcher.

| Example | Contents |
| --- | --- |
| [u-electromagnet.tex](examples/u-electromagnet.tex) | Saturating U core, winding and two air gaps |
| [curved-core.tex](examples/curved-core.tex) | Rounded paths, rotated armature and mesh overlay |
| [alnico-u.tex](examples/alnico-u.tex) | Alnico / soft-iron magnet with local gap refinement |
| [materials-profiles.tex](examples/materials-profiles.tex) | Catalogue materials, straight and curved field profiles |
| [tutorial-coil.tex](examples/tutorial-coil.tex) | Tutorial: winding, field lines and a profile across both gaps |
| [tutorial-alnico.tex](examples/tutorial-alnico.tex) | Tutorial: permanent magnet, straight and curved profiles |
| [plain.tex](examples/plain.tex) | Original command interface in plain LuaTeX |

## Repository layout

- `tex/generic/luafemm/`: generic TeX libraries and Lua modules.
- `tex/latex/luafemm/`, `tex/context/third/luafemm/`: thin format wrappers.
- `doc/`: complete English reference manual; `docs/`: supporting technical notes.
- `examples/`, `tests/`, `scripts/`: documents, regression cases and development tools.
- `data/femm/`: unchanged upstream material data and notices.
- `build/`: ignored generated PDFs, logs, caches and local release archives.

See [CONTRIBUTING.md](CONTRIBUTING.md) for coding conventions and release checks,
and [CHANGELOG.md](CHANGELOG.md) for migration from the prototypes.

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
