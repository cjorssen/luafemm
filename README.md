# luafemm

**Version 0.11.0-dev — 2026-11-12**

A generic LuaTeX library for **TikZ, PGF and PGFPlots**: describe magnetic
materials on ordinary paths, solve in Lua, and draw field lines, meshes,
quantitative profiles and integral results from the same model. It works with
**plain LuaTeX, LuaLaTeX and ConTeXt MkIV**, without a native solver or shell escape.

[Read the manual](doc/luatex/luafemm/luafemm-manual.pdf) ·
[Report a bug](https://github.com/cjorssen/luafemm/issues)

| Curved iron and a tilted armature | Slotted machine and its air-gap field |
| :---: | :---: |
| [![Computed magnetic field and mesh of a U electromagnet with a tilted armature.](images/tilted-u.png)](examples/curved-core.tex) | [![Computed field in a slotted stator surrounding a smooth rotor, with an angular gap profile.](images/slotted-machine-field.png)](examples/tutorial-machine-sine.tex) |

These are actual figures extracted from the manual. Click a preview for the
complete example. The machine tutorial progressively distributes slots and
winding currents to approach a sinusoidal gap field.

## First picture

```tex
\documentclass{article}
\usepackage{tikz,pgfplots}
\usetikzlibrary{femm}
\usepgfplotslibrary{femm}
\pgfplotsset{compat=1.18}
\begin{document}
\begin{tikzpicture}[scale=.65,
  femm/problem={xmin=-85,xmax=85,ymin=-65,ymax=85,mesh size=4}]
  \node[femm/u electromagnet={gap=3,ampere turns=6000}] (M) {};
  \path[femm/profile={name=gap,samples=101}]
    (M.left pole)--(M.left armature);
  \pic[femm/lines=19] {femm field};
\end{tikzpicture}
\par\medskip
\begin{tikzpicture}
\begin{axis}[width=10cm,height=5cm,
  xlabel={Distance (mm)},ylabel={Field (T)},grid=major]
  \femmplot[profile=gap,component=B,plot options={blue,thick}]
\end{axis}
\end{tikzpicture}
\end{document}
```

Compile with `lualatex`. Begin with a coarse mesh; refine the gaps and check
mesh and exterior-domain sensitivity before relying on a quantitative result.

Declare **all materials, physical paths, component nodes, circuits and boundaries
before the first mesh, field or measurement request**. A profile declaration does
not solve; plotting it does. The manual explains operation order and scope.
Options use PGF keys and ordinary `.style`/`.append style` definitions; TikZ paths,
transformations, nodes and drawing options retain their usual roles.

## Supported physics and figures

- Nonlinear static B–H curves, permanent magnets, linear anisotropy, planar and
  axisymmetric formulations, tagged and periodic boundaries, and a circular
  exterior approximation within their documented restrictions.
- A catalogue of 246 FEMM magnetic and nonmagnetic material records. Having a
  catalogue record does not imply support for all its frequency-dependent laws.
- Pure-Lua constrained triangular meshing, local spacing and bounded optional
  angle/area refinement. Arbitrary TikZ curves are flattened with a geometric
  tolerance; the original drawing remains a TikZ path.
- U and E cores, toroids and tapered toroids, permanent-magnet sections, and
  configurable slotted rotor/stator shapes. Ideal confined-flux models preserve
  branch flux and provide mean magnetic paths.
- Linear harmonic eddy currents, skin/proximity redistribution, coupled voltage
  and current ports, series/parallel branches and isolated solid conductors.
  Peak phasors use `exp(+i omega t)`.
- Complex linear permeability and insulated-sheet effective laws, with separate
  intrinsic magnetic, microscopic sheet and resolved Joule losses.
- Static operating-point captures, frozen/incremental responses, small harmonic
  responses about a fixed DC tangent, and independent static moving-band positions.
- Field contours and current views, scalar/vector plots, curved profiles and
  spatial harmonics, selections and integrals, static force/torque methods,
  circuit observables and PGFPlots parameter series.
- Passive versioned calculation caches and the documented static FEM/ANS subset.
  FEM export approximates arbitrary TikZ curves by controlled straight segments.

These capabilities have explicit combination and geometry limits. Large-signal
nonlinear harmonics, general harmonic force/torque observables, detailed wire
models, public harmonic FEM/ANS exchange and full FEMM parity are not supported.
The manual distinguishes the axial long-conductor reduction from the planar and
axisymmetric potential formulations. This is not a complete Triangle or Tangle
port. A small algebraic residual does not establish discretisation accuracy.

## Examples and manual

The English manual uses PGF's documentation environments, with guided tutorials,
full key paths/defaults, units, operation order, numerical methods and limitations.
Its LuaLaTeX examples execute from their displayed sources inside the manual;
three explicitly identified native plain/ConTeXt demonstrations are compiled
separately in their actual formats.

| Example | Question |
| --- | --- |
| [U electromagnet](examples/tutorial-coil.tex) | How do winding current, field lines and gap profiles connect? |
| [Permanent magnet](examples/tutorial-alnico.tex) | How does an alnico section drive a soft-iron return path? |
| [Slotted machine](examples/tutorial-machine-sine.tex) | How can distributed currents shape the angular gap field? |
| [Nondestructive testing](examples/tutorial-ndt-scan.tex) | How do groove depth and probe position change coil impedance? |
| [Sheet losses](examples/harmonic-materials.tex) | How does sheet thickness affect losses versus frequency? |
| [Static forces](examples/forces.tex) | How does armature gap affect attraction? |

## Installation and dependencies

Install the contents of `tex/` into a TEXMF tree, keeping the paths:

```text
tex/luatex/luafemm/          TeX libraries and format loaders
tex/luatex/luafemm/lua/      Lua modules
```

Keep all files from one version together. Missing components produce an error;
modules from another installation are not used to complete a partial package.
Refresh the tree's filename database if its configuration requires it.
The optional TDS ZIP can be extracted at a TEXMF root.

In plain LuaTeX, load `tikz.tex` and use `\usetikzlibrary{femm}`; load
`pgfplots.tex` and `\usepgfplotslibrary{femm}` for profiles. In LuaLaTeX the
same libraries work after loading `tikz,pgfplots`; `\usepackage{luafemm}` is
also available. In ConTeXt MkIV, use `\usemodule[luafemm]`, and load the
`pgfplots` module/library for plots. Select LuaTeX explicitly when the default
ConTeXt engine is LuaMetaTeX. Complete wrappers are in [examples/formats](examples/formats).

Runtime dependencies are LuaTeX and the chosen format's PGF/TikZ/PGFPlots
installation. Neither FEMM, xfemm, Triangle, Tangle, Python nor a compiler is
required to use the library. Observed platform and engine versions are recorded
in [validation notes](doc/luatex/luafemm/validation.md); Windows FEMM GUI
acceptance is not claimed.

To rebuild the manual from the browsable source archive, run from its root:

```sh
python3 doc/luatex/luafemm/build-manual.py
```

This needs Python 3, LuaLaTeX, plain LuaTeX, ConTeXt MkIV, the installed PGF manual
macros, the LaTeX packages named in the manual preamble, MakeIndex, and
pdfcrop/Ghostscript for the three native-format demonstrations. It writes only
under `build/`. It does not download software or require shell escape.
The requested CTAN category is `/macros/luatex/generic/luafemm`; this archive
category is distinct from the installed TDS paths and TeX Live import rules.

## Author, licence and validation

**Christophe Jorssen <christophe.jorssen@gmail.com>** is the author and Current
Maintainer; LPPL maintenance status: **maintained**. Original code and documentation
use **LPPL 1.3c or later**. See [LICENSE](LICENSE), [LICENSES.md](LICENSES.md) and
[MANIFEST.md](MANIFEST.md). FEMM catalogue data and their generated derivatives
retain separate upstream terms, including the Aladdin Free Public License.
They are not relicensed as LPPL.

The implementation draws on studying FEMM/xfemm and Triangle formulations and
algorithms, reimplemented in Lua/LuaTeX. No native solver/mesher source or binary
is included. See the retained data provenance and notices.

The package was **vibe-coded with OpenAI Codex and the Astra-6 model**, under
Christophe Jorssen's direction. Neither AI assistance nor publication certifies
scientific correctness. Independently validate computed results and material
assumptions, including mesh and domain convergence.
