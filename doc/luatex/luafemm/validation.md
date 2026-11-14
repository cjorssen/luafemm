# Validation scope

Version: 0.11.0-dev. Runtime cache revision: 13.

Observed execution on 2026-11-12 used macOS 27.2.1 on arm64, TeX Live 2026,
LuaTeX/LuaHBTeX 1.24.0, PGF 3.1.12, PGFPlots 1.18.3 and
ConTeXt MkIV 2026.09.10 with the LuaTeX engine. These are observed versions,
not a claim that all other engines/platforms have been tested.

The numerical regression comprises 26 Lua suites covering analytic/manufactured
fields, nonlinear material laws, permanent magnets, meshing/topology, boundaries,
profiles/integrals, static forces, harmonic diffusion/circuits and sheet losses.
A separate rational-arithmetic oracle checks geometric predicates. Invalid TeX
inputs and operation order are checked. Cold, warm and shared frozen calculations
agree across plain LuaTeX, LuaLaTeX and ConTeXt MkIV to the tested 1e-10 tolerance.
The isolated installation check copies runtime files into a private TEXMF,
compiles consumers outside the checkout, records loaded TeX/Lua paths and fresh
PDFs, and checks failure for missing siblings even with a complete shadow copy.

The manual contains 44 complete source/result examples. Its 41 LuaLaTeX examples
execute directly within the manual; the three plain/ConTeXt document-wrapper
exceptions are explicitly identified, rebuilt and checked in their native formats.
Reference fragments execute in supplied model/axis/loading contexts. The README
example is also compiled. The index, cross-references and rendered figures are
checked; compilation alone is not evidence of scientific accuracy.

Independent material-loss comparisons performed on 2026-10-22 used xfemm revision
`34ac6cc586a346066984526620dfb67383c8e3c5`: 15 samples across intrinsic complex
anisotropy, ordinary/deep-skin sheets, axisymmetry and a nonuniform core. The
largest relative differences were 0.255% for potential, 2.898% for B/H, 0.280%
for energy and 1.304% for magnetic loss. These comparisons apply to those cases,
not to all magnetic physics. They were not rerun for documentation/layout changes.

The general regression and optional native oracle harnesses are development
infrastructure and are not required to install or rebuild the public manual.
Users can reproduce the displayed examples with the included documentation
builder, compare physical limiting cases and refine meshes/domains for their
own problem. Windows FEMM GUI acceptance, all-platform support and full FEMM
parity are not certified. Nonlinear harmonic laws and detailed wire models remain
outside the documented supported subset.
