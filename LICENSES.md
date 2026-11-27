<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Licence, authorship and third-party material

## Original luafemm work

Copyright (C) 2026 Christophe Jorssen.

**Author and Current Maintainer:** Christophe Jorssen
<christophe.jorssen@gmail.com>. This work has the LPPL maintenance status
**maintained**. Development: <https://github.com/cjorssen/luafemm>.
Bug tracker: <https://github.com/cjorssen/luafemm/issues>.

The original luafemm code, tests, tooling and documentation may be distributed
and/or modified under the conditions of the **LaTeX Project Public License,
version 1.3c or later**, dated 2008-05-04. The full, unmodified licence is in `LICENSE`;
the [official text](https://www.latex-project.org/lppl/lppl-1-3c.txt) is also
available online. Source headers use `SPDX-License-Identifier: LPPL-1.3c`.
The copyright on the licence text itself belongs to the LaTeX3 Project.

The Work consists of the original files listed in `MANIFEST.md` and their
compiled documentation. The separately identified FEMM data and their generated
derivatives below are an aggregation with the Work, not a grant of LPPL rights
over those third-party materials. In the compiled manual, the material catalogue
appendix retains the data's upstream terms; the original explanatory text is LPPL.
For this generic package, the Base Interpreters include LuaTeX/LuaHBTeX with
plain TeX, LaTeX or ConTeXt MkIV, and the Lua interpreter for the Lua modules.
Development tools additionally use Python 3 and Make.

## FEMM 5 magnetic and dielectric catalogues: MIT terms

The bundled catalogue is pinned to femm-qt revision
`31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3`. Its material-library notice
explicitly licenses the collection, curve fits and data organisation under MIT.
The full permission notice is retained with these separately licensed files:

- `data/femm5/matlib.dat`: unchanged upstream magnetic library (282 records).
- `data/femm5/statlib.dat`: unchanged upstream electrostatic library (95 records).
- `tex/lua/luafemm-dielectrics-data.lua`: generated dielectric transcription.
- `doc/dielectrics-table.tex`: generated dielectric table.
- `data/femm5/LICENSE.txt`: unchanged MIT notice.
- `data/femm5/UPSTREAM.md`: unchanged upstream explanation of material provenance.
- `tex/lua/luafemm-materials-data.lua`: generated transcription.
- `doc/materials-catalog.md` and the manual appendix: generated index.

These data are an MIT-licensed aggregation with the LPPL Work. See
[data provenance](data/femm5/README.md) for pinned sources and checksums.
FEMM 5 solver code is separately GPL-3.0-or-later and is not bundled. Upstream
provenance papers have their own terms and are not included.

FEMM 4 uses the Aladdin Free Public License, with restrictions including paid
redistribution. No separate permissive grant was found for its old material
library. That library and its generated derivatives are therefore excluded from
the curated GitHub, CTAN and TDS distributions. The optional `material edition=femm4`
reader accepts a user-supplied, explicitly identified legacy file; it neither
relicenses nor distributes that file. Published archives contain only the FEMM 5
catalogue. This policy concerns the distributed collection and notices, not an
assertion of copyright over physical facts.

Publication on CTAN is planned, not yet completed. The distribution declares
LPPL 1.3c or later for the original Work and MIT for the bundled catalogue.

## Technical provenance and build dependencies

The work is based on studying the source code of FEMM/xfemm and Triangle, then
reimplementing selected formulations and algorithms in Lua/LuaTeX. No native
C/C++ source or binary from those projects is bundled into the runtime.
This describes the technical approach, not a legal determination that a
language translation removes upstream rights. The mesher has the limitations
stated in the manual and does not claim full Triangle quality guarantees.

The native air-gap module independently constructs Hermite interpolation and
a rectangular energy factorisation; it contains no native expanded matrix
tables. Local validation extracts upstream arithmetic into private build-only
oracles, excluded from the public archives. The native cubic and quintic bases
are identified separately in the interface and manual.

The effective round-wire, square-array and copper-clad aluminium models follow David Meeker's mathematical skin/proximity
formulation and fitted coefficients, audited against FEMM/xfemm `GetFillFactor`
([model description](https://www.femm.info/doku/doku.php?id=proximityloss)).
The Lua module evaluates those equations independently; no native implementation
is distributed. The nonlinear fundamental module projects H(B) and does not
claim equivalence to FEMM's B(H)-based effective curve construction.

The manual uses PGF documentation macros supplied by the user's TeX installation.
Its presentation follows tikz-ext; its original prose and examples are written
for luafemm. Neither PGF's nor tikz-ext's source is redistributed here. Build
tools and TeX dependencies retain their own licences and are not in the archives.
