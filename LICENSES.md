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
version 1.3c**, dated 2008-05-04. The full, unmodified licence is in `LICENSE`;
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

## FEMM catalogue: separate upstream terms

The following are upstream data, notices or data derivatives:

- `data/femm/matlib.dat`: unchanged FEMM material library.
- `data/femm/LICENSE.txt`: unchanged upstream distribution terms and notices.
- `tex/generic/luafemm/luafemm-materials-data.lua`: generated transcription.
- `docs/materials-catalog.md` and the manual's material appendix: generated index.

These retain their upstream FEMM terms, including the **Aladdin Free Public
License**; they are not relicensed as LPPL. Keep the upstream notice and
provenance with both source and TDS distributions. See `data/femm/README.md`
for the exact source revision and regeneration procedure.

Publication on CTAN is planned, not yet completed. The current archives include
these separately licensed data and must not be described as entirely LPPL.
Their redistribution status must be settled or explicitly declared before
submission; an LPPL notice on the original code does not settle it. CTAN asks
uploaders to identify the licences of their material; see its
[upload guidance](https://ctan.org/upload/) and
[licence catalogue](https://ctan.org/license).

## Technical provenance and build dependencies

The work is based on studying the source code of FEMM/xfemm and Triangle, then
reimplementing selected formulations and algorithms in Lua/LuaTeX. No native
C/C++ source or binary from those projects is bundled into the runtime.
This describes the technical approach, not a legal determination that a
language translation removes upstream rights. The mesher has the limitations
stated in the manual and does not claim full Triangle quality guarantees.

The manual uses PGF documentation macros supplied by the user's TeX installation.
Its presentation follows tikz-ext; its original prose and examples are written
for luafemm. Neither PGF's nor tikz-ext's source is redistributed here. Build
tools and TeX dependencies retain their own licences and are not in the archives.

The development used Codex/GPT-6-Astra for "vibe coding". Christophe Jorssen is
the package author and maintainer. The manual and README describe the need for
independent validation of AI-assisted numerical code.
