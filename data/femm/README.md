<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# FEMM data provenance

`matlib.dat` is an unchanged copy of `bin/matlib.dat` from
[cenit/FEMM](https://github.com/cenit/FEMM/blob/7d9e8ed90772ddb939a3256f976e231bcb5cca97/bin/matlib.dat),
revision `7d9e8ed90772ddb939a3256f976e231bcb5cca97`, obtained from the local
`../femm` checkout. FEMM copyright: David C. Meeker.

`LICENSE.txt` is an unchanged copy of that revision's `license.txt`. It contains
FEMM's distribution terms (Aladdin Free Public License) and upstream dependency
notices. Those dependencies' algorithms and binaries are not bundled in luafemm.

`../../tex/generic/luafemm/luafemm-materials-data.lua` is a transcription of this
catalogue made on 2026-09-25: identifiers are added, and folders, names, all
numeric fields and B–H points are retained. `../../docs/materials-catalog.md`
indexes it. These data remain separate from luafemm's original LPPL code and
retain their upstream terms; they are not presented as LPPL data.

Run `make materials` from the repository root to regenerate the derivatives.
The importer is written in Lua and needs neither the network nor FEMM.
