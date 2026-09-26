<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# FEMM material library

The package ships all **246 records** of FEMM's magnetic `matlib.dat`: permanent
magnets, soft magnetic alloys, steels, nonmagnetic metals and conductor records.
This is the magnetic solver's catalogue; separate thermal and electrostatic
libraries are outside the scope. [The complete index](materials-catalog.md)
lists exact identifiers, original names, folders, B–H counts and source currents.

```tex
\filldraw[femm/region={material=pure-iron},fill=gray!20] ... -- cycle;
\filldraw[femm/region={material=cast-alnico-5-lng37,
  magnetization angle=180},fill=red!20] ... -- cycle;
\path[femm/region={material=copper,current density=0}] ... -- cycle;
```

`/femm/region/material` is a PGF `.is choice` selector; `femm/material` forwards
to it. ASCII identifiers avoid punctuation in original names. Two records called
Supermalloy have distinct identifiers and remain separate. Unknown choices only
fall back to an already-defined custom style under `/tikz/femm/materials`.

A derived material can inherit a catalogue record:

```tex
\tikzset{femm/materials/my-magnet/.style={
  library=cast-alnico-5-lng37,coercivity=48000}}
```

Custom property keys are `library`, `mu r`, `bh`, `remanence` and `coercivity`.
The manual documents combinations and defaults. Keep custom names distinct from
catalogue identifiers, which select their catalogue records directly.

All original numeric fields and B–H points are retained in the data module.
The planar static solver interprets the supported subset: nonlinear curves,
coercivity, permeability, in-plane lamination fill factors and material source
current. Ten nonlinear curves have nonunit fill factors; 69 records have a
nonzero source current. Region `current density=0` explicitly overrides an
inherited material current. Conductivity, frequency losses, wire geometry and
other retained metadata do not imply eddy-current, harmonic or winding models.
Unsupported anisotropy/lamination combinations are rejected. See the manual for
the precise static transformations.

Lua access uses `require('luafemm-materials').get(id)`, `.list()` and
`.properties(id)`. Returned records are independent copies. `make materials`
regenerates the module and this index from preserved upstream data without FEMM.
The [upstream notices](../data/femm/README.md) remain separate from the package's
original LPPL code.
