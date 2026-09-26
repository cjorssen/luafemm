<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Permanent magnets

The static constitutive law is **H = ν(|B|) B − Hc m**, with the unit vector m
set by `magnetization angle` in degrees. A linear magnet can use `mu r` with
`remanence` in tesla, which implies `Hc = Br / (μ0 μr)`. A nonlinear magnet uses
`bh` together with `coercivity` in A/m. The B–H values must be the shifted
monotone table expected by FEMM, not the raw negative-H demagnetization curve.

```tex
\tikzset{femm/materials/my-alnico/.style={
  library=cast-alnico-5-lng37}}
\filldraw[femm/region={material=my-alnico,magnetization angle=180},
  fill=red!20] (-25,-30) rectangle (25,-15);
```

Cast Alnico 5 (LNG37) has Hc = 48,600 A/m and its shifted table gives
Br = 1.204 T. `examples/alnico-u.tex` declares that table explicitly, combines
it with a nonlinear illustrative soft iron, and refines two 1 mm gaps locally.
`examples/materials-profiles.tex` instead uses both materials from the catalogue.
The magnetization direction is measured in the model frame: rotating a path's
geometry does not automatically rotate the material's prescribed angle.

This is a prescribed static branch. There is no hysteresis memory, recoil-line
tracking, irreversible demagnetization or magnetic-domain simulation. In
particular, a field solve does not certify that a hard magnet remains on a
physically attainable branch after its field history changes.

Regression tests cover the exact piecewise-affine magnet/air solution, a
uniformly magnetized disk (infinite-domain reference B = Br/(μr+1)), reversal
of magnetization and the zero-coercivity limit. The disk tolerance includes
polygonal approximation, mesh error and finite exterior domain effects.
See [validation](validation.md) for the U magnet's mesh/domain sensitivity.
