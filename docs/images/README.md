<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# README previews

These images are rendered from the first TikZ picture in the corresponding
example, using its original geometry, material settings and field calculation.
They inherit the example's LPPL 1.3c licence.

| Image | Source |
| --- | --- |
| `tilted-u.png` | [`curved-core.tex`](../../examples/curved-core.tex) |
| `slotted-machine-field.png` | [`tutorial-machine-sine.tex`](../../examples/tutorial-machine-sine.tex) |
| `slotted-rotor.png` | [`machine-nodes.tex`](../../examples/machine-nodes.tex) |

Regenerate from the repository root with:

```sh
python3 scripts/readme-images.py
```

This documentation-only command needs LuaLaTeX and Poppler's `pdftoppm`.
It compiles fresh standalone wrappers under `build/readme/`, then renders
1280-pixel PNGs on a white background for readable previews in GitHub themes.
The original example sources remain the authoritative, editable figures.
