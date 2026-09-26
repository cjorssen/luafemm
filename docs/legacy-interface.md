<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Prototype command compatibility

The original polygon interface remains available in plain LuaTeX, LaTeX and
ConTeXt MkIV. It defaults to the structured `grid` mesher and axis-aligned
physical edges; new documents should use native TikZ paths. The manual appendix
lists every legacy option and initial value.

| Original command | Native-path equivalent |
| --- | --- |
| `\femmbegin[options]` | Picture option `femm/problem={options}` |
| `\femmmaterial[options]{name}` | `femm/materials/name/.style={options}` |
| `\femmregion[options]{x/y,...}` | `\filldraw[femm/region={options}] ... -- cycle;` |
| `\femmsolve` | `\path[femm/solve];` or implicit solving by `femm field` |
| `\femmfieldlines[lines=35]` | `\pic[femm/lines=35] {femm field};` |
| `\femmend` | End the TikZ picture |

`femm region` (a style name with a space) controls legacy polygon appearance.
`femm/region` (a key with a slash) captures a native path as a physical region.
The default contour counts differ: 35 in the legacy command, 29 in the modern
picture. `examples/plain.tex` exercises the retained legacy commands.
