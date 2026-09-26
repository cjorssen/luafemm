<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Field profiles with PGFPlots

Attach a named measurement path to ordinary TikZ drawing:

```tex
\draw[femm/profile={name=gap,samples=301},dashed]
  (-48,30.5)--(48,30.5);
\draw[femm/profile={name=arc,samples=201,curve tolerance=.01}]
  (-18,-8) .. controls (-22,22) and (22,22) .. (18,-8);
```

Then, after loading PGFPlots and `\usepgfplotslibrary{femm}`, use an axis:

```tex
\femmplot[profile=gap,component=By,abscissa=x,
  plot options={blue,thick}]
\addlegendentry{$B_y$}
```

`\femmplot` behaves like `\addplot+`: cycle lists and legends remain available.
The manual gives complete plain, LaTeX and ConTeXt loading examples.

| Capture key | Initial value | Meaning |
| --- | --- | --- |
| `name` | empty; required | Document-wide identifier |
| `samples` | 201 | Integer from 2 to 100000, endpoints included |
| `curve tolerance` | .03 | Flattening tolerance in model units |
| `offset` | 0 | Signed displacement along the left normal, in model units |
| `outside` | error | `error` or `nan` for points outside the rectangular domain |

| Plot key | Reset value | Meaning |
| --- | --- | --- |
| `profile` | empty; required | Registered identifier |
| `component` | B | `B,Bx,By,Bt,Bn,H,Hx,Hy,Ht,Hn,Az` |
| `abscissa` | s | `s,t,x,y` |
| `plot options` | empty | Ordinary PGFPlots options |

B components are in tesla, H components in A/m, Az in T m. Position and arc
length use model units (mm in TikZ). `t=s/L` is normalized arc length, not a
Bezier polynomial parameter. The unit normal is left of the traversal direction;
reversing traversal reverses tangential/normal signs at the same point. At an
internal corner, the next segment's tangent is used; the final endpoint uses
the preceding segment. A closed path repeats its first point at its last sample.
With an offset, s remains the original path length, while x/y are displaced.

Sampling is uniform along the flattened polyline. B and H are constant per P1
triangle; Az is linearly interpolated. No field smoothing is applied. At an
interface the first located triangle supplies the value; use a small signed
offset to select the desired side. More samples improve the plotted sampling,
not the underlying solution. A narrow feature can fall between samples.

A profile keeps its associated model after another picture becomes current.
The first sample request solves that model if needed and caches the results.
Reusing a name replaces the profile. Lua clients who modify a solution manually
must register the profile again to invalidate the cache. `outside=nan` breaks
plots outside the rectangle; it does not hide solver failures.
