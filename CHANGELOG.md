<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Changelog

## 0.5.0-rc.1 — 2026-09-25

First release candidate of the generic TikZ and PGFPlots libraries.

- Opening tutorial with two complete U-magnet examples, field lines and
  contour profiles, following the progressive presentation of the PGF tutorials.
  Coarse drafting meshes, local gap refinement and final convergence checks
  are explained separately from graphical detail and profile sampling.

- Original work licensed under LPPL 1.3c, with Christophe Jorssen as author and
  Current Maintainer; upstream FEMM data keep their separate terms.
- Author contact, GitHub development/issue links and an explicit Codex/GPT-6-Astra
  "vibe coding" notice with guidance on independent numerical validation.

- Installable TEXMF layout, with plain LuaTeX, LuaLaTeX and ConTeXt MkIV loaders.
- Explicit `\usepgfplotslibrary{femm}` entry point; PGFPlots stays optional.
- Lua numerical modules separated from their TeX output boundary, documented
  public API, independent input tables, finite numeric-key validation, consistent
  formatting and namespaced private TeX macros.
- English reference manual with indexed keys, defaults, units, executed examples,
  numerical limitations, Lua API and the full 246-record material catalogue.
- English examples and technical notes, automated format comparisons and local
  source/TDS packaging with provenance and installation verification.

### Migration from the prototypes

Install the `tex/` tree instead of copying files from the repository root.
The documented TikZ path interface and `\femmplot` syntax are preserved.
The original polygon commands remain available for existing plain examples;
new documents should use `femm/region` on native TikZ paths.

Private macros are now under `\luafemm@...`; old undocumented implementation
macros are not compatibility interfaces. Numeric literal keys reject non-finite
values or executable Lua text; keys explicitly documented as PGF expressions
continue to accept expressions. Lua model constructors and declarations copy
input tables, so changing a caller's table no longer changes an existing model.

## Prototype history

- **0.4:** FEMM material catalogue and named straight/curved PGFPlots profiles.
- **0.3:** alnico/soft-iron permanent magnet, local mesh density, filtered exact
  orientation/incircle predicates and mesh stress tests.
- **0.2:** native TikZ paths, transformations, curved geometry and constrained
  Delaunay meshing written in Lua.
- **0.1:** nonlinear planar P1 solver, oriented potential contours and a U-core
  electromagnet through an initial generic command interface.
