<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# File inventory and LPPL work definition

The original Work is luafemm, copyright (C) 2026 Christophe Jorssen,
licensed under LPPL 1.3c with maintenance status `maintained`.
Christophe Jorssen <christophe.jorssen@gmail.com> is the author and
Current Maintainer. See `LICENSES.md` for the licence scope, Base Interpreters,
upstream exceptions and the planned CTAN distribution.

## Original source components (LPPL-1.3c)

- `.editorconfig`
- `.gitattributes`
- `.github/workflows/check.yml`
- `.gitignore`
- `.stylua.toml`
- `.styluaignore`
- `CHANGELOG.md`
- `CONTRIBUTING.md`
- `LICENSES.md`
- `MANIFEST.md`
- `Makefile`
- `README.md`
- `data/femm/README.md`
- `doc/luafemm-manual.tex`
- `doc/manual-examples.tex`
- `doc/sections/boundaries.tex`
- `doc/sections/cache.tex`
- `doc/sections/components.tex`
- `doc/sections/development.tex`
- `doc/sections/domains.tex`
- `doc/sections/field.tex`
- `doc/sections/format-examples.tex`
- `doc/sections/getting-started.tex`
- `doc/sections/ideal.tex`
- `doc/sections/interchange.tex`
- `doc/sections/machines.tex`
- `doc/sections/materials.tex`
- `doc/sections/mesh-quality.tex`
- `doc/sections/model.tex`
- `doc/sections/numerics.tex`
- `doc/sections/profiles.tex`
- `doc/sections/project.tex`
- `doc/sections/tutorial-machines.tex`
- `doc/sections/tutorial.tex`
- `docs/femm-interoperability-audit.md`
- `docs/interoperability.md`
- `docs/materials.md`
- `docs/meshing.md`
- `docs/method.md`
- `docs/permanent-magnets.md`
- `docs/profiles.md`
- `docs/validation.md`
- `examples/alnico-u.tex`
- `examples/boundary-conditions.tex`
- `examples/curved-core.tex`
- `examples/formats/context.tex`
- `examples/formats/latex.tex`
- `examples/formats/plain.tex`
- `examples/formats/scene.tex`
- `examples/ideal-e.tex`
- `examples/ideal-paths.tex`
- `examples/ideal-toroids.tex`
- `examples/ideal-u.tex`
- `examples/interchange-alnico.tex`
- `examples/interchange-coil.tex`
- `examples/interchange-curves.tex`
- `examples/machine-nodes.tex`
- `examples/materials-profiles.tex`
- `examples/mesh-quality.tex`
- `examples/nested-domains.tex`
- `examples/plain.tex`
- `examples/tutorial-alnico.tex`
- `examples/tutorial-coil.tex`
- `examples/tutorial-machine-positioned.tex`
- `examples/tutorial-machine-sine.tex`
- `examples/tutorial-machine-single.tex`
- `examples/tutorial-machine-uniform.tex`
- `examples/u-electromagnet-node.tex`
- `examples/u-electromagnet.tex`
- `scripts/build.py`
- `scripts/import-materials.lua`
- `scripts/interop.py`
- `scripts/manual-catalog.lua`
- `scripts/release.py`
- `tests/alnico_convergence.lua`
- `tests/alnico_model.lua`
- `tests/bootstrap.lua`
- `tests/check_predicates.py`
- `tests/components-invalid.tex`
- `tests/convergence.lua`
- `tests/fixtures/uniform.fem`
- `tests/formats/boundaries.tex`
- `tests/formats/components.tex`
- `tests/formats/context.tex`
- `tests/formats/ideal.tex`
- `tests/formats/interchange.tex`
- `tests/formats/latex.tex`
- `tests/formats/machines.tex`
- `tests/formats/plain.tex`
- `tests/formats/scene.tex`
- `tests/interop_cases.lua`
- `tests/interop_machines.lua`
- `tests/machine-model.lua`
- `tests/test_api.lua`
- `tests/test_boundary.lua`
- `tests/test_cache.lua`
- `tests/test_components.lua`
- `tests/test_domains.lua`
- `tests/test_ideal.lua`
- `tests/test_interchange.lua`
- `tests/test_machines.lua`
- `tests/test_magnets.lua`
- `tests/test_materials.lua`
- `tests/test_mesh.lua`
- `tests/test_mesh_stress.lua`
- `tests/test_profiles.lua`
- `tests/test_quality.lua`
- `tests/test_solver.lua`
- `tests/tikz-boundaries.tex`
- `tests/tikz-components.tex`
- `tests/tikz-domains.tex`
- `tests/tikz-ideal.tex`
- `tests/tikz-interchange.tex`
- `tests/tikz-interface.tex`
- `tests/tikz-machines.tex`
- `tests/tikz-profiles.tex`
- `tests/tikz-quality.tex`
- `tests/u_model.lua`
- `tests/xfemm_probe.cpp`
- `tex/context/third/luafemm/t-luafemm.mkiv`
- `tex/generic/luafemm/luafemm-ans.lua`
- `tex/generic/luafemm/luafemm-bh.lua`
- `tex/generic/luafemm/luafemm-boundary.lua`
- `tex/generic/luafemm/luafemm-boundary.tex`
- `tex/generic/luafemm/luafemm-cache.lua`
- `tex/generic/luafemm/luafemm-components.lua`
- `tex/generic/luafemm/luafemm-components.tex`
- `tex/generic/luafemm/luafemm-fem.lua`
- `tex/generic/luafemm/luafemm-geometry.lua`
- `tex/generic/luafemm/luafemm-ideal.lua`
- `tex/generic/luafemm/luafemm-ideal.tex`
- `tex/generic/luafemm/luafemm-integrals.lua`
- `tex/generic/luafemm/luafemm-interchange.tex`
- `tex/generic/luafemm/luafemm-machines.lua`
- `tex/generic/luafemm/luafemm-machines.tex`
- `tex/generic/luafemm/luafemm-materials.lua`
- `tex/generic/luafemm/luafemm-mesh.lua`
- `tex/generic/luafemm/luafemm-path.lua`
- `tex/generic/luafemm/luafemm-predicates.lua`
- `tex/generic/luafemm/luafemm-profiles.lua`
- `tex/generic/luafemm/luafemm-profiles.tex`
- `tex/generic/luafemm/luafemm-refine.lua`
- `tex/generic/luafemm/luafemm-shapes.lua`
- `tex/generic/luafemm/luafemm-tex.lua`
- `tex/generic/luafemm/luafemm-tikz.tex`
- `tex/generic/luafemm/luafemm-topology.lua`
- `tex/generic/luafemm/luafemm.lua`
- `tex/generic/luafemm/luafemm.tex`
- `tex/generic/luafemm/tikzlibraryfemm.code.tex`
- `tex/generic/luafemm/tikzlibrarypgfplots.femm.code.tex`
- `tex/latex/luafemm/luafemm.sty`

## Aggregated upstream notices and data derivatives (separate FEMM terms)

- `data/femm/LICENSE.txt`
- `data/femm/matlib.dat`
- `docs/materials-catalog.md`
- `tex/generic/luafemm/luafemm-materials-data.lua`

## Licence text (verbatim LaTeX3 Project text)

- `LICENSE`

## Compiled and generated components

The manual is compiled as `build/doc/luafemm-manual.pdf` and installed as
`doc/generic/luafemm/luafemm-manual.pdf`. Its original text is LPPL; the
material appendix generated from the FEMM catalogue retains upstream terms.
Example PDFs are compiled from the files under `examples/` and inherit the
corresponding original-work licence, subject to any included third-party data.
`build/doc/materials-catalog.tex` is an intermediate catalogue derivative.
Logs, caches, tool downloads and release ZIP containers are build products,
not additional original source components. Source and TDS archives preserve
the notices applicable to each included component.
