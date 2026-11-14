# Public distribution manifest and LPPL work definition

The original Work is luafemm, copyright (C) 2026 Christophe Jorssen,
licensed under LPPL 1.3c or later, with maintenance status `maintained`.
Christophe Jorssen <christophe.jorssen@gmail.com> is the author and Current
Maintainer. See LICENSES.md for the scope, Base Interpreters and upstream terms.

This is the explicit distribution inventory. The GitHub tree and browsable CTAN
archive contain the same README and preview assets. The compiled
`doc/luatex/luafemm/luafemm-manual.pdf` is added by the documentation build.
Internal tools, tests, assessments, roadmaps, configuration and work history are
excluded. The documentation builder and its inputs are included solely to make
the published manual reproducible. No runtime native solver is distributed.

## Original source and documentation (LPPL 1.3c or later)

- `LICENSE`
- `LICENSES.md`
- `MANIFEST.md`
- `README.md`
- `RELEASE.md`
- `data/femm/README.md`
- `doc/luatex/luafemm/build-manual.py`
- `doc/luatex/luafemm/luafemm-manual.tex`
- `doc/luatex/luafemm/manual-catalog.lua`
- `doc/luatex/luafemm/manual-examples.tex`
- `doc/luatex/luafemm/manual-fragments.py`
- `doc/luatex/luafemm/sections/ac-circuits.tex`
- `doc/luatex/luafemm/sections/axisymmetric.tex`
- `doc/luatex/luafemm/sections/biased-harmonic.tex`
- `doc/luatex/luafemm/sections/boundaries.tex`
- `doc/luatex/luafemm/sections/cache.tex`
- `doc/luatex/luafemm/sections/components.tex`
- `doc/luatex/luafemm/sections/currents.tex`
- `doc/luatex/luafemm/sections/development.tex`
- `doc/luatex/luafemm/sections/domains.tex`
- `doc/luatex/luafemm/sections/field.tex`
- `doc/luatex/luafemm/sections/format-examples.tex`
- `doc/luatex/luafemm/sections/getting-started.tex`
- `doc/luatex/luafemm/sections/harmonic.tex`
- `doc/luatex/luafemm/sections/ideal.tex`
- `doc/luatex/luafemm/sections/interchange.tex`
- `doc/luatex/luafemm/sections/losses.tex`
- `doc/luatex/luafemm/sections/machines.tex`
- `doc/luatex/luafemm/sections/materials.tex`
- `doc/luatex/luafemm/sections/mesh-quality.tex`
- `doc/luatex/luafemm/sections/model.tex`
- `doc/luatex/luafemm/sections/numerics.tex`
- `doc/luatex/luafemm/sections/operating.tex`
- `doc/luatex/luafemm/sections/postprocessing.tex`
- `doc/luatex/luafemm/sections/profiles.tex`
- `doc/luatex/luafemm/sections/project.tex`
- `doc/luatex/luafemm/sections/sliding.tex`
- `doc/luatex/luafemm/sections/static.tex`
- `doc/luatex/luafemm/sections/status.tex`
- `doc/luatex/luafemm/sections/tutorial-machines.tex`
- `doc/luatex/luafemm/sections/tutorial-ndt.tex`
- `doc/luatex/luafemm/sections/tutorial.tex`
- `doc/luatex/luafemm/validation.md`
- `examples/alnico-u.tex`
- `examples/axisymmetric-magnet.tex`
- `examples/axisymmetric-solenoid.tex`
- `examples/biased-harmonic.tex`
- `examples/boundary-conditions.tex`
- `examples/curved-core.tex`
- `examples/eddy-current-bundle.tex`
- `examples/eddy-current-views.tex`
- `examples/forces.tex`
- `examples/formats/context.tex`
- `examples/formats/latex.tex`
- `examples/formats/plain.tex`
- `examples/formats/scene.tex`
- `examples/harmonic-circuits.tex`
- `examples/harmonic-materials.tex`
- `examples/harmonic-shield.tex`
- `examples/harmonic-skin.tex`
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
- `examples/operating-point.tex`
- `examples/plain.tex`
- `examples/sliding-band.tex`
- `examples/static-circuits.tex`
- `examples/static-frontiers.tex`
- `examples/static-open.tex`
- `examples/tutorial-alnico.tex`
- `examples/tutorial-coil.tex`
- `examples/tutorial-machine-positioned.tex`
- `examples/tutorial-machine-sine.tex`
- `examples/tutorial-machine-single.tex`
- `examples/tutorial-machine-uniform.tex`
- `examples/tutorial-ndt-fields.tex`
- `examples/tutorial-ndt-scan.tex`
- `examples/tutorial-ndt-setup.tex`
- `examples/u-electromagnet-node.tex`
- `examples/u-electromagnet.tex`
- `images/slotted-machine-field.png`
- `images/tilted-u.png`
- `tex/luatex/luafemm/lua/luafemm-ac-circuits.lua`
- `tex/luatex/luafemm/lua/luafemm-ans.lua`
- `tex/luatex/luafemm/lua/luafemm-axial.lua`
- `tex/luatex/luafemm/lua/luafemm-axisymmetric.lua`
- `tex/luatex/luafemm/lua/luafemm-bh.lua`
- `tex/luatex/luafemm/lua/luafemm-boundary.lua`
- `tex/luatex/luafemm/lua/luafemm-cache.lua`
- `tex/luatex/luafemm/lua/luafemm-complex.lua`
- `tex/luatex/luafemm/lua/luafemm-components.lua`
- `tex/luatex/luafemm/lua/luafemm-currents.lua`
- `tex/luatex/luafemm/lua/luafemm-fem.lua`
- `tex/luatex/luafemm/lua/luafemm-frontiers.lua`
- `tex/luatex/luafemm/lua/luafemm-geometry.lua`
- `tex/luatex/luafemm/lua/luafemm-harmonic.lua`
- `tex/luatex/luafemm/lua/luafemm-ideal.lua`
- `tex/luatex/luafemm/lua/luafemm-integrals.lua`
- `tex/luatex/luafemm/lua/luafemm-loader.lua`
- `tex/luatex/luafemm/lua/luafemm-losses.lua`
- `tex/luatex/luafemm/lua/luafemm-machines.lua`
- `tex/luatex/luafemm/lua/luafemm-materials.lua`
- `tex/luatex/luafemm/lua/luafemm-mesh.lua`
- `tex/luatex/luafemm/lua/luafemm-operating.lua`
- `tex/luatex/luafemm/lua/luafemm-path.lua`
- `tex/luatex/luafemm/lua/luafemm-post.lua`
- `tex/luatex/luafemm/lua/luafemm-predicates.lua`
- `tex/luatex/luafemm/lua/luafemm-profiles.lua`
- `tex/luatex/luafemm/lua/luafemm-refine.lua`
- `tex/luatex/luafemm/lua/luafemm-shapes.lua`
- `tex/luatex/luafemm/lua/luafemm-sliding.lua`
- `tex/luatex/luafemm/lua/luafemm-static-exchange.lua`
- `tex/luatex/luafemm/lua/luafemm-static.lua`
- `tex/luatex/luafemm/lua/luafemm-tex.lua`
- `tex/luatex/luafemm/lua/luafemm-topology.lua`
- `tex/luatex/luafemm/lua/luafemm.lua`
- `tex/luatex/luafemm/luafemm-boundary.tex`
- `tex/luatex/luafemm/luafemm-components.tex`
- `tex/luatex/luafemm/luafemm-currents.tex`
- `tex/luatex/luafemm/luafemm-ideal.tex`
- `tex/luatex/luafemm/luafemm-interchange.tex`
- `tex/luatex/luafemm/luafemm-machines.tex`
- `tex/luatex/luafemm/luafemm-operating.tex`
- `tex/luatex/luafemm/luafemm-post.tex`
- `tex/luatex/luafemm/luafemm-profiles.tex`
- `tex/luatex/luafemm/luafemm-sliding.tex`
- `tex/luatex/luafemm/luafemm-static.tex`
- `tex/luatex/luafemm/luafemm-tikz.tex`
- `tex/luatex/luafemm/luafemm.sty`
- `tex/luatex/luafemm/luafemm.tex`
- `tex/luatex/luafemm/t-luafemm.mkiv`
- `tex/luatex/luafemm/tikzlibraryfemm.code.tex`
- `tex/luatex/luafemm/tikzlibrarypgfplots.femm.code.tex`

## Separately licensed upstream data and derivatives

These files retain their FEMM/AFPL terms; they are not relicensed as LPPL.

- `data/femm/LICENSE.txt`
- `data/femm/matlib.dat`
- `docs/materials-catalog.md`
- `tex/luatex/luafemm/lua/luafemm-materials-data.lua`
