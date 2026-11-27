# luafemm 0.36.0-rc.1

This release candidate combines magnetic and electrostatic field diagrams with
one generic TikZ/PGF/PGFPlots interface, for plain LuaTeX, LuaLaTeX and ConTeXt MkIV.
The calculations run in Lua without a native solver or shell escape.

- Planar and axisymmetric magnetic finite elements: nonlinear static materials,
  permanent magnets, supported harmonic and transient fields, winding ports,
  electrical networks, forces, machine studies and bounded FEMM interchange.
- Analytical electrostatics, finite-element dielectrics and electrodes, floating
  conductors, capacitance matrices and singular-source correction. Oriented
  field lines, equipotentials and profiles use the same declared model.
- Prescribed 3D wires and coils use Biot--Savart integration. Owned magnetic and
  electric scenes share cameras, axes, reference planes and vector visibility.
  Curved PGF arrows follow the displayed paths in both 2D and 3D.
- The MIT-licensed FEMM 5 catalogues supply 282 magnetic and 95 dielectric records.
- An extensive illustrated English manual executes its displayed PGF examples.
  The public tree uses shallow `tex/`, `doc/` and `examples/` directories; a
  separate TDS archive supplies installation paths.

Scientific scope, unsupported combinations, conventions and measured validation
are described in the manual and [validation notes](doc/validation.md). This is
not full FEMM parity or a general 3D finite-element solver. A catalogue entry
does not imply support for every constitutive law. Refine geometry, mesh, domain
and field tracing independently for quantitative work. Native cache revision 25
requires recomputing incompatible cached solutions.

Download the browsing ZIP for documentation and sources, or extract the TDS ZIP
at a TEXMF root for installation. Keep all runtime files from this version
together. The PDF manual and SHA-256 checksums are separate release assets.
CTAN submission has been prepared but has not been made.

Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>.
Maintenance status: maintained. Original work: LPPL 1.3c or later. Bundled FEMM 5
catalogue data retain their MIT terms; see [LICENSES.md](LICENSES.md).

The package was vibe-coded with OpenAI Codex and the Astra-6 model, under
Christophe Jorssen's direction. Neither AI assistance nor publication certifies
scientific correctness; independent validation remains necessary.
