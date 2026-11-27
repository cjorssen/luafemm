# FEMM 5 material data provenance

The default catalogue is an unchanged copy of `matlibs/matlib.dat` from
[dcm3c/femm-qt](https://github.com/dcm3c/femm-qt/tree/31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3),
revision `31ee4bf6a8cae1a1860c108e7934ab3403bfd6e3` (282 magnetic records).
`UPSTREAM.md` is that revision's unchanged `matlibs/README.md`.
`LICENSE.txt` is its unchanged `LICENSE-MIT`, copyright 2026 David Meeker and
contributors. References within the unchanged upstream explanation are relative
to the upstream repository; the local permission notice is `LICENSE.txt`.
No provenance papers, solver source or native binaries are copied.

The collection and its generated transcription/index retain MIT terms. The
original Lua reader and package code are LPPL 1.3c or later; see `../../LICENSES.md`.
The transcription retains names, folders, numeric fields and B--H samples, adding
ASCII identifiers. Different editions need not have identical names or curves.
`armco-pure-iron` is the default soft-iron component material. The former
`pure-iron` identifier is not an alias: it belongs to another data edition.
FEMM 5 CCA records have zero default imposed current; their conductor properties
are available for DC and the documented effective harmonic CCA models.
Their fit domain and homogenisation assumptions remain explicit limitations.

The optional FEMM 4 edition requires an exact, user-supplied legacy `matlib.dat`
from cenit/FEMM revision `7d9e8ed90772ddb939a3256f976e231bcb5cca97`.
It is not included in the distribution. Choose the edition and file in
`femm/problem` before resolving material styles. FEM/ANS files carry their own
material values independently of the catalogue choice. Catalogue provenance is
part of native solution-cache identity; a cache from another edition is not reused.

## Byte checksums (SHA-256)

- `matlib.dat`: `aa637415d661cce0e5da8d8947926a5a248e05abc47fcf182963b5e90b3030f3`
- `LICENSE.txt`: `7282de29815adde5d613ca6fcef01ea640e73bfffb82d332f853d65d6c4fe97e`
- `UPSTREAM.md`: `5efd36dc80cfba19b7c5c81d21cb035562f1ed785595bf6ebb801ce0bd26c7a6`
- Optional legacy FEMM 4 `matlib.dat`: `ba0ca80ffae56440a9ad785c9089b0e7f63b45e1c7a73dd5c0bdc338beeb63b7`

The runtime transcription preserves the pinned records and their upstream
notices. Catalogue lookups require no network access at document compilation.

## Electrostatic catalogue

`statlib.dat` is an unchanged copy of `matlibs/statlib.dat` at the same pinned
revision, with 95 dielectric records (94 isotropic, one anisotropic). The
MIT notice above covers this library explicitly. SHA-256: `95d733c22e21477d0d7a2395c86a8e4dac51234c5b20d157287a32785e178694`.
The finite-element electric provider supports the anisotropic silicon-carbide
entry. Analytical homogeneous scalar providers reject it instead of averaging it.
All supplied volume charge densities are zero. No FEMM 4 dielectric data are bundled.
