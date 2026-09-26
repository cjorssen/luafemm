<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Contributing

The author and Current Maintainer is Christophe Jorssen
<christophe.jorssen@gmail.com>; the LPPL maintenance status is `maintained`.
Original contributions use LPPL 1.3c. See `LICENSES.md` for the scope and
third-party exceptions. Development and bug reports use
[GitHub](https://github.com/cjorssen/luafemm) and its
[issue tracker](https://github.com/cjorssen/luafemm/issues).

All source, comments, examples, commit messages and documentation are written
in English. Keep the installed package generic: no LaTeX-only command in the
`tex/generic` tree, no shell command or native-library dependency at runtime.

## Source conventions

Lua targets Lua 5.3 and returns local module tables. Only the generic TeX
loader creates the documented `luafemm` global. Keep TeX output in
`luafemm-tex.lua`; numerical modules must also work under `texlua`.
Do not change `package.path` in installed modules. Public entry points have
LDoc-style comments describing units, ownership, arguments and results;
private helpers explain numerical invariants and non-obvious decisions.
Model declarations copy caller-owned tables. Report invalid input rather
than silently substituting a physically different result.

Use StyLua 2.5.2 with `.stylua.toml`: four spaces and a 100-column target.
Generated `luafemm-materials-data.lua` is excluded; change the importer instead.
With StyLua installed, run `make format` or `make check-format`. Set
`STYLUA=build/tools/stylua` if using a local downloaded binary.

TeX private macros use `\luafemm@...`; public keys belong to the documented
PGF families. Preserve incoming catcodes, use grouping for temporary state,
and use package-owned registers instead of global scratch assignments.
Pass user strings through `\luaescapestring` and validate numbers in Lua;
do not splice user input into executable Lua expressions. Keep wrappers thin.
Document each new key's default, units, inheritance/reset behaviour and errors
in the PGF manual environments, with an executable example where useful.

## Validation

Run these commands from the repository root:

```sh
make test
make examples
make manual
make check-format
make dist
```

`make test` runs analytic and manufactured solutions, nonlinear and permanent
magnet regressions, mesh stress cases, an independent exact-predicate oracle,
catalogue and profile tests, and TeX integration. The same magnetic scene
and curved profile are compiled under all three supported formats; the
numerical snapshots must agree. Results are written to `build/`.

The MkIV test invokes TeX Live's `mtxrun.lua` and `mtx-context.lua` explicitly
under `texlua` and selects `--luatex`; it must not silently test LMTX instead.
ConTeXt uses an absolute cache directory because format creation changes its
working directory. LaTeX font caches use a relative path for TeX's file-opening
rules. No global TeX configuration is changed.

`make convergence` and `make alnico-convergence` run longer sensitivity studies.
A small algebraic residual is not a discretization-error estimate. For physical
changes, compare against an analytic result or an independent formulation,
and check both mesh spacing and the artificial boundary location.

`make materials` regenerates the Lua catalogue and Markdown index from the
preserved upstream file. Do not edit the upstream bytes or generated files by
hand. `make dist` checks that regeneration leaves both outputs identical.

## Local release checklist

1. Update the version/date in the generic loader, Lua module and LaTeX wrapper;
   update the manual, changelog and README together.
2. Run all checks above. Inspect the manual and representative example PDFs,
   including arrow directions, profile units, labels and the index.
3. Read `docs/validation.md` and keep evidence separate from physical guarantees.
4. Review `LICENSES.md` and the file scope in `MANIFEST.md`; retain the FEMM
   data's upstream notices in every archive. Before the planned CTAN submission,
   resolve or explicitly declare those separate data terms. Do not describe
   the current mixed-licence archives as an entirely LPPL distribution.
5. Inspect the local source/TDS archives and their SHA-256 sums. The distribution
   check extracts the TDS archive, then recompiles all three format tests using
   that installation from outside the checkout's runtime tree.
6. Review the files to be staged. Build products, editor files and private agent
   configuration do not belong in the public repository.

Archive creation never commits, tags, uploads or publishes. Those are separate
maintainer actions after review. CI repeats the portable build checks on Linux;
a local success does not imply that a remote workflow has already run.
