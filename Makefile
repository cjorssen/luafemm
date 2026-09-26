# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
PYTHON ?= python3
STYLUA ?= stylua
LUA_SOURCES = tex/generic/luafemm scripts tests

.PHONY: all test formats examples manual materials convergence alnico-convergence format check-format dist
all:
	$(PYTHON) scripts/build.py all
test formats examples manual:
	$(PYTHON) scripts/build.py $@
materials:
	texlua scripts/import-materials.lua data/femm/matlib.dat
convergence:
	texlua tests/convergence.lua
alnico-convergence:
	texlua tests/alnico_convergence.lua
format:
	$(STYLUA) $(LUA_SOURCES)
check-format:
	$(STYLUA) --check $(LUA_SOURCES)
dist:
	$(PYTHON) scripts/release.py
