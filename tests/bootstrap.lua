-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Development scripts run from the repository root. Installed documents use
-- the engine's normal TEXMF resolver instead of changing package.path.
package.path = "./tex/generic/luafemm/?.lua;./tests/?.lua;" .. package.path
