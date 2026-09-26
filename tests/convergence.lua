-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local model = require("u_model")
local mesher = arg[1] or "delaunay"
print("# mesher=" .. mesher)
print("h_mm,half_domain_mm,nodes,triangles,By_left_T,By_right_T,newton,residual,seconds")
for _, case in ipairs({ { 3, 110 }, { 2, 110 }, { 1.5, 110 }, { 2, 160 } }) do
    local m = model(case[1], case[2])
    m.options.mesher = mesher
    f.solve(m)
    local _, left = f.sample(m, -32.5, 31.5)
    local _, right = f.sample(m, 32.5, 31.5)
    print(
        string.format(
            "%g,%g,%d,%d,%.8f,%.8f,%d,%.3g,%.2f",
            case[1],
            case[2],
            m.stats.nodes,
            m.stats.triangles,
            left,
            right,
            m.stats.newton,
            m.stats.residual,
            m.stats.seconds
        )
    )
end
