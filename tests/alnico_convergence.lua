-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
-- Local gap refinement and finite-domain sensitivity, same geometry as the PDF.
local f = require("luafemm")
local model = dofile("tests/alnico_model.lua")
print("h,gap_h,domain_width,nodes,triangles,By_left,By_right,newton,residual,min_angle")
for _, c in ipairs({ { 2, 0.5, 100 }, { 1.5, 0.25, 100 }, { 2, 0.5, 140 } }) do
    local m = model(c[1])
    local extra = c[3] - 100
    m.options.xmin = -c[3]
    m.options.xmax = c[3]
    m.options.ymin = -85 - extra
    m.options.ymax = 100 + extra
    f.region(m, "air", { { -40, 30 }, { -25, 30 }, { -25, 31 }, { -40, 31 } }, 0, 0, c[2])
    f.region(m, "air", { { 25, 30 }, { 40, 30 }, { 40, 31 }, { 25, 31 } }, 0, 0, c[2])
    f.solve(m)
    local _, left = f.sample(m, -32.5, 30.5)
    local _, right = f.sample(m, 32.5, 30.5)
    print(
        string.format(
            "%g,%g,%g,%d,%d,%.6f,%.6f,%d,%.3g,%.3g",
            c[1],
            c[2],
            2 * c[3],
            #m.nodes,
            #m.triangles,
            left,
            right,
            m.stats.newton,
            m.stats.residual,
            m.mesh_stats.min_angle
        )
    )
end
