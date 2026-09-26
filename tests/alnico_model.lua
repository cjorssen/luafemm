-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
return function(h, angle, gap)
    local m =
        f.new({ mesher = "delaunay", h = h or 2, xmin = -100, xmax = 100, ymin = -85, ymax = 100 })
    f.material(m, "iron", {
        bh = {
            { 0, 0 },
            { 0.5, 80 },
            { 1, 180 },
            { 1.3, 350 },
            { 1.5, 900 },
            { 1.65, 2500 },
            {
                1.8,
                7000,
            },
            {
                2,
                20000,
            },
        },
    })
    -- FEMM material database: Cast Alnico 5 (LNG37), shifted H(B).
    f.material(m, "alnico", {
        coercivity = 48600,
        bh = {
            { 0, 0 },
            { 0.2781, 1600 },
            { 0.6383, 4100 },
            { 0.8234, 6300 },
            { 0.9214, 8600 },
            { 1.013, 13400 },
            { 1.064, 18200 },
            { 1.106, 24300 },
            { 1.149, 33000 },
            { 1.184, 42530 },
            { 1.204, 48600 },
        },
    })
    f.region(m, "iron", {
        { -40, -30 },
        { 40, -30 },
        { 40, 30 },
        { 25, 30 },
        { 25, -15 },
        { -25, -15 },
        { -25, 30 },
        {
            -40,
            30,
        },
    })
    f.region(m, "alnico", { { -25, -30 }, { 25, -30 }, { 25, -15 }, { -25, -15 } }, 0, angle or 180)
    gap = gap or 1
    f.region(
        m,
        "iron",
        { { -40, 30 + gap }, { 40, 30 + gap }, { 40, 45 + gap }, { -40, 45 + gap } }
    )
    return m
end
