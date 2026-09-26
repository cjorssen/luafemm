-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local femm = require("luafemm")
return function(h, extent, ampere_turns, linear)
    extent = extent or 110
    ampere_turns = ampere_turns or 6000
    local m = femm.new({
        h = h or 2,
        xmin = -extent,
        xmax = extent,
        ymin = 15 - extent,
        ymax = 5 + extent,
    })
    femm.material(m, "iron", linear and { mur = 2000 } or {
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
    femm.region(m, "iron", {
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
    femm.region(m, "iron", { { -40, 33 }, { 40, 33 }, { 40, 48 }, { -40, 48 } })
    local j = ampere_turns / (6 * 26 * 1e-6)
    femm.region(m, "air", { { -48, -10 }, { -42, -10 }, { -42, 16 }, { -48, 16 } }, j)
    femm.region(m, "air", { { -23, -10 }, { -17, -10 }, { -17, 16 }, { -23, 16 } }, -j)
    return m
end
