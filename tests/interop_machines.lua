-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local make = require("machine-model")
local file = assert(io.open("build/interop/machine-samples.tsv", "w"))
for _, case in ipairs({
    {
        name = "machine-coil",
        slots = { count = 2, start_angle = 90 },
        winding = {},
        extra = { coils = { { slots = { 1, 2 }, ampere_turns = 200 } } },
    },
    {
        name = "machine-sine",
        slots = { count = 12, start_angle = 15 },
        winding = { distribution = "sine", ampere_turns = 200 },
    },
}) do
    local m = make(case.slots, case.winding, 2, false, case.extra)
    f.solve(m)
    local base = "build/interop/" .. case.name
    f.export_fem(m, base .. ".fem")
    f.export_ans(m, base .. "-lua.ans")
    for i = 0, 11 do
        local a = math.rad(i * 30 + 7)
        local x, y = 29 * math.cos(a), 29 * math.sin(a)
        local bx, by, az, hx, hy = f.sample(m, x, y)
        file:write(
            string.format(
                "%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
                case.name,
                x,
                y,
                az,
                bx,
                by,
                hx,
                hy
            )
        )
    end
end
file:close()
