-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Optional external-validation fixtures, generated from original project models.
dofile("tests/bootstrap.lua")
local M = require("luafemm")
local make_coil = dofile("tests/u_model.lua")
local make_magnet = dofile("tests/alnico_model.lua")
local manifest = assert(io.open("build/interop/samples.tsv", "w"))
for _, kind in ipairs({ "coil", "alnico" }) do
    for _, h in ipairs({ 6, 3 }) do
        local m = kind == "coil" and make_coil(h, 75, 3000) or make_magnet(h, 180, 3)
        m.options.depth = 20
        m.options.interpolation = "femm"
        m.options.mesher = "delaunay"
        for name, p in pairs(m.materials) do
            p.interpolation = "femm"
            M.material(m, name, p)
        end
        M.solve(m)
        local base = "build/interop/" .. kind .. "-" .. h
        M.export_fem(m, base .. ".fem")
        M.export_ans(m, base .. "-lua.ans")
        local count = {}
        for _, e in ipairs(m.triangles) do
            if (count[e.material] or 0) < 12 and e.material ~= "air" then
                local p, q, r = m.nodes[e.ids[1]], m.nodes[e.ids[2]], m.nodes[e.ids[3]]
                local x, y =
                    (p[1] + q[1] + r[1]) / (3 * m.options.unit),
                    (p[2] + q[2] + r[2]) / (3 * m.options.unit)
                local bx, by, a, hx, hy = M.sample(m, x, y)
                manifest:write(
                    string.format(
                        "%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
                        base .. "-lua.ans",
                        x,
                        y,
                        a,
                        bx,
                        by,
                        hx,
                        hy
                    )
                )
                count[e.material] = (count[e.material] or 0) + 1
            end
        end
        local bx, by, a = M.sample(m, -32.5, 31.5)
        manifest:write(
            string.format(
                "%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
                base .. ".ans",
                -32.5,
                31.5,
                a,
                bx,
                by
            )
        )
        print(kind, h, #m.triangles, bx, by)
    end
end
manifest:close()
