-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Fixtures for the optional independent xfemm reader/solver.
dofile("tests/bootstrap.lua")
require("lfs").mkdir("build/static")
local M, F, S = require("luafemm"), require("luafemm-frontiers"), require("luafemm-static")
local rows = {}
local function model()
    return M.new({
        xmin = -12,
        xmax = 12,
        ymin = -12,
        ymax = 12,
        h = 0.4,
        depth = 10,
        mesher = "delaunay",
        tolerance = 1e-10,
    })
end
local function save(m, name, points)
    M.solve(m)
    M.export_fem(m, "build/static/" .. name .. ".fem")
    M.export_ans(m, "build/static/" .. name .. "-lua.ans")
    for _, p in ipairs(points) do
        local bx, by, a, hx, hy = M.sample(m, p[1], p[2])
        local topology = require("luafemm-topology")
        local face = topology.locate(m.topology or topology.native(m), p[1], p[2])
        rows[#rows + 1] = string.format(
            "%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g",
            name,
            p[1],
            p[2],
            a,
            bx,
            by,
            hx,
            hy,
            face.current
        )
    end
end
do
    local m = model()
    M.material(m, "tensor", { mur_x = 3, mur_y = 9 })
    M.region(m, "tensor", { { -12, -12 }, { 12, -12 }, { 12, 12 }, { -12, 12 } }, 0, 0, 0.4)
    M.boundary(m, "all", { dx = -0.3, dy = 0.2 })
    save(m, "tensor", { { 1.3, 2.7 }, { -5.1, 3.4 } })
end
do
    local m = model()
    M.boundary(m, "all", { type = "neumann" })
    F.path(m, "left", { { -12, -12 }, { -12, 12 } }, false, { type = "periodic" })
    F.path(m, "right", { { 12, -12 }, { 12, 12 } }, false, { type = "periodic" })
    F.pair(m, "left", "right", { x = 24, sign = -1 })
    M.region(m, "air", { { -12, -12 }, { 12, -12 }, { 12, 12 }, { -12, 12 } }, 0, 0, 0.4)
    M.material(m, "wire", { conductivity = 5.8e7 })
    M.region(m, "wire", { { -4, -5 }, { 4, -5 }, { 4, 5 }, { -4, 5 } }, 2e6, 0, 0.2)
    save(m, "antiperiodic", { { 1.3, 2.7 }, { -7.1, 3.4 } })
end
do
    local m = model()
    M.region(m, "air", { { -12, -12 }, { 12, -12 }, { 12, 12 }, { -12, 12 } }, 0, 0, 0.4)
    M.material(m, "wire", { conductivity = 5.8e7 })
    S.circuit(m, "coil", { current = 2 })
    for k, x in ipairs({ -7, 4 }) do
        local r =
            M.region(m, "wire", { { x, -5 }, { x + 3, -5 }, { x + 3, 5 }, { x, 5 } }, 0, 0, 0.2)
        S.region(r, { circuit = "coil", turns = k == 1 and 50 or -50, fill = 0.7 })
    end
    save(m, "circuit", { { 1.3, 2.7 }, { -5.1, 3.4 } })
end
local function circle(r, n)
    local p = {}
    for i = 1, n do
        local a = 2 * math.pi * (i - 1) / n
        p[i] = { r * math.cos(a), r * math.sin(a) }
    end
    return p
end
do
    local m = model()
    local outer = circle(11, 64)
    M.region(m, "air", outer, 0, 0, 0.4)
    F.domain(m, { outer })
    F.path(m, "open", outer, true, { type = "open", radius = 11 })
    M.material(m, "pm", { mur = 1.05, coercivity = 8e5 })
    local r = M.region_contours(m, "pm", { circle(7, 64), circle(4, 64) }, 0, 0, 0.2, "even odd")
    S.region(r, { magnetization = "angular", order = 2 })
    save(m, "angular", { { 1.3, 2.7 }, { -5.1, 3.4 } })
end
local file = assert(io.open("build/static/reference.tsv", "w"))
file:write(table.concat(rows, "\n"), "\n")
file:close()
print("Static interchange fixtures generated")
