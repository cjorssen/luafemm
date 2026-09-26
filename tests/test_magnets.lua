-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local function near(a, b, t)
    assert(math.abs(a - b) < t, string.format("%g != %g", a, b))
end
local function test(name, fn)
    io.write(name .. " ... ")
    fn()
    print("OK")
end
test("exact permanent-magnet interface patch", function()
    local m = f.new({
        unit = 1,
        mesher = "delaunay",
        xmin = -1,
        xmax = 1,
        ymin = -1,
        ymax = 1,
        h = 0.2,
        boundary = function(x, y)
            return 0.8 * math.min(y, 0)
        end,
    })
    f.material(m, "magnet", { mur = 3.7, remanence = 0.8 })
    f.region(m, "magnet", { { -1, -1 }, { 1, -1 }, { 1, 0 }, { -1, 0 } })
    f.solve(m)
    for i, p in ipairs(m.nodes) do
        near(m.A[i], 0.8 * math.min(p[2], 0), 1e-9)
    end
    local bx, by, _, hx, hy = f.sample(m, 0.23, -0.37)
    near(bx, 0.8, 1e-8)
    near(by, 0, 1e-8)
    near(hx, 0, 0.002)
    near(hy, 0, 0.002)
end)
test("uniformly magnetized disk versus analytic infinite-domain solution", function()
    local m = f.new({
        unit = 0.001,
        mesher = "delaunay",
        xmin = -15,
        xmax = 15,
        ymin = -15,
        ymax = 15,
        h = 0.8,
    })
    f.material(m, "magnet", { mur = 3, remanence = 1.2 })
    f.region(m, "air", { { -2, -2 }, { 2, -2 }, { 2, 2 }, { -2, 2 } }, 0, 0, 0.08)
    local p = {}
    for i = 0, 95 do
        local a = i * 2 * math.pi / 96
        p[#p + 1] = { math.cos(a), math.sin(a) }
    end
    f.region(m, "magnet", p, 0, 90)
    f.solve(m)
    local bx, by = f.sample(m, 0.07, 0.11)
    near(bx, 0, 0.005)
    near(by, 1.2 / (3 + 1), 0.006)
end)
test("invalid permanent material combinations", function()
    assert(not pcall(f.material, f.new(), "bad", { remanence = 1, coercivity = 1000 }))
    assert(not pcall(f.material, f.new(), "bad", { remanence = 1, bh = { { 0, 0 }, { 1, 100 } } }))
    assert(not pcall(f.material, f.new(), "bad", { coercivity = -1 }))
end)
test("Alnico / soft iron, zero current and magnetization reversal", function()
    local model = require("alnico_model")
    local m = model(3)
    f.solve(m)
    local _, left = f.sample(m, -32.5, 30.5)
    local _, right = f.sample(m, 32.5, 30.5)
    assert(left > 0 and right < 0 and left < 1.5)
    near(m.stats.net_current, 0, 1e-15)
    local reverse = model(3, 0)
    f.solve(reverse)
    for i, a in ipairs(m.A) do
        near(a, -reverse.A[i], 1e-8)
    end
    local nul = model(3)
    nul.materials.alnico.hc = 0
    f.solve(nul)
    near(nul.stats.peak, 0, 1e-15)
    print(string.format("\n  gap By %.6f / %.6f T; Newton %d", left, right, m.stats.newton))
end)
print("All permanent-magnet tests passed.")
