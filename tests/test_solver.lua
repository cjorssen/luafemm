-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local function near(a, b, t)
    assert(math.abs(a - b) < t, string.format("%.12g != %.12g (tol %g)", a, b, t))
end
local function test(name, fn)
    io.write(name .. " ... ")
    fn()
    print("OK")
end
test("affine patch, nonzero Dirichlet data", function()
    local m = f.new({
        unit = 1,
        h = 0.2,
        xmin = 0,
        xmax = 1,
        ymin = 0,
        ymax = 1,
        boundary = function(x, y)
            return x + 2 * y
        end,
    })
    f.solve(m)
    for i, p in ipairs(m.nodes) do
        near(m.A[i], p[1] + 2 * p[2], 1e-8)
    end
    local bx, by = f.sample(m, 0.37, 0.41)
    near(bx, 2, 1e-8)
    near(by, -1, 1e-8)
end)
local function manufactured(h)
    local pi = math.pi
    local m = f.new({
        unit = 1,
        h = h,
        xmin = 0,
        xmax = 1,
        ymin = 0,
        ymax = 1,
        source = function(x, y)
            return 2 * pi * pi * math.sin(pi * x) * math.sin(pi * y) / f.mu0
        end,
    })
    f.solve(m)
    local sum = 0
    for i, p in ipairs(m.nodes) do
        sum = sum + (m.A[i] - math.sin(pi * p[1]) * math.sin(pi * p[2])) ^ 2
    end
    return math.sqrt(sum / #m.nodes)
end
test("manufactured Poisson solution, second-order convergence", function()
    local a, b = manufactured(0.1), manufactured(0.05)
    assert(a / b > 3.3 and a / b < 4.8)
    assert(b < 0.003)
end)
test("B-H secant and differential reluctivity", function()
    local p = { bh = { { 0, 0 }, { 1, 100 }, { 2, 1000 } } }
    local nu, d = f.constitutive(p, 1.5)
    near(nu / f.mu0, 550 / 1.5, 1e-9)
    near(d / f.mu0, 900, 1e-9)
    local _, sat = f.constitutive(p, 3)
    near(sat, 1, 1e-12)
    assert(not pcall(f.material, f.new(), "bad", { bh = { { 0, 0 }, { 1, 3 }, { 2, 2 } } }))
end)
test("zero excitation", function()
    local m = f.new({ h = 10 })
    f.solve(m)
    near(m.stats.peak, 0, 1e-15)
    assert(#f.contours(m, 12) == 0)
end)
test("invalid geometry and malformed B-H input are rejected", function()
    assert(not pcall(f.parse_pairs, "0/0,invalid,1/100"))
    local m = f.new()
    f.region(m, "air", { { 0, 0 }, { 2, 0 }, { 1, 1 }, { 0, 1 } })
    assert(not pcall(f.mesh, m), "the structured grid must reject oblique edges")
    assert(
        not pcall(
            f.region,
            f.new(),
            "air",
            { { 0, 0 }, { 3, 0 }, { 3, 3 }, { 1, 3 }, { 1, -1 }, { 0, -1 } }
        )
    )
end)
test("nonlinear U magnet, continuity, contour direction and current reversal", function()
    local model = require("u_model")
    local m = model(3, 85)
    f.solve(m)
    assert(m.stats.residual < 1e-7)
    near(m.stats.net_current, 0, 1e-8)
    local bx, by = f.sample(m, -32.5, 31.5)
    assert(by > 0.5 and by < 2)
    local _, right = f.sample(m, 32.5, 31.5)
    assert(right < -0.5)
    local paths = f.contours(m, 21)
    assert(#paths > 15)
    for _, path in ipairs(paths) do
        near(path[1].x, path[#path].x, 1e-7)
        near(path[1].y, path[#path].y, 1e-7)
        for i = 1, #path - 1 do
            local a, b = path[i], path[i + 1]
            local dx, dy = b.x - a.x, b.y - a.y
            if dx * dx + dy * dy > 1e-14 then
                local x, y, v = f.sample(m, (a.x + b.x) / 2, (a.y + b.y) / 2)
                near(v, path.level, 1e-10)
                assert(dx * x + dy * y >= -1e-10, "reversed field line")
            end
        end
    end
    local reverse = model(3, 85, -6000)
    f.solve(reverse)
    for i, a in ipairs(m.A) do
        near(a, -reverse.A[i], 1e-9)
    end
    local linear = model(3, 85, 6000, true)
    f.solve(linear)
    local _, bl = f.sample(linear, -32.5, 31.5)
    assert(bl > by, "saturation should reduce gap flux")
    print(
        string.format(
            "\n  gap By=%.5f T; linear %.5f T; Newton=%d; residual=%.3g",
            by,
            bl,
            m.stats.newton,
            m.stats.residual
        )
    )
end)
print("All solver tests passed.")
