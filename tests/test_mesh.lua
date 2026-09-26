-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local g = require("luafemm-mesh")
local function near(a, b, t)
    assert(math.abs(a - b) < t, string.format("%g != %g", a, b))
end
local function edge(a, b)
    return math.min(a, b) .. ":" .. math.max(a, b)
end
local function audit(m)
    local area, edges = 0, {}
    for _, t in ipairs(m.triangles) do
        assert(t.area > 0)
        area = area + t.area
        for i = 1, 3 do
            local k = edge(t.ids[i], t.ids[i % 3 + 1])
            edges[k] = (edges[k] or 0) + 1
            assert(edges[k] <= 2)
        end
    end
    near(
        area,
        (m.options.xmax - m.options.xmin) * (m.options.ymax - m.options.ymin) * m.options.unit ^ 2,
        1e-10
    )
    for _, s in ipairs(m.mesh_stats.segments) do
        assert(edges[edge(s[1], s[2])], "missing constraint")
    end
end
local function model(h)
    return f.new({
        mesher = "delaunay",
        unit = 1,
        xmin = -10,
        xmax = 10,
        ymin = -10,
        ymax = 10,
        h = h or 1.2,
    })
end
local function area(m, name)
    local s = 0
    for _, e in ipairs(m.triangles) do
        if e.material == name then
            s = s + e.area
        end
    end
    return s
end
local function test(name, fn)
    io.write(name .. " ... ")
    fn()
    print("OK")
end
test("oblique region, exact area and all constrained segments", function()
    local m = model()
    f.material(m, "iron", { mur = 100 })
    local poly = { { -5, -3 }, { 4, -4 }, { 3, 5 }, { -4, 4 } }
    f.region(m, "iron", poly)
    f.mesh(m)
    audit(m)
    local expected = 0
    for i, p in ipairs(poly) do
        local q = poly[i % #poly + 1]
        expected = expected + p[1] * q[2] - p[2] * q[1]
    end
    near(area(m, "iron"), expected / 2, 1e-10)
end)
test("crossings, overlap, shared edges and T junctions", function()
    local m = model()
    f.material(m, "a", { mur = 10 })
    f.material(m, "b", { mur = 20 })
    f.region(m, "a", { { -6, -4 }, { 2, -4 }, { 2, 4 }, { -6, 4 } })
    f.region(m, "b", { { -2, -2 }, { 6, -2 }, { 6, 6 }, { -2, 6 } })
    f.region(m, "air", { { 2, -2 }, { 4, -2 }, { 4, 2 }, { 2, 2 } })
    f.mesh(m)
    audit(m)
    near(area(m, "a"), 40, 1e-9)
    near(area(m, "b"), 56, 1e-9)
end)
test("circle approximation and nested air hole", function()
    local m = model()
    f.material(m, "iron", { mur = 100 })
    local function circle(r)
        local p = {}
        for i = 0, 63 do
            local a = 2 * math.pi * i / 64
            p[#p + 1] = { r * math.cos(a), r * math.sin(a) }
        end
        return p
    end
    f.region(m, "iron", circle(6))
    f.region(m, "air", circle(3))
    f.mesh(m)
    audit(m)
    near(area(m, "iron"), 32 * math.sin(2 * math.pi / 64) * (36 - 9), 1e-9)
end)
test("thin gap and concave interface", function()
    local m = model(0.8)
    f.material(m, "iron", { mur = 150 })
    f.region(
        m,
        "iron",
        { { -6, -6 }, { 6, -6 }, { 6, 5 }, { 4, 5 }, { 4, -4 }, { -4, -4 }, { -4, 5 }, { -6, 5 } }
    )
    f.region(m, "iron", { { -6, 5.1 }, { 6, 5.1 }, { 6, 7 }, { -6, 7 } })
    f.mesh(m)
    audit(m)
end)
test("affine solution on unstructured triangles", function()
    local m = model()
    m.options.boundary = function(x, y)
        return 2 * x - 3 * y
    end
    f.solve(m)
    audit(m)
    for i, p in ipairs(m.nodes) do
        near(m.A[i], 2 * p[1] - 3 * p[2], 2e-7)
    end
    local bx, by = f.sample(m, 0.13, 0.71)
    near(bx, -3, 1e-7)
    near(by, -2, 1e-7)
end)
test("manufactured solution converges under refinement", function()
    local function error(h)
        local pi = math.pi
        local m = f.new({
            mesher = "delaunay",
            unit = 1,
            xmin = 0,
            xmax = 1,
            ymin = 0,
            ymax = 1,
            h = h,
            source = function(x, y)
                return 2 * pi * pi * math.sin(pi * x) * math.sin(pi * y) / f.mu0
            end,
        })
        f.solve(m)
        audit(m)
        local sum = 0
        for i, p in ipairs(m.nodes) do
            sum = sum + (m.A[i] - math.sin(pi * p[1]) * math.sin(pi * p[2])) ^ 2
        end
        return math.sqrt(sum / #m.nodes)
    end
    local a, b = error(0.12), error(0.06)
    assert(a / b > 3 and b < 0.004)
end)
test("nonlinear electromagnet on constrained Delaunay mesh", function()
    -- h=2 exercises an IC(0) pivot breakdown and the shifted retry.
    local m = require("u_model")(2, 110)
    m.options.mesher = "delaunay"
    f.solve(m)
    audit(m)
    assert(m.stats.residual < 1e-7)
    near(m.stats.net_current, 0, 1e-7)
    local _, left = f.sample(m, -32.5, 31.5)
    local _, right = f.sample(m, 32.5, 31.5)
    assert(left > 1 and left < 1.2 and right < -0.85 and right > -1.05)
    print(
        string.format(
            "\n  %d nodes, %d triangles; gap By %.6f / %.6f T",
            #m.nodes,
            #m.triangles,
            left,
            right
        )
    )
end)
print("All constrained-mesh tests passed.")
