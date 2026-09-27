-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Independent analytic patches exercise orientations, SI scaling and weak loads.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local function near(a, b, tolerance)
    assert(math.abs(a - b) < tolerance, string.format("%.15g != %.15g", a, b))
end
local function fails(fn, message)
    local ok, why = pcall(fn)
    assert(not ok and tostring(why):find(message, 1, true), tostring(why))
end
local sides = { "left", "right", "bottom", "top" }
local function model(mesher, unit, cache)
    return f.new({
        xmin = 0,
        xmax = 20,
        ymin = 0,
        ymax = 12,
        h = 2.5,
        unit = unit,
        mesher = mesher,
        tolerance = 1e-10,
        cache = cache,
    })
end
local function patch(mesher, unit, kind, nonlinear)
    local m = model(mesher, unit)
    local a, b, c = kind == "neumann" and 0 or 0.004, 0.3, -0.7
    local p = { mur = 7, coercivity = 23000 }
    if nonlinear then
        p.bh = { { 0, 0 }, { 0.4, 50 }, { 1, 15000 }, { 2, 200000 } }
    end
    f.material(m, "medium", p)
    local angle = 31
    f.region(m, "medium", { { 0, 0 }, { 20, 0 }, { 20, 12 }, { 0, 12 } }, 0, angle)
    p = m.materials.medium
    local nu = f.constitutive(p, math.sqrt(b * b + c * c)) / f.mu0
    local hx = nu * c - p.hc * math.cos(math.rad(angle))
    local hy = -nu * b - p.hc * math.sin(math.rad(angle))
    local ht = { left = -hy, right = hy, bottom = hx, top = -hx }
    for _, side in ipairs(sides) do
        local this = kind
        if kind == "mixed" then
            this = side == "left" and "dirichlet" or side == "top" and "robin" or "neumann"
        end
        if this == "dirichlet" then
            f.boundary(m, side, { value = a, dx = b, dy = c })
        else
            local coefficient = this == "robin" and 30000 or 0
            f.boundary(m, side, {
                type = this,
                value = ht[side] - coefficient * a,
                dx = -coefficient * b,
                dy = -coefficient * c,
                coefficient = coefficient,
            })
        end
    end
    f.solve(m)
    assert((m.gauge_node ~= nil) == (kind == "neumann"))
    for id, point in ipairs(m.nodes) do
        near(m.A[id], a + b * point[1] + c * point[2], 2e-8 * unit)
    end
    for _, triangle in ipairs(m.triangles) do
        near(triangle.bx, c, 1e-7)
        near(triangle.by, -b, 1e-7)
        near(triangle.hx, hx, 0.02)
        near(triangle.hy, hy, 0.02)
    end
end
for _, mesher in ipairs({ "grid", "delaunay" }) do
    for _, unit in ipairs({ 1, 0.001 }) do
        for _, kind in ipairs({ "dirichlet", "neumann", "robin", "mixed" }) do
            patch(mesher, unit, kind, false)
        end
    end
    patch(mesher, 0.001, "mixed", true)
    patch(mesher, 0.001, "neumann", true)
end
print(
    "Boundary affine patches: both meshers, SI units, four orientations, permanent/nonlinear media OK"
)

-- Geometric tolerance must not swap top and bottom on a very thin rectangle.
local thin = f.new({ xmin = 0, xmax = 1e5, ymin = 0, ymax = 1e-6, unit = 1, h = 5e4 })
f.boundary(thin, "all", { type = "neumann" })
f.boundary(thin, "bottom", {})
f.boundary(thin, "top", { value = 1e-6 })
f.solve(thin)
for _, e in ipairs(thin.triangles) do
    near(e.bx, 1, 1e-12)
    near(e.by, 0, 1e-12)
end

-- Manufactured quadratic potential with a mixed boundary. Robin includes
-- nonconstant edge data; refinement must reduce the field error.
local function manufactured(h)
    local m = f.new({
        xmin = 0,
        xmax = 1,
        ymin = 0,
        ymax = 1,
        unit = 1,
        h = h,
        source = function()
            return -2 / f.mu0
        end,
    })
    f.boundary(m, "left", { type = "neumann" })
    f.boundary(m, "right", { type = "robin", coefficient = 1 / f.mu0, value = -3 / f.mu0 })
    -- A=x^2 is not affine on the top/bottom, so use natural zero data there.
    f.boundary(m, "bottom", { type = "neumann" })
    f.boundary(m, "top", { type = "neumann" })
    f.solve(m)
    local error = 0
    for _, e in ipairs(m.triangles) do
        local x = 0
        for _, id in ipairs(e.ids) do
            x = x + m.nodes[id][1] / 3
        end
        -- Include exact within-element variance of the analytic gradient.
        local variance = 0
        for _, id in ipairs(e.ids) do
            variance = variance + (m.nodes[id][1] - x) ^ 2 / 12
        end
        error = error + e.area * (e.bx ^ 2 + (e.by + 2 * x) ^ 2 + 4 * variance)
    end
    return math.sqrt(error)
end
local coarse, fine = manufactured(0.125), manufactured(0.0625)
assert(coarse / fine > 1.8 and coarse / fine < 2.2)
print("Robin manufactured solution: first-order L2 field convergence OK")

-- A mirror-symmetric U device may be cut at x=0 with Ht=0. Compare its
-- nonlinear gap field with the full device, rather than just checking a residual.
local function symmetric_u(half)
    local m = f.new({ xmin = half and 0 or -85, xmax = 85, ymin = -75, ymax = 95, h = 2 })
    f.material(m, "iron", { library = "pure-iron" })
    local left = half and 0 or -40
    f.region(m, "iron", { { left, -30 }, { 40, -30 }, { 40, -15 }, { left, -15 } })
    f.region(m, "iron", { { 25, -15 }, { 40, -15 }, { 40, 30 }, { 25, 30 } })
    f.region(m, "iron", { { left, 33 }, { 40, 33 }, { 40, 48 }, { left, 48 } })
    if half then
        f.boundary(m, "left", { type = "neumann" })
    else
        f.region(m, "iron", { { -40, -15 }, { -25, -15 }, { -25, 30 }, { -40, 30 } })
        -- Preserve the same geometric anchor at the symmetry plane in both grids.
        f.region(m, "iron", { { 0, -30 }, { 40, -30 }, { 40, -15 }, { 0, -15 } })
    end
    local current = 3000 / (6 * 26 * 1e-6)
    local coils = { { 17, 23, -current }, { 42, 48, current } }
    if not half then
        coils[#coils + 1] = { -23, -17, -current }
        coils[#coils + 1] = { -48, -42, current }
    end
    for _, coil in ipairs(coils) do
        f.region(
            m,
            "air",
            { { coil[1], -10 }, { coil[2], -10 }, { coil[2], 16 }, { coil[1], 16 } },
            coil[3]
        )
    end
    f.solve(m)
    return m
end
local full_u, half_u = symmetric_u(false), symmetric_u(true)
local _, full_gap = f.sample(full_u, 32.5, 31.5)
local _, half_gap = f.sample(half_u, 32.5, 31.5)
assert(
    full_gap < -0.5 and math.abs(half_gap / full_gap - 1) < 0.01,
    string.format("full gap %.9g, half gap %.9g", full_gap, half_gap)
)
assert(#half_u.triangles < 0.55 * #full_u.triangles)
print(
    string.format(
        "Symmetric U: full/half gap By=%.6f/%.6f T; half-domain agreement <1%%",
        full_gap,
        half_gap
    )
)

-- Current balance must be checked even for tiny nonzero incompatible loads.
for _, current in ipairs({ 1e6, 1e-20 }) do
    local m = model("grid", 0.001)
    f.boundary(m, "all", { type = "neumann" })
    f.region(m, "air", { { 0, 0 }, { 20, 0 }, { 20, 12 }, { 0, 12 } }, current)
    fails(function()
        f.solve(m)
    end, "incompatible Neumann")
    assert(not m.nodes, "a failed preparation must not leave a usable partial model")
    f.boundary(m, "left", { type = "neumann", value = current * 0.02 / 2 })
    f.boundary(m, "right", { type = "neumann", value = current * 0.02 / 2 })
    f.solve(m)
    assert(m.gauge_node)
end
local zero = model("grid", 0.001)
f.boundary(zero, "all", { type = "robin", coefficient = 0 })
f.solve(zero)
assert(zero.gauge_node and zero.stats.peak == 0)
fails(function()
    f.boundary(zero, "left", {})
end, "before meshing")
local bad = model("delaunay", 0.001)
fails(function()
    f.boundary(bad, "left", { coefficient = 1 })
end, "requires type=robin")
fails(function()
    f.boundary(bad, "left", { type = "robin", coefficient = -1 })
end, "nonnegative")
fails(function()
    f.boundary(bad, "left", { dx = 0 / 0 })
end, "non-finite")
fails(function()
    f.boundary(bad, "left", { typo = 1 })
end, "unknown option")
fails(function()
    f.boundary(bad, "inner", {})
end, "unknown side")
fails(function()
    f.boundary(bad, "left", { type = "periodic" })
end, "unknown type")
f.boundary(bad, "left", { value = 0.01 })
fails(function()
    f.solve(bad)
end, "conflicting Dirichlet")
assert(not bad.nodes)
f.boundary(bad, "all", {})
f.solve(bad)
local callback = f.new({
    boundary = function()
        return 0
    end,
})
fails(function()
    f.boundary(callback, "all", {})
end, "callback")
print("Gauge compatibility, corner conflicts and invalid declarations OK")

-- Boundary physics changes reuse geometry but never a stale solution.
local lfs = require("lfs")
lfs.mkdir("build")
lfs.mkdir("build/tests")
local stem = "build/tests/boundary-cache"
local function cached(mode, kind, value)
    local m = model("delaunay", 0.001, { file = stem, mode = mode })
    if kind == "dirichlet" then
        f.boundary(m, "all", { value = value, dx = 0.3, dy = -0.7 })
    elseif kind == "neumann" then
        f.boundary(m, "all", { type = "neumann" })
        f.boundary(m, "bottom", { type = "neumann", value = value / f.mu0 })
        f.boundary(m, "top", { type = "neumann", value = -value / f.mu0 })
    else
        f.boundary(m, "all", { type = "robin", coefficient = 1e5 * value, value = -700 * value })
    end
    return m
end
for _, kind in ipairs({ "dirichlet", "neumann", "robin" }) do
    local cold = cached("refresh", kind, 0.2)
    f.solve(cold)
    local warm = cached("auto", kind, 0.2)
    f.solve(warm)
    assert(warm.cache_info.status == "solution")
    if kind == "robin" then
        for _, a in ipairs(warm.A) do
            near(a, 0.007, 1e-10)
        end
    end
    for i, a in ipairs(cold.A) do
        assert(a == warm.A[i])
    end
    local frozen = cached("frozen", kind, 0.2)
    f.solve(frozen)
    assert(frozen.cache_info.status == "solution")
    fails(function()
        f.solve(cached("frozen", kind, 0.3))
    end, "frozen cache unavailable")
    local changed = cached("auto", kind, 0.3)
    f.solve(changed)
    assert(changed.cache_info.status == "mesh")
end
local switched = cached("auto", "neumann", 0.5)
f.solve(switched)
assert(switched.cache_info.status == "mesh" and switched.gauge_node)
near(select(1, f.sample(switched, 8, 7)), 0.5, 1e-8)
os.remove(stem .. ".lfc")
print("Boundary caches: solution hits, gauge restoration and physical invalidation OK")
