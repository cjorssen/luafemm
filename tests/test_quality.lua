-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local function near(a, b, tol)
    assert(math.abs(a - b) <= tol, string.format("%.16g != %.16g", a, b))
end
local function fails(fn, pattern)
    local ok, message = pcall(fn)
    assert(not ok and tostring(message):find(pattern, 1, true), tostring(message))
end
local function model(changes, scale, offset)
    scale, offset = scale or 1, offset or 0
    local o = { mesher = "delaunay", unit = 1, xmin = -2, xmax = 2, ymin = -2, ymax = 2, h = 0.8 }
    for k, v in pairs(changes or {}) do
        o[k] = v
    end
    for _, k in ipairs({ "xmin", "xmax", "ymin", "ymax" }) do
        o[k] = offset + o[k] * scale
    end
    o.h = o.h * scale
    o.max_area = o.max_area and o.max_area * scale * scale
    return f.new(o)
end
local polygon = { { -0.7, -0.2 }, { 0.9, -0.3 }, { 0.8, 0.4 }, { -0.6, 0.3 } }
local function region(m, scale, offset)
    local p = {}
    for i, q in ipairs(polygon) do
        p[i] = { (offset or 0) + q[1] * (scale or 1), (offset or 0) + q[2] * (scale or 1) }
    end
    f.material(m, "test", { mur = 1 })
    f.region(m, "test", p)
end
-- Recompute quality from the returned, renumbered physical mesh, independently
-- of the refinement queue and diagnostics; also integrate the material area.
local function audit(m)
    local smallest, largest, material, total = 180, 0, 0, 0
    for _, t in ipairs(m.triangles) do
        total = total + t.area
        if t.material == "test" then
            material = material + t.area
        end
        largest = math.max(largest, t.area / m.options.unit ^ 2)
        for j = 1, 3 do
            local a, b, c =
                m.nodes[t.ids[j]], m.nodes[t.ids[j % 3 + 1]], m.nodes[t.ids[(j + 1) % 3 + 1]]
            local x, y, u, v = b[1] - a[1], b[2] - a[2], c[1] - a[1], c[2] - a[2]
            local angle = math.deg(math.atan(math.abs(x * v - y * u), x * u + y * v))
            smallest = math.min(smallest, angle)
        end
    end
    assert(smallest >= m.options.min_angle - 1e-6)
    assert(m.options.max_area == 0 or largest <= m.options.max_area * (1 + 1e-8))
    near(smallest, m.mesh_stats.min_angle, 1e-6)
    near(
        total,
        (m.options.xmax - m.options.xmin) * (m.options.ymax - m.options.ymin) * m.options.unit ^ 2,
        total * 1e-10
    )
    assert(m.mesh_stats.steiner_points == m.mesh_stats.segment_splits + m.mesh_stats.circumcenters)
    return material
end
local base = model()
region(base)
f.mesh(base)
assert(base.mesh_stats.min_angle < 20 and base.mesh_stats.steiner_points == 0)
for _, options in ipairs({
    { min_angle = 20 },
    { max_area = 0.03 },
    { min_angle = 25, max_area = 0.03 },
}) do
    local m = model(options)
    region(m)
    f.mesh(m)
    near(audit(m), 0.9, 1e-12)
    assert(#m.nodes > #base.nodes and m.mesh_stats.steiner_points > 0)
    assert(m.mesh_stats.circumcenters > 0)
    if options.max_area then
        assert(m.mesh_stats.segment_splits > 0)
    end
end
for _, config in ipairs({ { 1e-6, 0 }, { 1e6, 0 }, { 1, 1e6 } }) do
    local scale, offset = table.unpack(config)
    local m = model({ min_angle = 20, max_area = 0.2 }, scale, offset)
    region(m, scale, offset)
    f.mesh(m)
    near(audit(m) / scale ^ 2, 0.9, 1e-8)
end
-- Compound cavity and a separately refined thin air gap preserve their areas.
local m = model({ min_angle = 20, h = 0.6 })
f.material(m, "test", { mur = 2 })
f.region_contours(m, "test", {
    { { -1.5, -1.5 }, { 1.5, -1.5 }, { 1.5, 1.5 }, { -1.5, 1.5 } },
    { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } },
}, 0, 0, 0, "even odd")
f.region(m, "air", { { -1, -0.05 }, { 1, -0.05 }, { 1, 0.05 }, { -1, 0.05 } }, 0, 0, 0.2)
f.mesh(m)
near(audit(m), 5, 1e-10)
-- Segment splitting must retain the exterior tags needed by the weak Robin
-- assembly. Az=.2*y gives Bx=.2, By=0 in the uniform test material and air.
local patch = model({ min_angle = 25, max_area = 0.03 })
region(patch)
f.boundary(patch, "left", { type = "dirichlet", dy = 0.2 })
f.boundary(patch, "right", { type = "dirichlet", dy = 0.2 })
f.boundary(patch, "bottom", { type = "robin", coefficient = 2, value = 0.2 / f.mu0 + 0.8 })
f.boundary(patch, "top", { type = "robin", coefficient = 2, value = -0.2 / f.mu0 - 0.8 })
f.solve(patch)
audit(patch)
for i, p in ipairs(patch.nodes) do
    near(patch.A[i], 0.2 * p[2], 2e-7)
end
-- A disk cache round-trip preserves the refined mesh. Every quality setting
-- participates in its descriptor, including the finite insertion budget.
local stem = "build/tests/quality-cache"
local function cached(mode, target)
    local q = model({
        min_angle = target or 20,
        max_area = 0.2,
        cache = { file = stem, mode = mode },
    })
    region(q)
    f.solve(q)
    return q
end
local cold, warm = cached("refresh"), cached("frozen")
assert(warm.cache_info.status == "solution")
for i, value in ipairs(cold.A) do
    assert(warm.A[i] == value)
end
assert(cold.mesh_stats.steiner_points == warm.mesh_stats.steiner_points)
fails(function()
    cached("frozen", 25)
end, "frozen cache unavailable")
for _, change in ipairs({ { max_area = 0.1 }, { max_steiner = 6000 } }) do
    local q = model({ min_angle = 20, max_area = 0.2, cache = { file = stem, mode = "frozen" } })
    for k, v in pairs(change) do
        q.options[k] = v
    end
    region(q)
    fails(function()
        f.solve(q)
    end, "frozen cache unavailable")
end
-- Hard limits never return an apparently successful under-refined result.
fails(function()
    local q = model({ max_area = 0.001, max_steiner = 0 })
    f.mesh(q)
end, "Steiner point budget exhausted")
fails(function()
    local q = model({ min_angle = 20, max_steiner = 20 })
    f.region(q, "air", { { -1, 0 }, { 1, -0.01 }, { 1, 0.01 } })
    f.mesh(q)
end, "Steiner point budget exhausted")
for _, options in ipairs({
    { min_angle = -1 },
    { min_angle = 31 },
    { max_area = -1 },
    { max_area = math.huge },
    { max_steiner = -1 },
    { max_steiner = 2.5 },
    { mesher = "grid", min_angle = 20 },
}) do
    assert(not pcall(model, options))
end
print("Quality targets, geometry, affine Robin patch, caches and explicit limits passed.")
