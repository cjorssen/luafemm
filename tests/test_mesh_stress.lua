-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
package.path = "./?.lua;" .. package.path
local f = require("luafemm")
local cases = 0
local fallbacks = 0
local function run(name, regions, h, scale, offset)
    scale, offset = scale or 1, offset or 0
    local m = f.new({
        mesher = "delaunay",
        unit = 1,
        xmin = offset - 5 * scale,
        xmax = offset + 5 * scale,
        ymin = offset - 5 * scale,
        ymax = offset + 5 * scale,
        h = (h or 0.7) * scale,
    })
    for _, r in ipairs(regions) do
        local p = {}
        for _, q in ipairs(r) do
            p[#p + 1] = { offset + q[1] * scale, offset + q[2] * scale }
        end
        f.region(m, "air", p, 0, 0, r.mesh_size and r.mesh_size * scale or 0)
    end
    local ok, err = pcall(f.mesh, m)
    assert(ok, name .. ": " .. tostring(err))
    assert(m.mesh_stats.min_angle > 0)
    fallbacks = fallbacks + m.mesh_stats.exact_orientation + m.mesh_stats.exact_incircle
    cases = cases + 1
    return m
end
local rectangle = { { -3, -2 }, { 1, -2 }, { 1, 2 }, { -3, 2 } }
for _, scale in ipairs({ 1e-6, 1, 1e6 }) do
    run("scaled rectangle", { rectangle }, 0.7, scale)
end
run("translated rectangle", { rectangle }, 0.7, 1, 1e8)
for _, gap in ipairs({ 1e-2, 1e-5, 1e-8 }) do
    run("thin gap " .. gap, {
        { { -3, -3 }, { 0, -3 }, { 0, 3 }, { -3, 3 } },
        {
            { gap, -3 },
            { 3, -3 },
            { 3, 3 },
            {
                gap,
                3,
            },
        },
    })
end
run("almost straight vertex", { { { -3, -2 }, { 0, -2 + 1e-12 }, { 3, -2 }, { 3, 2 }, { -3, 2 } } })
run("shared diagonal", { { { -3, -3 }, { 3, 3 }, { -3, 3 } }, { { -3, -3 }, { 3, -3 }, { 3, 3 } } })
run(
    "T junction",
    { { { -3, -2 }, { 0, -2 }, { 0, 2 }, { -3, 2 } }, { { 0, -1 }, { 3, -1 }, { 3, 1 }, { 0, 1 } } }
)
run(
    "crossing diagonals",
    { { { -4, -1 }, { 4, 1 }, { 4, 2 }, { -4, 0 } }, { { -4, 1 }, { 4, -1 }, { 4, 0 }, { -4, 2 } } }
)
run("near-parallel intersections", {
    { { -4, 0 }, { 4, 1e-7 }, { 4, 1 }, { -4, 1 } },
    { { -4, 1e-7 }, { 4, 0 }, { 4, -1 }, { -4, -1 } },
})
local fine = { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } }
fine.mesh_size = 0.12
local refined = run("local refinement", { fine })
local total, n = 0, 0
for _, e in ipairs(refined.triangles) do
    local x, y = 0, 0
    for _, id in ipairs(e.ids) do
        x = x + refined.nodes[id][1] / 3
        y = y + refined.nodes[id][2] / 3
    end
    if math.abs(x) < 0.8 and math.abs(y) < 0.8 then
        total = total + e.area
        n = n + 1
    end
end
assert(n > 100 and total / n < 0.012, "local mesh size was ignored")
-- Reproducible convex/concave, partly overlapping polygons.
math.randomseed(73421)
for k = 1, 36 do
    local regions = {}
    for j = 1, 2 do
        local p = {}
        local n = 7 + k % 7
        for i = 1, n do
            local a = 2 * math.pi * i / n
            local r = 1 + math.random()
            p[i] = { r * math.cos(a) + (j - 1) * 0.8, r * math.sin(a) }
        end
        regions[#regions + 1] = p
    end
    run("random polygons " .. k, regions)
end
-- Invalid contours must fail explicitly, including adjacent retracing.
for _, p in ipairs({
    { { 0, 0 }, { 2, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } },
    {
        { 0, 0 },
        { 2, 2 },
        { 0, 2 },
        { 2, 0 },
    },
}) do
    assert(not pcall(f.region, f.new(), "air", p))
end
print(
    string.format(
        "%d mesh stress cases passed; %d exact-predicate fallbacks; invalid contours rejected",
        cases,
        fallbacks
    )
)
