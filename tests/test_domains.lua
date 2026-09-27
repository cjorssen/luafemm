-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local geometry = require("luafemm-geometry")
local function rectangle(x0, y0, x1, y1)
    return { { x0, y0 }, { x1, y0 }, { x1, y1 }, { x0, y1 } }
end
local outer, inner = rectangle(-4, -4, 4, 4), rectangle(-2, -2, 2, 2)
local function reverse(points)
    local out = {}
    for i = #points, 1, -1 do
        out[#out + 1] = points[i]
    end
    return out
end
local function near(a, b, tolerance)
    assert(math.abs(a - b) < (tolerance or 1e-8), string.format("%g != %g", a, b))
end
local function model(mesher, cache)
    local m = f.new({
        xmin = -5,
        xmax = 5,
        ymin = -5,
        ymax = 5,
        h = 0.7,
        mesher = mesher,
        cache = cache,
    })
    f.material(m, "shell", { mur = 100 })
    f.material(m, "core", { mur = 5, coercivity = 12000 })
    return m
end
local function areas(m)
    f.mesh(m)
    local area = { shell = 0, core = 0, air = 0 }
    for _, e in ipairs(m.triangles) do
        area[e.material] = area[e.material] + e.area / m.options.unit ^ 2
    end
    near(area.shell + area.core + area.air, 100)
    return area
end
for _, mesher in ipairs({ "grid", "delaunay" }) do
    for _, rule in ipairs({ "even odd", "nonzero" }) do
        for _, reversed in ipairs({ false, true }) do
            local m = model(mesher)
            -- Earlier core must survive the hole, independent of contour order.
            f.region(m, "core", inner, 0, 90)
            local hole = reversed and reverse(inner) or inner
            f.region_contours(m, "shell", { hole, outer }, 2e5, 0, 0.5, rule)
            local has_hole = rule == "even odd" or reversed
            local a = areas(m)
            near(a.shell, has_hole and 48 or 64)
            near(a.core, has_hole and 16 or 0)
            near(a.air, 36)
            f.solve(m)
            near(m.stats.net_current, a.shell * 1e-6 * 2e5)
            assert(m.stats.residual < 1e-7)
        end
    end
    -- Three nested loops alternate material/air/material under parity fill.
    local m = model(mesher)
    f.region_contours(
        m,
        "shell",
        { outer, inner, rectangle(-1, -1, 1, 1) },
        nil,
        nil,
        nil,
        "even odd"
    )
    near(areas(m).shell, 52)
    -- Two simple regions still use declaration order (no automatic subtraction).
    local painted = model(mesher)
    f.region(painted, "core", inner)
    f.region(painted, "shell", outer)
    near(areas(painted).core, 0)
    -- Distinct contours may intersect; parity removes their intersection.
    local intersecting = model(mesher)
    f.region_contours(intersecting, "shell", {
        rectangle(-4, -3, 1, 3),
        rectangle(-1, -3, 4, 3),
    }, nil, nil, nil, "even odd")
    near(areas(intersecting).shell, 36)
    local disjoint = model(mesher)
    f.region_contours(disjoint, "shell", {
        rectangle(-4, -1, -2, 1),
        rectangle(2, -1, 4, 1),
    })
    near(areas(disjoint).shell, 8)
end

-- Membership is invariant under translation, scale and reflection.
for _, scale in ipairs({ 1e-8, 1, 1e8, -3 }) do
    local loops = {}
    for i, loop in ipairs({ outer, inner }) do
        loops[i] = {}
        for j, p in ipairs(loop) do
            loops[i][j] = { 12 + p[1] * scale, -7 + p[2] * math.abs(scale) }
        end
    end
    local r = { contours = loops, fill_rule = "even odd" }
    assert(geometry.contains(r, 12 + 3 * scale, -7))
    assert(not geometry.contains(r, 12, -7))
end

-- Refinement belongs to the filled ring, not to its central air cavity.
local refined = model("delaunay")
f.region_contours(refined, "shell", { outer, inner }, 0, 0, 0.2, "even odd")
f.mesh(refined)
local cavity_nodes = 0
for _, p in ipairs(refined.nodes) do
    if math.abs(p[1]) < 0.0015 and math.abs(p[2]) < 0.0015 then
        cavity_nodes = cavity_nodes + 1
    end
end
assert(cavity_nodes > 0 and cavity_nodes < 40, "ring refinement leaked into the cavity")

-- Analytic shielding of a circular shell in a transverse uniform field.
-- Prescribe the exact exterior dipole potential on our rectangular boundary
-- to isolate polygonal/FE error from finite-domain truncation error.
local mur, a, b, imposed = 20, 0.002, 0.004, 0.5
local k = (mur - 1) / (mur + 1)
local d = 2 * mur * (mur + 1) * imposed / ((mur + 1) ^ 2 - (mur - 1) ^ 2 * (a / b) ^ 2)
local dipole = d * (b * b - k * a * a) - imposed * b * b
local shield = f.new({
    xmin = -10,
    xmax = 10,
    ymin = -10,
    ymax = 10,
    h = 0.5,
    mesher = "delaunay",
    boundary = function(x, y)
        return (imposed + dipole / (x * x + y * y)) * y
    end,
})
f.material(shield, "shell", { mur = mur })
local circles = {}
for i, r in ipairs({ 4, 2 }) do
    circles[i] = {}
    for j = 0, 95 do
        local angle = 2 * math.pi * j / 96
        circles[i][#circles[i] + 1] = { r * math.cos(angle), r * math.sin(angle) }
    end
end
f.region_contours(shield, "shell", circles, 0, 0, 0.35, "even odd")
f.solve(shield)
local bx, by = f.sample(shield, 0.13, 0.17)
local expected = 2 * d / (mur + 1)
assert(math.abs(bx / expected - 1) < 0.03)
assert(math.abs(by) < 0.003)
print(string.format("Circular shell: cavity Bx=%.6g T; analytic %.6g T (error <3%%)", bx, expected))

local paths = require("luafemm-path")
local frame = { 1, 0, 0, 1, 0, 0 }
local open = "\\pgfsyssoftpath@movetotoken{0pt}{0pt}\\pgfsyssoftpath@linetotoken{1pt}{0pt}"
assert(not pcall(paths.read_contours, open, frame, 0.01))
local triangle = open
    .. "\\pgfsyssoftpath@linetotoken{0pt}{1pt}\\pgfsyssoftpath@closepathtoken{0pt}{0pt}"
assert(#paths.read_contours(triangle .. triangle, frame, 0.01) == 2)
assert(
    not pcall(paths.read, triangle .. triangle, frame, 0.01, true),
    "profiles must remain connected"
)
assert(not pcall(paths.read_contours, triangle .. open, frame, 0.01), "every subpath must close")
local zero_open = "\\pgfsyssoftpath@movetotoken{0pt}{0pt}\\pgfsyssoftpath@linetotoken{0pt}{0pt}"
assert(
    not pcall(paths.read_contours, triangle .. zero_open, frame, 0.01),
    "do not discard open degenerate subpaths"
)

local lfs = require("lfs")
lfs.mkdir("build")
lfs.mkdir("build/tests")
local stem = "build/tests/domains-cache"
local function cached(mode, rule, size)
    local m = model("delaunay", { file = stem, mode = mode })
    f.region(m, "core", inner)
    f.region_contours(m, "shell", { outer, rectangle(-size, -size, size, size) }, 0, 0, 0.5, rule)
    return m
end
local cold = cached("refresh", "even odd", 2)
f.solve(cold)
local warm = cached("frozen", "even odd", 2)
f.solve(warm)
assert(warm.cache_info.status == "solution")
for id, a in ipairs(cold.A) do
    assert(a == warm.A[id])
end
local changed = cached("auto", "nonzero", 2)
f.solve(changed)
assert(changed.cache_info.status == "miss")
near(areas(changed).shell, 64)
local resized = cached("auto", "nonzero", 1)
f.solve(resized)
assert(resized.cache_info.status == "miss")
os.remove(stem .. ".lfc")

local m = model("grid")
assert(not pcall(f.region_contours, m, "shell", {}, nil, nil, nil, "even odd"))
assert(not pcall(f.region_contours, m, "shell", { outer }, nil, nil, nil, "unknown"))
assert(
    not pcall(f.region_contours, m, "shell", { outer, { { 0, 0 }, { 1, 1 }, { 0, 1 }, { 1, 0 } } })
)
local owned = rectangle(-1, -1, 1, 1)
f.region_contours(m, "shell", { owned })
owned[1][1] = -100
assert(m.regions[1].contours[1][1][1] == -1)
print("Compound domains: winding/parity, nesting, overlap, refinement, current and cache OK")
