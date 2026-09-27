-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Flux identities are discrete invariants, not mesh-convergence assertions.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local C = require("luafemm-components")
local I = require("luafemm-ideal")
local P = require("luafemm-profiles")
local function near(a, b, tolerance)
    assert(
        math.abs(a - b) <= tolerance * math.max(math.abs(a), math.abs(b), 1e-20),
        string.format("%g != %g", a, b)
    )
end
local function model(kind, options, problem, nonlinear)
    options = options or {}
    options.ideal = 1
    local c = C.make(kind, options)
    local o = { field_model = "ideal", mesher = "delaunay", h = 3, depth = 10 }
    for k, v in pairs(problem or {}) do
        o[k] = v
    end
    local m = f.new(o)
    for _, role in ipairs({ "core", "armature" }) do
        f.material(m, role, nonlinear and { library = "pure-iron" } or { mur = 1000 })
    end
    f.material(m, "coil", { mur = 1 })
    f.material(m, "gap", { mur = 1 })
    f.material(m, "magnet", { library = "cast-alnico-5-lng37" })
    local center = c.anchors["core center"]
    for _, r in ipairs(c.regions) do
        local p = {}
        for _, v in ipairs(r.points) do
            p[#p + 1] = { v[1] - center[1], v[2] - center[2] }
        end
        local j = r.turns and C.current(p, m.options.unit, r.turns) or 0
        f.region(m, r.role, p, j, r.angle or 0, r.mesh_size).ideal_domain = r.role ~= "coil"
    end
    m.component_cycles = { M = {} }
    for name, points in pairs(c.cycles) do
        local p = {}
        for _, v in ipairs(points) do
            p[#p + 1] = { v[1] - center[1], v[2] - center[2] }
        end
        m.component_cycles.M[name] = p
    end
    return m, c
end
local function confinement(m)
    I.validate_potential(
        m,
        m.nodes,
        (function()
            local t = {}
            for _, e in ipairs(m.triangles) do
                t[#t + 1] = e.ids
            end
            return t
        end)(),
        m.A
    )
    for _, path in ipairs(f.contours(m, 13)) do
        local p, q = path[1], path[#path]
        near((p.x - q.x) ^ 2 + (p.y - q.y) ^ 2, 0, 1e-8)
    end
end
local toroid = model("toroid", { gap_angle = 0, ampere_turns = 100, arc_step = 3 }, { h = 1.3 })
f.solve(toroid)
local expected = f.mu0 * 1000 * 100 / (2 * math.pi) * math.log(35 / 25) * 0.01
local flux = f.flux(toroid, { { 25, 0 }, { 35, 0 } })
near(-flux, expected, 0.012)
for _, angle in ipairs({ 30, 75, 120, 210, 315 }) do
    -- Angles coincide with polygon vertices; endpoints lie on both walls.
    local a = math.rad(angle)
    near(
        f.flux(
            toroid,
            { { 25 * math.cos(a), 25 * math.sin(a) }, { 35 * math.cos(a), 35 * math.sin(a) } }
        ),
        flux,
        1e-9
    )
end
local mean = I.mean(toroid, "M", "main", "geometric")
P.register(toroid, "ideal-ampere", mean, { samples = 31 })
P.coordinates("ideal-ampere", "circulation")
local rows = P.sample("ideal-ampere")
near(rows[#rows].circulation, 100, 0.025)
assert(#I.mean(toroid, "M", "main", "flux") > 10)
confinement(toroid)
print(
    "Confined annulus: analytical flux, Ampere circulation and exact section conservation passed."
)
near(f.flux(toroid, { { 25.0005, 0 }, { 35.0005, 0 } }, 0.001), flux, 1e-9)
assert(not pcall(f.flux, toroid, { { 25, 0 }, { 35.01, 0 } }, 0.001))
assert(not pcall(f.flux, toroid, { { 25, 0 }, { 35, 0 } }, -1))
local seeded = model("toroid", { gap_angle = 0, ampere_turns = 0, arc_step = 3 }, { h = 1.3 })
f.excitation(seeded, 0, 0, 40)
f.excitation(seeded, 1, 1, 60)
f.solve(seeded)
near(-f.flux(seeded, { { 25, 0 }, { 35, 0 } }), expected, 0.012)
local u = model(
    "u electromagnet",
    { ampere_turns = 1200, left_thickness = 12, right_thickness = 20, yoke_thickness = 10 },
    nil,
    true
)
f.solve(u)
local a = f.flux(u, { { -28, 0 }, { -40, 0 } })
near(a, f.flux(u, { { 0, -20 }, { 0, -30 } }), 1e-9)
near(a, f.flux(u, { { 20, 0 }, { 40, 0 } }), 1e-9)
near(a, f.flux(u, { { -28, 31 }, { -40, 31 } }), 1e-9)
confinement(u)
print("Nonlinear unequal-section U: iron, yoke and gap carry identical flux.")
local e = model(
    "e electromagnet",
    { ampere_turns = 900, left_gap = 1, right_gap = 4, center_gap = 2, coil_width = 3 },
    nil,
    true
)
f.solve(e)
local left = f.flux(e, { { -40, 0 }, { -25, 0 } })
local right = f.flux(e, { { 25, 0 }, { 40, 0 } })
local center = f.flux(e, { { 7.5, 0 }, { -7.5, 0 } })
near(center, left + right, 1e-9)
assert(left > right and right > 0)
near(left, f.flux(e, { { -40, 33.5 }, { -25, 33.5 } }), 1e-9)
near(right, f.flux(e, { { 25, 33.5 }, { 40, 33.5 } }), 1e-9)
confinement(e)
assert(#I.mean(e, "M", "left", "geometric") == 5)
assert(#I.mean(e, "M", "right", "flux") > 5)
print("Asymmetric nonlinear E: branch flux balance and both confined gaps passed.")
local magnet =
    model("tapered toroid", { ampere_turns = 0, magnet_span = 30, arc_step = 5 }, nil, true)
f.solve(magnet)
assert(magnet.stats.peak > 0)
local a = math.rad(29)
local phi = f.flux(
    magnet,
    { { 25 * math.cos(a), 25 * math.sin(a) }, { 35 * math.cos(a), 35 * math.sin(a) } }
)
local theta = math.rad(4)
near(
    f.flux(magnet, {
        { 27 * math.cos(theta), 27 * math.sin(theta) },
        { 33 * math.cos(theta), 33 * math.sin(theta) },
    }),
    phi,
    1e-9
)
near(f.flux(magnet, { { 27 * math.cos(theta), 0 }, { 33 * math.cos(theta), 0 } }), phi, 1e-9)
confinement(magnet)
print("Tapered permanent-magnet toroid: pole, gap and core flux agree.")
-- Reuse preserves condensed boundary degrees; changed excitation re-solves.
local stem = "build/ideal-regression"
local m = model("toroid", { ampere_turns = 100 }, { cache = { mode = "refresh", file = stem } })
f.solve(m)
local saved = f.flux(m, { { -25, 0 }, { -35, 0 } })
local warm = model("toroid", { ampere_turns = 100 }, { cache = { mode = "frozen", file = stem } })
f.solve(warm)
assert(warm.cache_info.status == "solution")
near(f.flux(warm, { { -25, 0 }, { -35, 0 } }), saved, 1e-12)
local changed = model("toroid", { ampere_turns = 200 }, { cache = { mode = "auto", file = stem } })
f.solve(changed)
assert(changed.cache_info.status == "mesh")
near(f.flux(changed, { { -25, 0 }, { -35, 0 } }), 2 * saved, 1e-8)
assert(not pcall(f.export_fem, changed, "build/unsupported-ideal.fem"))
assert(not pcall(f.export_ans, changed, "build/unsupported-ideal.ans"))
assert(not pcall(f.import_fem, "tests/fixtures/uniform.fem", { field_model = "ideal" }))
local empty = f.new({ field_model = "ideal", mesher = "delaunay" })
assert(not pcall(f.solve, empty))
local solid = f.new({ field_model = "ideal", mesher = "delaunay", h = 5 })
f.region(solid, "air", { { -10, -10 }, { 10, -10 }, { 10, 10 }, { -10, 10 } }).ideal_domain = true
f.excitation(solid, 0, 0, 10)
assert(not pcall(f.solve, solid))
local bad = model("toroid", { ampere_turns = 0 })
f.excitation(bad, 80, 80, 100)
assert(not pcall(f.solve, bad))
print("Ideal caches, excitation validation and unsupported-export rejection passed.")
