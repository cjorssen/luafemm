-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local C = require("luafemm-components")
local function near(a, b)
    assert(math.abs(a - b) < 1e-8, tostring(a) .. " != " .. tostring(b))
end
local function area(points)
    return 1 / C.current(points, 1, 1)
end
local core = C.make("u core")
assert(#core.regions == 1 and not core.anchors["left gap"])
near(area(core.regions[1].points), 2550)
near(core.anchors["left pole"][1], -32.5)
near(core.anchors["left pole"][2], 30)
local options = { ampere_turns = -1234, gap_mesh_size = 0.5 }
local m = C.make("u electromagnet", options)
options.ampere_turns = 0
assert(m.options.ampere_turns == -1234 and #m.regions == 6)
near(m.half_width, 44)
near(m.half_height, 39)
near(m.anchors["left armature"][2] - m.anchors["left pole"][2], 3)
near(area(m.regions[2].points), 1200)
near(area(m.regions[3].points), 156)
near(area(m.regions[4].points), 156)
near(area(m.regions[5].points), 45)
assert(m.regions[3].turns == -1234 and m.regions[4].turns == 1234)
for _, angle in ipairs({ 0, 30, 90, 179 }) do
    local a = math.rad(angle)
    for _, scale in ipairs({ 0.5, 1, 2 }) do
        local points = {}
        for i, p in ipairs(m.regions[3].points) do
            points[i] = {
                100 + scale * (math.cos(a) * p[1] - math.sin(a) * p[2]),
                -70 + scale * (math.sin(a) * p[1] + math.cos(a) * p[2]),
            }
        end
        near(C.current(points, 0.001, 1234) * area(points) * 1e-6, 1234)
    end
end
assert(#C.make("u electromagnet", { gap_mesh_size = 0 }).regions == 4)
for _, o in ipairs({
    { width = 30 },
    { height = 15 },
    { leg_thickness = -1 },
    { gap = 0 },
    { coil_clearance = -1 },
    { coil_width = 100 },
    { coil_height = 100 },
    { coil_center_y = 40 },
    { armature_thickness = 0 },
    { width = math.huge },
    { width = 0 / 0 },
    { ampere_turns = "6000" },
    { mesh_size = -1 },
    { gap_mesh_size = -1 },
    { typo = 3 },
}) do
    assert(not pcall(C.make, "u electromagnet", o), "invalid component was accepted")
end
assert(not pcall(C.make, "rotor"))
assert(not pcall(C.current, { { 0, 0 }, { 1, 0 }, { 2, 0 } }, 0.001, 1))
print("Component areas, anchors, signed ampere-turns, transformations and validation passed.")

-- Extended shapes must retain valid simple part polygons and explicit cycles.
local mesh = require("luafemm-mesh")
for _, kind in ipairs({ "e core", "e electromagnet", "toroid", "tapered toroid" }) do
    local c = C.make(kind, { ideal = 1, magnet_span = 30 })
    for _, r in ipairs(c.regions) do
        mesh.validate(r.points)
        assert(area(r.points) > 0)
    end
    for _, cycle in pairs(c.cycles) do
        near(cycle[1][1], cycle[#cycle][1])
        near(cycle[1][2], cycle[#cycle][2])
    end
end
for _, options in ipairs({
    { radius = 4 },
    { gap_angle = -1 },
    { gap_angle = 90 },
    { pole_thickness = 11 },
    { taper_angle = 90 },
    { magnet_span = 180 },
    { arc_step = 1e-12 },
}) do
    assert(not pcall(C.make, "tapered toroid", options), "invalid toroid accepted")
end
for _, options in ipairs({
    { width = 25 },
    { left_gap = -0.5 },
    { height = 10 },
    { left_thickness = -1 },
    { coil_width = 50 },
}) do
    assert(not pcall(C.make, "e electromagnet", options), "invalid E accepted")
end
local closed = C.make("e electromagnet", { gap = 0 })
assert(#closed.regions == 4)
local permanent = C.make("toroid", { ampere_turns = 0, magnet_span = 30 })
assert(not permanent.anchors["coil positive"])
print("Extended E/toroid polygons, closed cycles, magnet-only shapes and bounds passed.")
