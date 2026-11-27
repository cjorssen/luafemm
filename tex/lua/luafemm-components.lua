-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Parametric components: one geometric description for paths and anchors.
-- No TeX state, materials, mesh or solver is accessed here. Coordinates are
-- model units; the TeX interface fixes those units to millimetres.
-- @module luafemm-components
local C = {}
local defaults = {
    width = 80,
    left_thickness = 0,
    right_thickness = 0,
    yoke_thickness = 0,
    ideal = 0,
    height = 60,
    leg_thickness = 15,
    gap = 3,
    armature_thickness = 15,
    coil_width = 6,
    coil_height = 26,
    coil_clearance = 2,
    coil_center_y = 3,
    ampere_turns = 6000,
    mesh_size = 0,
    gap_mesh_size = 1,
    magnetization_angle = 0,
}
for key, value in pairs(require("luafemm-shapes").defaults) do
    defaults[key] = value
end
local function check(value, message)
    assert(value, "luafemm component: " .. message)
end
local function rectangle(x0, y0, x1, y1)
    return { { x0, y0 }, { x1, y0 }, { x1, y1 }, { x0, y1 } }
end

--- Construct a U core or a complete left-wound electromagnet.
-- Every region has one simple closed polygon. Later air refinement regions
-- occupy only the gaps, so no compound paths or domain subtraction is needed.
-- @tparam string kind "u core" or "u electromagnet".
-- @tparam[opt] table options Numeric geometry, mesh and excitation parameters.
-- @treturn table Independent options, regions, anchors and centred bounds.
function C.make(kind, options)
    if kind == "rotor" or kind == "stator" then
        return require("luafemm-machines").make(kind, options)
    end
    if kind ~= "u core" and kind ~= "u electromagnet" then
        return require("luafemm-shapes").make(kind, options)
    end
    local o = {}
    for key, value in pairs(options or {}) do
        check(defaults[key] ~= nil, "unknown option " .. tostring(key))
        check(
            type(value) == "number" and value == value and math.abs(value) < math.huge,
            "expected a finite number for " .. key
        )
    end
    for key, value in pairs(defaults) do
        o[key] = options and options[key] or value
    end
    local w, h, t = o.width, o.height, o.leg_thickness
    local tl = o.left_thickness > 0 and o.left_thickness or t
    local tr = o.right_thickness > 0 and o.right_thickness or t
    local tb = o.yoke_thickness > 0 and o.yoke_thickness or t
    check(
        o.left_thickness >= 0 and o.right_thickness >= 0 and o.yoke_thickness >= 0,
        "section thickness overrides must be nonnegative"
    )
    check(
        t > 0 and w > tl + tr and h > tb,
        "require leg thickness > 0, width > 2*thickness and height > thickness"
    )
    check(o.mesh_size >= 0 and o.gap_mesh_size >= 0, "mesh sizes must be nonnegative")
    local x0, x1, y0, y1 = -w / 2, w / 2, -h / 2, h / 2
    local left, right = x0 + tl / 2, x1 - tr / 2
    local regions = {
        {
            role = "core",
            mesh_size = o.mesh_size,
            points = {
                { x0, y0 },
                { x1, y0 },
                { x1, y1 },
                { x1 - tr, y1 },
                { x1 - tr, y0 + tb },
                { x0 + tl, y0 + tb },
                { x0 + tl, y1 },
                { x0, y1 },
            },
        },
    }
    local anchors = {
        ["core center"] = { 0, 0 },
        ["window center"] = { 0, t / 2 },
        ["left pole"] = { left, y1 },
        ["right pole"] = { right, y1 },
        ["left pole outer"] = { x0, y1 },
        ["left pole inner"] = { x0 + tl, y1 },
        ["right pole inner"] = { x1 - tr, y1 },
        ["right pole outer"] = { x1, y1 },
        ["yoke center"] = { 0, y0 + tb / 2 },
    }
    if kind == "u electromagnet" then
        local g, a = o.gap, o.armature_thickness
        local cw, ch, c, cy = o.coil_width, o.coil_height, o.coil_clearance, o.coil_center_y
        check(g > 0 and a > 0, "gap and armature thickness must be positive")
        check(
            cw > 0 and ch > 0 and c >= 0,
            "coil dimensions must be positive; clearance must be nonnegative"
        )
        check(
            cy - ch / 2 >= y0 + tb and cy + ch / 2 <= y1,
            "coil must fit along the free part of the left leg"
        )
        check(c + cw <= w - tl - tr, "inner coil section does not fit in the window")
        regions[#regions + 1] = {
            role = "armature",
            mesh_size = o.mesh_size,
            points = rectangle(x0, y1 + g, x1, y1 + g + a),
        }
        regions[#regions + 1] = {
            role = "coil",
            mesh_size = o.mesh_size,
            turns = o.ampere_turns,
            points = rectangle(x0 - c - cw, cy - ch / 2, x0 - c, cy + ch / 2),
        }
        regions[#regions + 1] = {
            role = "coil",
            mesh_size = o.mesh_size,
            turns = -o.ampere_turns,
            points = rectangle(x0 + tl + c, cy - ch / 2, x0 + tl + c + cw, cy + ch / 2),
        }
        if o.gap_mesh_size > 0 or o.ideal == 1 then
            regions[#regions + 1] = {
                role = "gap",
                mesh_size = o.gap_mesh_size,
                points = rectangle(x0, y1, x0 + tl, y1 + g),
            }
            regions[#regions + 1] = {
                role = "gap",
                mesh_size = o.gap_mesh_size,
                points = rectangle(x1 - tr, y1, x1, y1 + g),
            }
        end
        anchors["left armature"], anchors["right armature"] = { left, y1 + g }, { right, y1 + g }
        anchors["left gap"], anchors["right gap"] = { left, y1 + g / 2 }, { right, y1 + g / 2 }
        anchors["armature center"] = { 0, y1 + g + a / 2 }
        anchors["coil positive"] = { x0 - c - cw / 2, cy }
        anchors["coil negative"] = { x0 + tl + c + cw / 2, cy }
        x0, y1 = x0 - c - cw, y1 + g + a
    end
    -- Standard compass anchors describe the envelope of all component parts.
    -- Physical anchors retain their meaning even when that envelope changes.
    local cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    anchors.center = { cx, cy }
    anchors.north, anchors.south = { cx, y1 }, { cx, y0 }
    anchors.east, anchors.west = { x1, cy }, { x0, cy }
    anchors["north east"], anchors["north west"] = { x1, y1 }, { x0, y1 }
    anchors["south east"], anchors["south west"] = { x1, y0 }, { x0, y0 }
    local cycles = {}
    if kind == "u electromagnet" then
        cycles.main = {
            { left, y0 + tb / 2 },
            { left, h / 2 + o.gap + o.armature_thickness / 2 },
            { right, h / 2 + o.gap + o.armature_thickness / 2 },
            { right, y0 + tb / 2 },
            { left, y0 + tb / 2 },
        }
    end
    for _, points in pairs(cycles) do
        for _, p in ipairs(points) do
            p[1], p[2] = p[1] - cx, p[2] - cy
        end
    end
    for _, r in ipairs(regions) do
        for _, p in ipairs(r.points) do
            p[1], p[2] = p[1] - cx, p[2] - cy
        end
    end
    for _, p in pairs(anchors) do
        p[1], p[2] = p[1] - cx, p[2] - cy
    end
    return {
        kind = kind,
        options = o,
        regions = regions,
        cycles = cycles,
        anchors = anchors,
        half_width = (x1 - x0) / 2,
        half_height = (y1 - y0) / 2,
    }
end

--- Convert ampere-turns to uniform Jz using the actual captured section.
-- @tparam table points Polygon in model coordinates, after transformations.
-- @tparam number unit Metres per model unit.
-- @tparam number turns Signed ampere-turns (integral of Jz over the section).
-- @treturn number Current density in A/m^2.
function C.current(points, unit, turns)
    -- Translate before summation to limit cancellation far from the origin.
    local x, y, twice = points[1][1], points[1][2], 0
    for i, p in ipairs(points) do
        local q = points[i % #points + 1]
        twice = twice + (p[1] - x) * (q[2] - y) - (q[1] - x) * (p[2] - y)
    end
    check(math.abs(twice) > 0, "zero-area winding section")
    return turns / (math.abs(twice) * 0.5 * unit * unit)
end
return C
