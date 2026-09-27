-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- E and toroidal component geometry, with named cycles and pole anchors.
-- Toroidal components describe a planar annulus of constant extrusion depth.
-- Tapered poles are straight-sided in this plane, not an axisymmetric 3D solve.
-- @module luafemm-shapes
local S = {}
S.defaults = {
    radius = 30,
    core_thickness = 10,
    pole_thickness = 6,
    gap_angle = 8,
    taper_angle = 25,
    magnet_span = 0,
    arc_step = 5,
    left_gap = -1,
    right_gap = -1,
    center_gap = -1,
    center_thickness = 0,
}
local common = {
    width = 80,
    height = 60,
    leg_thickness = 15,
    left_thickness = 0,
    right_thickness = 0,
    yoke_thickness = 0,
    ideal = 0,
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
local function check(v, s)
    assert(v, "luafemm component: " .. s)
end
local function rect(a, b, c, d)
    return { { a, b }, { c, b }, { c, d }, { a, d } }
end
local function polar(r, a)
    return { r * math.cos(math.rad(a)), r * math.sin(math.rad(a)) }
end
local function finish(kind, o, regions, anchors, cycles, bounds, sources)
    local x0, y0, x1, y1 = table.unpack(bounds)
    local cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    anchors.center = { cx, cy }
    anchors.north, anchors.south = { cx, y1 }, { cx, y0 }
    anchors.east, anchors.west = { x1, cy }, { x0, cy }
    anchors["north east"], anchors["north west"] = { x1, y1 }, { x0, y1 }
    anchors["south east"], anchors["south west"] = { x1, y0 }, { x0, y0 }
    local function shift(p)
        p[1], p[2] = p[1] - cx, p[2] - cy
    end
    for _, r in ipairs(regions) do
        for _, p in ipairs(r.points) do
            shift(p)
        end
    end
    for _, points in pairs(cycles) do
        for _, p in ipairs(points) do
            shift(p)
        end
    end
    for _, p in pairs(anchors) do
        shift(p)
    end
    for _, p in ipairs(sources or {}) do
        shift(p)
    end
    return {
        kind = kind,
        options = o,
        regions = regions,
        anchors = anchors,
        cycles = cycles,
        sources = sources,
        half_width = (x1 - x0) / 2,
        half_height = (y1 - y0) / 2,
    }
end
local function ecore(kind, o)
    local w, h, t = o.width, o.height, o.leg_thickness
    local tl = o.left_thickness > 0 and o.left_thickness or t
    local tr = o.right_thickness > 0 and o.right_thickness or t
    local tc = o.center_thickness > 0 and o.center_thickness or t
    local tb = o.yoke_thickness > 0 and o.yoke_thickness or t
    local gl = o.left_gap >= 0 and o.left_gap or o.gap
    local gr = o.right_gap >= 0 and o.right_gap or o.gap
    local gc = o.center_gap >= 0 and o.center_gap or o.gap
    check(
        tl > 0
            and tr > 0
            and tc > 0
            and tb > 0
            and h > tb
            and w / 2 - tl > tc / 2
            and w / 2 - tr > tc / 2,
        "E core requires positive legs and two open windows"
    )
    check(
        gl >= 0 and gr >= 0 and gc >= 0 and o.armature_thickness > 0,
        "invalid E gaps or armature"
    )
    local x0, x1, y0 = -w / 2, w / 2, -h / 2
    local top = h / 2 + math.max(gl, gr, gc)
    local yl, yr, yc = top - gl, top - gr, top - gc
    check(math.min(yl, yr, yc) > y0 + tb, "gap consumes a complete E leg")
    local regions = {
        {
            role = "core",
            mesh_size = o.mesh_size,
            points = {
                { x0, y0 },
                { x1, y0 },
                { x1, yr },
                { x1 - tr, yr },
                { x1 - tr, y0 + tb },
                { tc / 2, y0 + tb },
                { tc / 2, yc },
                { -tc / 2, yc },
                { -tc / 2, y0 + tb },
                { x0 + tl, y0 + tb },
                { x0 + tl, yl },
                { x0, yl },
            },
        },
    }
    local left, right = x0 + tl / 2, x1 - tr / 2
    local anchors = {
        ["core center"] = { 0, 0 },
        ["window center"] = { (x0 + tl - tc / 2) / 2, 0 },
        ["left window"] = { (x0 + tl - tc / 2) / 2, 0 },
        ["right window"] = { (tc / 2 + x1 - tr) / 2, 0 },
        ["left pole"] = { left, yl },
        ["right pole"] = { right, yr },
        ["center pole"] = { 0, yc },
        ["left pole outer"] = { x0, yl },
        ["left pole inner"] = { x0 + tl, yl },
        ["right pole outer"] = { x1, yr },
        ["right pole inner"] = { x1 - tr, yr },
        ["yoke center"] = { 0, y0 + tb / 2 },
    }
    local cycles = {}
    local ymax = math.max(yl, yr, yc)
    if kind == "e electromagnet" then
        ymax = top + o.armature_thickness
        regions[#regions + 1] =
            { role = "armature", mesh_size = o.mesh_size, points = rect(x0, top, x1, ymax) }
        local cw, ch, cl, cy = o.coil_width, o.coil_height, o.coil_clearance, o.coil_center_y
        check(
            cw > 0 and ch > 0 and cl >= 0 and tc / 2 + cl + cw < math.min(w / 2 - tl, w / 2 - tr),
            "E winding does not fit in the windows"
        )
        check(
            cy - ch / 2 >= y0 + tb and cy + ch / 2 <= math.min(yl, yr, yc),
            "E winding exceeds leg height"
        )
        regions[#regions + 1] = {
            role = "coil",
            mesh_size = o.mesh_size,
            turns = o.ampere_turns,
            points = rect(-tc / 2 - cl - cw, cy - ch / 2, -tc / 2 - cl, cy + ch / 2),
        }
        regions[#regions + 1] = {
            role = "coil",
            mesh_size = o.mesh_size,
            turns = -o.ampere_turns,
            points = rect(tc / 2 + cl, cy - ch / 2, tc / 2 + cl + cw, cy + ch / 2),
        }
        for _, part in ipairs({
            { x0, yl, x0 + tl, gl, "left", left },
            { -tc / 2, yc, tc / 2, gc, "center", 0 },
            { x1 - tr, yr, x1, gr, "right", right },
        }) do
            if part[4] > 0 and (o.gap_mesh_size > 0 or o.ideal == 1) then
                regions[#regions + 1] = {
                    role = "gap",
                    mesh_size = o.gap_mesh_size,
                    points = rect(part[1], part[2], part[3], top),
                }
            end
            anchors[part[5] .. " gap"] = { part[6], (part[2] + top) / 2 }
            anchors[part[5] .. " armature"] = { part[6], top }
        end
        anchors["coil positive"] = { -tc / 2 - cl - cw / 2, cy }
        anchors["coil negative"] = { tc / 2 + cl + cw / 2, cy }
        anchors["armature center"] = { 0, (top + ymax) / 2 }
        local low, high = y0 + tb / 2, (top + ymax) / 2
        cycles.left = { { left, low }, { 0, low }, { 0, high }, { left, high }, { left, low } }
        cycles.right = { { 0, low }, { 0, high }, { right, high }, { right, low }, { 0, low } }
    end
    return finish(kind, o, regions, anchors, cycles, { x0, y0, x1, ymax })
end
local function toroid(kind, o)
    local r, t, g = o.radius, o.core_thickness, o.gap_angle
    local taper = kind == "tapered toroid" and o.taper_angle or 0
    local pole = taper > 0 and o.pole_thickness or t
    local span = o.magnet_span
    check(
        r > t / 2 and t > 0 and pole > 0 and pole <= t,
        "require radius > thickness/2 and 0 < pole thickness <= core thickness"
    )
    check(
        g >= 0 and g < 90 and taper >= 0 and taper < 90 and o.arc_step > 0 and o.arc_step <= 30,
        "invalid toroidal gap, taper or arc step"
    )
    check(span >= 0 and span < 180 and g + 2 * taper + span < 360, "magnet overlaps tapered poles")
    local a, b = g / 2, 360 - g / 2
    check(math.ceil(360 / o.arc_step) <= 4096, "arc step exceeds the 4096-edge geometry budget")
    local regions, mean = {}, {}
    local function section(angle, width)
        return polar(r + width / 2, angle), polar(r - width / 2, angle)
    end
    local function sector(first, last, role)
        local n = math.max(1, math.ceil((last - first) / o.arc_step))
        -- Magnet slices carry tangent magnetisation. The core remains one
        -- polygon per arc to avoid unnecessary internal mesh constraints.
        if role == "magnet" then
            for j = 1, n do
                local u, v = first + (last - first) * (j - 1) / n, first + (last - first) * j / n
                local ao, ai = section(u, t)
                local bo, bi = section(v, t)
                regions[#regions + 1] = {
                    role = role,
                    mesh_size = o.mesh_size,
                    angle = (u + v) / 2 + 90 + o.magnetization_angle,
                    points = { ao, bo, bi, ai },
                }
            end
        else
            local outer, inner = {}, {}
            for j = 0, n do
                local angle = first + (last - first) * j / n
                local po, pi = section(angle, t)
                outer[#outer + 1], inner[#inner + 1] = po, pi
            end
            for j = #inner, 1, -1 do
                outer[#outer + 1] = inner[j]
            end
            regions[#regions + 1] = { role = role, mesh_size = o.mesh_size, points = outer }
        end
        for j = 0, n do
            mean[#mean + 1] = polar(r, first + (last - first) * j / n)
        end
    end
    if taper > 0 then
        local ao, ai = section(a, pole)
        local bo, bi = section(a + taper, t)
        regions[#regions + 1] =
            { role = "core", mesh_size = o.mesh_size, points = { ao, bo, bi, ai } }
        mean[#mean + 1] = polar(r, a)
    end
    local start, stop = a + taper, b - taper
    if span > 0 then
        sector(start, 180 - span / 2, "core")
        sector(180 - span / 2, 180 + span / 2, "magnet")
        sector(180 + span / 2, stop, "core")
    else
        -- Split a closed annulus at pi to retain simple polygons.
        sector(start, 180, "core")
        sector(180, stop, "core")
    end
    if taper > 0 then
        local ao, ai = section(b - taper, t)
        local bo, bi = section(b, pole)
        regions[#regions + 1] =
            { role = "core", mesh_size = o.mesh_size, points = { ao, bo, bi, ai } }
        mean[#mean + 1] = polar(r, b)
    end
    if g > 0 and (o.gap_mesh_size > 0 or o.ideal == 1) then
        local ao, ai = section(b, pole)
        local bo, bi = section(a, pole)
        regions[#regions + 1] =
            { role = "gap", mesh_size = o.gap_mesh_size, points = { ao, bo, bi, ai } }
    end
    if (mean[#mean][1] - mean[1][1]) ^ 2 + (mean[#mean][2] - mean[1][2]) ^ 2 < 1e-20 then
        mean[#mean] = { mean[1][1], mean[1][2] }
    else
        mean[#mean + 1] = { mean[1][1], mean[1][2] }
    end
    local anchors = {
        ["core center"] = { 0, 0 },
        ["window center"] = { 0, 0 },
        ["upper pole"] = polar(r, a),
        ["lower pole"] = polar(r, b),
        ["upper pole inner"] = polar(r - pole / 2, a),
        ["upper pole outer"] = polar(r + pole / 2, a),
        ["lower pole inner"] = polar(r - pole / 2, b),
        ["lower pole outer"] = polar(r + pole / 2, b),
        ["gap center"] = { r * math.cos(math.rad(a)), 0 },
        ["magnet center"] = { -r, 0 },
    }
    local radius = r + t / 2
    local xmax = radius
    if o.ampere_turns ~= 0 then
        local cw, ch, cl = o.coil_width, o.coil_height, o.coil_clearance
        check(
            cw > 0 and ch > 0 and cl > 0 and math.sqrt((cw / 2) ^ 2 + (ch / 2) ^ 2) < r - t / 2,
            "toroidal winding must fit strictly inside the window"
        )
        regions[#regions + 1] = {
            role = "coil",
            mesh_size = o.mesh_size,
            turns = o.ampere_turns,
            points = rect(-cw / 2, -ch / 2, cw / 2, ch / 2),
        }
        xmax = radius + cl + cw
        anchors["coil positive"] = { 0, 0 }
        anchors["coil negative"] = { radius + cl + cw / 2, 0 }
        regions[#regions + 1] = {
            role = "coil",
            mesh_size = o.mesh_size,
            turns = -o.ampere_turns,
            points = rect(radius + cl, -ch / 2, xmax, ch / 2),
        }
    end
    return finish(kind, o, regions, anchors, { main = mean }, { -radius, -radius, xmax, radius })
end

--- Construct an E component or planar toroid with optional permanent-magnet arc.
-- All numeric dimensions use model units; angles are degrees.
function S.make(kind, options)
    check(
        kind == "e core"
            or kind == "e electromagnet"
            or kind == "toroid"
            or kind == "tapered toroid",
        "unknown component " .. tostring(kind)
    )
    local o = {}
    for k, v in pairs(common) do
        o[k] = v
    end
    for k, v in pairs(S.defaults) do
        o[k] = v
    end
    for k, v in pairs(options or {}) do
        check(o[k] ~= nil, "unknown option " .. tostring(k))
        check(type(v) == "number" and v == v and math.abs(v) < math.huge, "invalid " .. k)
        o[k] = v
    end
    check(o.mesh_size >= 0 and o.gap_mesh_size >= 0, "negative mesh spacing")
    for _, key in ipairs({
        "left_thickness",
        "right_thickness",
        "center_thickness",
        "yoke_thickness",
    }) do
        check(o[key] >= 0, "negative section thickness")
    end
    for _, key in ipairs({ "left_gap", "right_gap", "center_gap" }) do
        check(o[key] == -1 or o[key] >= 0, "gap override must be -1 or nonnegative")
    end
    if kind == "e core" or kind == "e electromagnet" then
        return ecore(kind, o)
    end
    return toroid(kind, o)
end
return S
