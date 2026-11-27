-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Owned spatial primitives and orthographic layered scenes, independent of physics.
-- Coordinates stay in model units until projection. Camera/layout operations never
-- invoke a field provider. Layers are explicit drawing order, NOT depth visibility.
local S = {}
local copy = require("luafemm-physical").copy
local function check(ok, message)
    assert(ok, "luafemm 3D view: " .. message)
end
local function number(x, name, default)
    if x == nil or x == "" then
        x = default
    end
    x = tonumber(x)
    check(x and x == x and math.abs(x) < math.huge, "finite " .. name .. " required")
    return x
end
local function choice(x, values, name, default)
    x = x or default
    check(values[x], "invalid " .. name)
    return x
end
local function vector(x, name)
    if type(x) == "string" then
        local a, b, c = x:match("^%s*([^/]+)/([^/]+)/([^/]+)%s*$")
        check(a, name .. " must be x/y/z")
        x = { a, b, c }
    end
    check(type(x) == "table" and #x == 3, name .. " must have three components")
    return { number(x[1], name), number(x[2], name), number(x[3], name) }
end
local function point(o, u, v, a, b)
    return { o[1] + a * u[1] + b * v[1], o[2] + a * u[2] + b * v[2], o[3] + a * u[3] + b * v[3] }
end
local function integer(x, name, default, lo, hi)
    x = number(x, name, default)
    check(x % 1 == 0 and x >= lo and x <= hi, name .. " outside its integer budget")
    return x
end
S.check, S.number, S.copy = check, number, copy
function S.new(unit)
    check(({ m = true, cm = true, mm = true, um = true })[unit], "unknown scene length unit")
    return { unit = unit, paths = {}, planes = {}, point_count = 0 }
end
--- Append owned paths transactionally; no caller-owned point array is retained.
function S.add_paths(scene, role, paths)
    choice(
        role,
        { geometry = true, field = true, vectors = true, source = true, equipotential = true },
        "primitive role"
    )
    local added, count = {}, scene.point_count
    for _, path in ipairs(paths) do
        check(type(path.points) == "table" and #path.points >= 2, "path needs two points")
        local p = { role = role, points = {} }
        if role == "geometry" then
            p.part = choice(
                path.part,
                { winding = true, leads = true, ["return"] = true },
                "source part",
                "winding"
            )
            p.radius = number(path.radius, "physical wire radius", 0)
            check(p.radius >= 0, "physical wire radius must be nonnegative")
            p.source = path.source and integer(path.source, "source index", nil, 1, 2000)
        end
        count = count + #path.points
        check(count <= 2000000, "scene point budget exceeded")
        for _, q in ipairs(path.points) do
            p.points[#p.points + 1] = vector(q, "path point")
        end
        added[#added + 1] = p
    end
    check(#scene.paths + #added <= 20000, "scene path budget exceeded")
    for _, p in ipairs(added) do
        scene.paths[#scene.paths + 1] = p
    end
    scene.point_count = count
end
--- A reference plane is a finite annotation o+a*u+b*v, never a boundary condition.
-- Styles/labels are caller-supplied TeX display payloads, never native file data.
function S.add_plane(scene, o)
    check(#scene.planes < 32, "plane budget exceeded")
    local mode =
        choice(o.plane, { xy = true, xz = true, yz = true, custom = true }, "reference plane", "xy")
    local bases = {
        xy = { { 1, 0, 0 }, { 0, 1, 0 } },
        xz = { { 1, 0, 0 }, { 0, 0, 1 } },
        yz = {
            { 0, 1, 0 },
            { 0, 0, 1 },
        },
    }
    local p = {
        origin = vector(o.origin or "0/0/0", "plane origin"),
        u = mode == "custom" and vector(o.u, "plane u") or bases[mode][1],
        v = mode == "custom" and vector(o.v, "plane v") or bases[mode][2],
        label = o.label or "",
        style = o.style or "",
        grid_style = o.grid_style or "",
        layer = integer(o.layer, "plane layer", 0, -1000, 1000),
        occlusion = choice(
            o.occlusion,
            { opaque = true, none = true, translucent = true },
            "plane occlusion",
            "opaque"
        ),
        opacity = number(o.opacity, "plane opacity", 0.25),
        grid = choice(o.grid, { none = true, major = true, both = true }, "plane grid", "none"),
        step = number(o.step, "grid step", 1),
        subdivisions = integer(o.subdivisions, "minor subdivisions", 5, 2, 10),
    }
    check(p.opacity > 0 and p.opacity <= 1, "plane opacity must be in (0,1]")
    check(p.occlusion ~= "translucent" or p.grid == "none", "translucent planes require grid=none")
    for _, a in ipairs({ "u", "v" }) do
        p[a .. "min"] = number(o[a .. "min"], a .. "min", -2)
        p[a .. "max"] = number(o[a .. "max"], a .. "max", 2)
        check(p[a .. "max"] > p[a .. "min"], "reversed plane extent")
    end
    local u, v = p.u, p.v
    local un = math.sqrt(u[1] ^ 2 + u[2] ^ 2 + u[3] ^ 2)
    local vn = math.sqrt(v[1] ^ 2 + v[2] ^ 2 + v[3] ^ 2)
    local cross = math.sqrt(
        (u[2] * v[3] - u[3] * v[2]) ^ 2
            + (u[3] * v[1] - u[1] * v[3]) ^ 2
            + (u[1] * v[2] - u[2] * v[1]) ^ 2
    )
    check(
        un > 0 and vn > 0 and un * vn < math.huge and cross > 1e-10 * un * vn,
        "degenerate plane basis"
    )
    check(p.step > 0, "grid step must be positive")
    local fine = p.grid == "both" and p.subdivisions or 1
    check(
        p.grid == "none" or ((p.umax - p.umin) + (p.vmax - p.vmin)) * fine / p.step + 2 <= 2000,
        "plane grid budget exceeded"
    )
    scene.planes[#scene.planes + 1] = p
end
function S.options(o)
    o = o or {}
    local r = {}
    for k, d in pairs({
        azimuth = -45,
        elevation = 25,
        coordinate_scale = 1,
        xmin = -3,
        xmax = 3,
        ymin = -3,
        ymax = 3,
        zmin = -3,
        zmax = 3,
        triad_size = 0.65,
        triad_margin = 0.45,
        tick_step = 1,
        tick_size = 0.04,
        label_distance = 0.16,
        geometry_layer = 20,
        source_layer = 20,
        equipotential_layer = 5,
        field_layer = 10,
        vector_layer = 30,
        axes_layer = 40,
        visibility_tolerance = 1e-8,
        visibility_budget = 1000000,
        fragment_budget = 50000,
        crossing_gap = 0.06,
        mark_position = 0.55,
        mark_spacing = 0,
        mark_clearance = 0.16,
        mark_length = 0.13,
        mark_tolerance = 0.002,
        mark_minimum = 0.25,
    }) do
        r[k] = number(o[k], k, d)
    end
    r.conductors = choice(o.conductors, { line = true, tube = true }, "conductors", "line")
    r.tube_radius =
        choice(o.tube_radius, { physical = true, schematic = true }, "tube radius", "physical")
    r.tube_color = o.tube_color or "orange!70!black"
    r.tube_size = number(o.tube_size, "tube size", 0.04)
    check(r.tube_size > 0, "tube size must be positive")
    r.tube_sides = integer(o.tube_sides, "tube sides", 8, 4, 32)
    r.tube_budget = integer(o.tube_budget, "tube budget", 8000, 100, 20000)
    for _, role in ipairs({
        "geometry",
        "field",
        "vectors",
        "source",
        "equipotential",
        "planes",
        "winding",
        "leads",
        "return",
    }) do
        local key = "show_" .. role
        r[key] = choice(o[key], { ["true"] = true, ["false"] = true }, "show " .. role, "true")
            == "true"
    end
    r.axes =
        choice(o.axes, { none = true, triad = true, origin = true, box = true }, "axes", "triad")
    r.projection = choice(o.projection, { orthographic = true }, "projection", "orthographic")
    r.visibility = choice(o.visibility, { layers = true, vector = true }, "visibility", "layers")
    r.mark_mode =
        choice(o.mark_mode, { bend = true, adaptive = true, tangent = true }, "arrow mode", "bend")
    check(
        r.mark_length > 0 and r.mark_tolerance > 0 and r.mark_minimum > 0 and r.mark_minimum <= 1,
        "invalid arrow fit controls"
    )
    r.hidden = choice(o.hidden, { omit = true, dashed = true }, "hidden lines", "omit")
    r.field_marks =
        choice(o.field_marks, { ["true"] = true, ["false"] = true }, "field marks", "true")
    check(
        r.visibility_tolerance > 0 and r.visibility_tolerance <= 1e-3,
        "invalid visibility tolerance"
    )
    r.visibility_budget = integer(r.visibility_budget, "visibility budget", nil, 100, 10000000)
    r.fragment_budget = integer(r.fragment_budget, "fragment budget", nil, 100, 1000000)
    check(
        r.crossing_gap >= 0
            and r.mark_position >= 0
            and r.mark_position <= 1
            and r.mark_spacing >= 0
            and r.mark_clearance >= 0,
        "invalid visibility layout"
    )
    r.clip = choice(o.clip, { ["true"] = true, ["false"] = true }, "clip", "true") == "true"
    r.ticks = choice(o.ticks, { none = true, major = true }, "ticks", "major")
    r.tick_labels = choice(
        o.tick_labels,
        { ["true"] = true, ["false"] = true },
        "tick labels",
        "true"
    ) == "true"
    r.axis_units = choice(o.axis_units, { ["true"] = true, ["false"] = true }, "axis units", "true")
        == "true"
    r.origin = vector(o.origin or "0/0/0", "axis origin")
    check(
        r.coordinate_scale > 0
            and r.triad_size > 0
            and r.triad_margin >= 0
            and r.tick_step > 0
            and r.tick_size >= 0
            and r.label_distance >= 0,
        "invalid layout size"
    )
    for _, a in ipairs({ "x", "y", "z" }) do
        check(r[a .. "max"] > r[a .. "min"], "reversed view bounds")
        if (r.axes == "origin" or r.axes == "box") and r.ticks == "major" then
            check((r[a .. "max"] - r[a .. "min"]) / r.tick_step <= 200, "axis tick budget exceeded")
        end
        r[a .. "label"] = o[a .. "label"] or ("$" .. a .. "$")
    end
    for _, role in ipairs({ "geometry", "field", "vector", "source", "equipotential", "axes" }) do
        r[role .. "_layer"] = integer(r[role .. "_layer"], role .. " layer", nil, -1000, 1000)
    end
    if r.axes == "origin" then
        for k, a in ipairs({ "x", "y", "z" }) do
            check(
                r.origin[k] >= r[a .. "min"] and r.origin[k] <= r[a .. "max"],
                "axis origin outside view bounds"
            )
        end
    end
    local a, e = (r.azimuth % 360) * math.pi / 180, (r.elevation % 360) * math.pi / 180
    local c, s, ce, se = math.cos(a), math.sin(a), math.cos(e), math.sin(e)
    r.basis = { { c, s, 0 }, { -s * se, c * se, ce }, { s * ce, -c * ce, se } }
    return r
end
--- Project to drawing centimetres and retain view depth for future renderers.
function S.project(p, o)
    local q = {}
    for k, b in ipairs(o.basis) do
        q[k] = b[1] * p[1] + b[2] * p[2] + b[3] * p[3]
    end
    q[1], q[2] = q[1] * o.coordinate_scale, q[2] * o.coordinate_scale
    for k = 1, 3 do
        check(q[k] == q[k] and math.abs(q[k]) < math.huge, "non-finite projected coordinate")
    end
    check(
        math.abs(q[1]) < 400 and math.abs(q[2]) < 400,
        "projected coordinate exceeds TeX range; reduce coordinate scale"
    )
    return q
end
local function corners(o)
    local p = {}
    for _, x in ipairs({ o.xmin, o.xmax }) do
        for _, y in ipairs({ o.ymin, o.ymax }) do
            for _, z in ipairs({ o.zmin, o.zmax }) do
                p[#p + 1] = { x, y, z }
            end
        end
    end
    return p
end
--- Produce projected primitives with a stable explicit layer order.
-- Clip describes the projected bounding rectangle, not a 3D clipping volume.
function S.layout(scene, options)
    local o = S.options(options)
    local out =
        { primitives = {}, options = o, bounds = { math.huge, math.huge, -math.huge, -math.huge } }
    local function emit(p)
        p.sequence = #out.primitives + 1
        out.primitives[#out.primitives + 1] = p
    end
    local function path(points, style, layer, extra)
        local p = { kind = "path", points = {}, style = style, layer = layer }
        for _, v in ipairs(points) do
            p.points[#p.points + 1] = S.project(v, o)
        end
        for k, v in pairs(extra or {}) do
            p[k] = v
        end
        emit(p)
    end
    local box = corners(o)
    for _, p in ipairs(box) do
        local q = S.project(p, o)
        local b = out.bounds
        b[1], b[2], b[3], b[4] =
            math.min(b[1], q[1]), math.min(b[2], q[2]), math.max(b[3], q[1]), math.max(b[4], q[2])
    end
    local styles = {
        geometry = "femm 3d geometry",
        source = "femm 3d source",
        equipotential = "femm 3d equipotential",
        field = "femm 3d field",
        vectors = "femm 3d vector",
    }
    local face_count = 0
    for _, p in ipairs(scene.paths) do
        local shown = o["show_" .. p.role]
            and (p.role ~= "geometry" or o["show_" .. (p.part or "winding")])
        if shown then
            if p.role == "geometry" and o.conductors == "tube" then
                local radius = o.tube_radius == "physical" and p.radius or o.tube_size
                check(
                    radius and radius > 0,
                    "physical tube radius is zero; choose tube radius=schematic"
                )
                local count = 2 * o.tube_sides * (#p.points - 1) + 2
                check(
                    face_count + count <= o.tube_budget,
                    "tube face budget exceeded; reduce geometry samples or tube sides"
                )
                local faces = require("luafemm-tubes").faces(
                    p.points,
                    radius,
                    o.tube_sides,
                    function()
                        face_count = face_count + 1
                    end
                )
                for _, points in ipairs(faces) do
                    local a, b, c =
                        S.project(points[1], o), S.project(points[2], o), S.project(points[3], o)
                    local u, v = {}, {}
                    for k = 1, 3 do
                        u[k] = b[k] - a[k]
                        v[k] = c[k] - a[k]
                    end
                    u[3], v[3] = u[3] * o.coordinate_scale, v[3] * o.coordinate_scale
                    local nx, ny, nz =
                        u[2] * v[3] - u[3] * v[2],
                        u[3] * v[1] - u[1] * v[3],
                        u[1] * v[2] - u[2] * v[1]
                    local length = math.sqrt(nx * nx + ny * ny + nz * nz)
                    local light = length > 0
                            and math.abs((-0.3 * nx + 0.4 * ny + 0.8660254 * nz) / length)
                        or 0
                    local shade = math.floor(55 + 45 * light + 0.5)
                    if nz > o.visibility_tolerance * length then
                        path(
                            points,
                            "fill={" .. o.tube_color .. "!" .. shade .. "!black" .. "},femm 3d tube",
                            o.geometry_layer,
                            {
                                clip = o.clip,
                                role = "surface",
                                occlusion = "opaque",
                                closed = true,
                                tube = true,
                                plane_id = -face_count - #out.primitives - 1,
                            }
                        )
                    end
                end
            else
                path(
                    p.points,
                    styles[p.role],
                    o[(p.role == "vectors" and "vector" or p.role) .. "_layer"],
                    { clip = o.clip, role = p.role }
                )
            end
        end
    end
    for plane_id, p in ipairs(o.show_planes and scene.planes or {}) do
        local function at(a, b)
            return point(p.origin, p.u, p.v, a, b)
        end
        path(
            { at(p.umin, p.vmin), at(p.umax, p.vmin), at(p.umax, p.vmax), at(p.umin, p.vmax) },
            "femm 3d plane," .. p.style,
            p.layer,
            {
                closed = true,
                clip = o.clip,
                role = "surface",
                plane_id = plane_id,
                occlusion = p.occlusion,
                alpha = p.opacity,
            }
        )
        if p.grid ~= "none" then
            local n = p.grid == "both" and p.subdivisions or 1
            local step = p.step / n
            for _, axis in ipairs({ "u", "v" }) do
                for i = math.ceil(p[axis .. "min"] / step), math.floor(p[axis .. "max"] / step) do
                    local t = i * step
                    local pts = axis == "u" and { at(t, p.vmin), at(t, p.vmax) }
                        or { at(p.umin, t), at(p.umax, t) }
                    path(
                        pts,
                        (i % n == 0 and "femm 3d grid" or "femm 3d minor grid")
                            .. ","
                            .. p.grid_style,
                        p.layer,
                        {
                            clip = o.clip,
                            role = "plane",
                            plane_id = plane_id,
                            occlusion = p.occlusion,
                        }
                    )
                end
            end
        end
        if p.label ~= "" then
            emit({
                kind = "label",
                position = S.project(at(p.umax, p.vmax), o),
                text = p.label,
                style = "femm 3d plane label",
                layer = p.layer,
                clip = false,
            })
        end
    end
    if o.axes ~= "none" then
        local origin = o.origin
        if o.axes == "box" then
            -- Put the vertical axis on the left silhouette, not across the object.
            local a = (o.azimuth % 360) * math.pi / 180
            origin = {
                math.cos(a) >= 0 and o.xmin or o.xmax,
                math.sin(a) >= 0 and o.ymin or o.ymax,
                o.zmin,
            }
        end
        if o.axes == "box" then
            for i, p in ipairs(box) do
                for j = i + 1, #box do
                    local different = 0
                    for k = 1, 3 do
                        if p[k] ~= box[j][k] then
                            different = different + 1
                        end
                    end
                    if different == 1 then
                        path({ p, box[j] }, "femm 3d box", o.axes_layer)
                    end
                end
            end
        end
        for k, a in ipairs({ "x", "y", "z" }) do
            local axis_origin = copy(origin)
            if o.axes == "box" and k < 3 then
                local other = k == 1 and 2 or 1
                local key = other == 1 and "x" or "y"
                axis_origin[other] = o.basis[2][other] >= 0 and o[key .. "min"] or o[key .. "max"]
            end
            local first, last = copy(axis_origin), copy(axis_origin)
            first[k], last[k] = o[a .. "min"], o[a .. "max"]
            local q, r = S.project(first, o), S.project(last, o)
            if o.axes == "triad" then
                local unit = { 0, 0, 0 }
                unit[k] = 1
                local v = S.project(unit, o)
                q = { out.bounds[1] + o.triad_size, out.bounds[2] - o.triad_margin - o.triad_size }
                r = {
                    q[1] + o.triad_size * v[1] / o.coordinate_scale,
                    q[2] + o.triad_size * v[2] / o.coordinate_scale,
                }
            end
            local dx, dy = r[1] - q[1], r[2] - q[2]
            local length = math.sqrt(dx * dx + dy * dy)
            -- A view exactly along an axis has no projected axis arrow or ticks.
            if length > 1e-10 then
                dx, dy = dx / length, dy / length
                local nx, ny = dy, -dx
                if o.axes == "box" then
                    local b = out.bounds
                    if nx * (q[1] + r[1] - b[1] - b[3]) + ny * (q[2] + r[2] - b[2] - b[4]) < 0 then
                        nx, ny = -nx, -ny
                    end
                end
                emit({
                    kind = "path",
                    points = { q, r },
                    style = "femm 3d axis,femm 3d " .. a .. " axis",
                    layer = o.axes_layer,
                })
                local label = o[a .. "label"]
                if o.axes ~= "triad" and o.axis_units then
                    label = label .. " (" .. scene.unit .. ")"
                end
                local label_position =
                    { r[1] + o.label_distance * dx, r[2] + o.label_distance * dy }
                local label_anchor = dx > 0.2 and "west"
                    or dx < -0.2 and "east"
                    or dy > 0 and "south"
                    or "north"
                if o.axes == "box" and k < 3 then
                    label_position = {
                        (q[1] + r[1]) / 2 + 3 * o.label_distance * nx,
                        (q[2] + r[2]) / 2 + 3 * o.label_distance * ny,
                    }
                    label_anchor = "center"
                end
                emit({
                    kind = "label",
                    position = label_position,
                    text = label,
                    style = "femm 3d axis label,anchor=" .. label_anchor,
                    layer = o.axes_layer,
                })
                if o.axes ~= "triad" and o.ticks == "major" then
                    for i = math.ceil(o[a .. "min"] / o.tick_step), math.floor(o[a .. "max"] / o.tick_step) do
                        local v = copy(axis_origin)
                        v[k] = i * o.tick_step
                        local t = S.project(v, o)
                        emit({
                            kind = "path",
                            points = {
                                { t[1] - dy * o.tick_size, t[2] + dx * o.tick_size },
                                { t[1] + dy * o.tick_size, t[2] - dx * o.tick_size },
                            },
                            style = "femm 3d tick",
                            layer = o.axes_layer,
                        })
                        if o.tick_labels and math.abs(v[k] - axis_origin[k]) > 1e-10 then
                            emit({
                                kind = "label",
                                position = {
                                    t[1] + nx * o.label_distance,
                                    t[2] + ny * o.label_distance,
                                },
                                text = string.format("\\pgfmathprintnumber{%.10g}", v[k]),
                                style = "femm 3d tick label",
                                layer = o.axes_layer,
                            })
                        end
                    end
                end
            end
        end
    end
    table.sort(out.primitives, function(a, b)
        return a.layer == b.layer and a.sequence < b.sequence or a.layer < b.layer
    end)
    if o.visibility == "vector" then
        return require("luafemm-visibility").render(out)
    end
    return out
end
return S
