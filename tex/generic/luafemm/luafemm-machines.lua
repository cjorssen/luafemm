-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Planar slotted machines. Geometry and winding compilation are independent of TeX.
-- Lengths are model millimetres; angles are mechanical degrees counterclockwise.
-- Slots are numbered in declaration order, never renumbered by geometric sorting.
-- @module luafemm-machines
local M = {}
local pi, sin, cos, sqrt = math.pi, math.sin, math.cos, math.sqrt
local function check(ok, message)
    assert(ok, "luafemm machine: " .. message)
end
local function number(x, name)
    local n = tonumber(x)
    check(n and n == n and math.abs(n) < math.huge, "finite number required for " .. name)
    return n
end
local function clone(value)
    if type(value) ~= "table" then
        return value
    end
    local result = {}
    for k, v in pairs(value) do
        result[k] = clone(v)
    end
    return result
end
local function options(defaults, input)
    local out = {}
    for k in pairs(input or {}) do
        check(defaults[k] ~= nil, "unknown option " .. k)
    end
    for k, v in pairs(defaults) do
        local x = input and input[k]
        if x == nil or x == "" then
            x = v
        end
        out[k] = type(v) == "number" and number(x, k) or clone(x)
    end
    return out
end
local function list(v)
    if type(v) == "table" then
        return v
    end
    local t = {}
    for x in tostring(v or ""):gmatch("[^,%s]+") do
        t[#t + 1] = number(x, "list entry")
    end
    return t
end
local function polar(r, a)
    a = math.rad(a)
    return { r * cos(a), r * sin(a) }
end
local function append(t, p)
    t[#t + 1] = p
end
local function arc(t, r, a, b, tolerance, include_start)
    local step = math.deg(2 * math.acos(math.max(-1, 1 - tolerance / r)))
    local n = math.max(1, math.ceil(math.abs(b - a) / math.min(10, step)))
    check(n <= 10000, "curve tolerance requests too many arc vertices")
    for i = include_start and 0 or 1, n do
        append(t, polar(r, a + (b - a) * i / n))
    end
end
local function circle(r, tol)
    local t = {}
    arc(t, r, 0, 360, tol, false)
    return t
end
local function square(r)
    return { { -r, -r }, { r, -r }, { r, r }, { -r, r } }
end
local function copy(t)
    local out = {}
    for k, v in pairs(t or {}) do
        out[k] = v
    end
    return out
end
-- Round only internal notch corners. Both core and slot reuse these vertices.
local function rounded(path, radius, tolerance)
    if radius == 0 then
        return path
    end
    local out = { path[1] }
    for i = 2, #path - 1 do
        local a, p, b = path[i - 1], path[i], path[i + 1]
        local ux, uy = a[1] - p[1], a[2] - p[2]
        local vx, vy = b[1] - p[1], b[2] - p[2]
        local l, h = sqrt(ux * ux + uy * uy), sqrt(vx * vx + vy * vy)
        ux, uy, vx, vy = ux / l, uy / l, vx / h, vy / h
        local dot = math.max(-1, math.min(1, ux * vx + uy * vy))
        local half = math.acos(dot) / 2
        if math.abs(pi / 2 - half) < 1e-8 then
            append(out, p)
        else
            local distance = radius / math.tan(half)
            check(distance < 0.49 * math.min(l, h), "corner radius exceeds adjacent notch edges")
            local cx, cy =
                p[1] + (ux + vx) * radius / (2 * sin(half) * cos(half)),
                p[2] + (uy + vy) * radius / (2 * sin(half) * cos(half))
            local start = math.atan(p[2] + uy * distance - cy, p[1] + ux * distance - cx)
            local stop = math.atan(p[2] + vy * distance - cy, p[1] + vx * distance - cx)
            local sweep = (stop - start + pi) % (2 * pi) - pi
            local step = 2 * math.acos(math.max(-1, 1 - tolerance / radius))
            local n = math.max(1, math.ceil(math.abs(sweep) / math.min(pi / 12, step)))
            check(n < 10000, "corner tolerance requests too many vertices")
            for j = 0, n do
                append(out, {
                    cx + radius * cos(start + sweep * j / n),
                    cy + radius * sin(start + sweep * j / n),
                })
            end
        end
    end
    append(out, path[#path])
    return out
end
local defaults = {
    radius = 29,
    bore_radius = 30,
    outer_radius = 45,
    outer_width = 90,
    outer_shape = "circle",
    shaft_radius = 0,
    housing_thickness = 0,
    minimum_yoke = 1,
    core_material = "pure-iron",
    shaft_material = "air",
    housing_material = "air",
    mesh_size = 0,
    slot_mesh_size = 0,
    conductor_mesh_size = 0,
    curve_tolerance = 0.03,
    magnetization_angle = 0,
    slots = {},
    winding = {},
    overrides = {},
    coils = {},
    conductors = {},
}
local slot_defaults = {
    placement = "uniform",
    count = 0,
    start_angle = 0,
    pitch = 0,
    angles = "",
    shape = "trapezoid",
    depth = 8,
    neck_depth = 1,
    opening_width = 2,
    width = 4,
    bottom_width = 5,
    opening_angle = 0,
    offset = 0,
    corner_radius = 0,
    material = "air",
    conductor_material = "copper",
    clearance = 0.4,
    layers = 1,
    layer_gap = 0.3,
    current_density = 0,
    mesh_size = 0,
    conductor_mesh_size = 0,
    style = "",
    conductor_style = "",
}
--- Compile a rotor or stator into compound core and simple conductor regions.
-- Input slots may have uniform or explicit positions. overrides is an array of
-- indexed slot records; coils is an array of balanced slot-pair excitations.
-- @tparam string kind rotor or stator.
-- @tparam[opt] table input Geometry, slots, overrides, coils and winding records.
-- @treturn table Component compatible with the generic PGF component emitter.
function M.make(kind, input)
    check(kind == "rotor" or kind == "stator", "unknown machine kind")
    local o = options(defaults, input)
    local s = options(slot_defaults, o.slots)
    check(o.outer_shape == "circle" or o.outer_shape == "square", "invalid outer shape")
    check(
        o.curve_tolerance > 0 and o.minimum_yoke > 0,
        "positive tolerance and minimum yoke required"
    )
    for _, k in ipairs({
        "mesh_size",
        "slot_mesh_size",
        "conductor_mesh_size",
        "shaft_radius",
        "housing_thickness",
    }) do
        check(o[k] >= 0, k .. " must be nonnegative")
    end
    local rotor = kind == "rotor"
    local R = rotor and o.radius or o.bore_radius
    local outer = o.outer_shape == "circle" and o.outer_radius or o.outer_width / 2
    check(R > 0 and (rotor or outer > R + o.minimum_yoke), "invalid air-gap or outer radius")
    check(not rotor or o.shaft_radius + o.minimum_yoke < R, "shaft consumes rotor")
    check(s.placement == "uniform" or s.placement == "explicit", "invalid slot placement")
    local angles = s.placement == "explicit" and list(s.angles) or {}
    if s.placement == "uniform" then
        check(s.count >= 0 and s.count % 1 == 0 and s.count <= 512, "slot count must be 0..512")
        for i = 1, s.count do
            angles[i] = s.start_angle + (i - 1) * (s.pitch == 0 and 360 / s.count or s.pitch)
        end
    end
    check(#angles <= 512, "at most 512 slots supported")
    local overrides = {}
    for _, record in ipairs(o.overrides) do
        local v = copy(record)
        local i = number(v.index, "slot index")
        v.index = nil
        check(
            i % 1 == 0 and i >= 1 and i <= #angles and not overrides[i],
            "invalid or duplicate slot override"
        )
        for _, k in ipairs({ "placement", "count", "start_angle", "pitch", "angles" }) do
            check(
                v[k] == nil or v[k] == "",
                "array placement belongs in slots, not a slot override"
            )
        end
        overrides[i] = v
    end
    local slots, sorted, anchors, conductors =
        {}, {}, { center = { 0, 0 }, ["core center"] = { 0, 0 }, ["window center"] = { 0, 0 } }, {}
    local sign = rotor and -1 or 1
    for i, a in ipairs(angles) do
        local values = copy(s)
        for k, v in pairs(overrides[i] or {}) do
            if v ~= "" then
                values[k] = v
            end
        end
        local v = options(slot_defaults, values)
        local theta = (number(a, "slot angle") + v.offset) % 360
        local mouth = v.opening_angle > 0 and 2 * R * sin(math.rad(v.opening_angle / 2))
            or v.opening_width
        check(v.opening_angle >= 0 and v.opening_angle < 90, "opening angle must be in [0,90)")
        check(
            v.shape == "rectangle" or v.shape == "trapezoid" or v.shape == "radial",
            "invalid slot shape"
        )
        check(
            v.depth > v.neck_depth and v.neck_depth > 0 and mouth > 0 and mouth < R,
            "invalid slot opening/depth"
        )
        check(v.width >= mouth and v.bottom_width > 0, "slot body must contain its opening")
        local w0 = v.width
        local w1 = v.shape == "rectangle" and w0 or v.bottom_width
        if v.shape == "radial" then
            w1 = w0 * (R + sign * v.depth) / (R + sign * v.neck_depth)
        end
        check(
            w1 > 0 and R + sign * v.depth > (rotor and o.shaft_radius + o.minimum_yoke or 0),
            "slot reaches shaft or machine centre"
        )
        check(
            v.corner_radius >= 0 and v.clearance >= v.corner_radius,
            "corner radius must not exceed conductor clearance"
        )
        local a0 = math.deg(math.asin(mouth / (2 * R)))
        local ca, sa = cos(math.rad(theta)), sin(math.rad(theta))
        local function point(d, x)
            local r = R + sign * d
            return { r * ca - x * sa, r * sa + x * ca }
        end
        local edge = sign * (sqrt(R * R - mouth * mouth / 4) - R)
        check(not rotor or v.neck_depth > edge, "rotor neck must lie inside its nominal circle")
        local notch = { point(edge, -mouth / 2), point(v.neck_depth, -mouth / 2) }
        if w0 > mouth then
            append(notch, point(v.neck_depth, -w0 / 2))
        end
        append(notch, point(v.depth, -w1 / 2))
        append(notch, point(v.depth, w1 / 2))
        if w0 > mouth then
            append(notch, point(v.neck_depth, w0 / 2))
        end
        append(notch, point(v.neck_depth, mouth / 2))
        append(notch, point(edge, mouth / 2))
        notch = rounded(notch, v.corner_radius, o.curve_tolerance)
        local span = a0
        for _, p in ipairs(notch) do
            local r = sqrt(p[1] ^ 2 + p[2] ^ 2)
            if rotor then
                check(r >= o.shaft_radius + o.minimum_yoke, "slot violates shaft clearance")
                check(r <= R + 1e-10, "rotor slot body extends beyond nominal surface")
            elseif o.outer_shape == "circle" then
                check(r <= outer - o.minimum_yoke, "slot violates minimum yoke")
            else
                check(
                    math.max(math.abs(p[1]), math.abs(p[2])) <= outer - o.minimum_yoke,
                    "slot violates square yoke"
                )
            end
            local x, y = p[1] * ca + p[2] * sa, -p[1] * sa + p[2] * ca
            span = math.max(span, math.abs(math.deg(math.atan(y, x))))
        end
        local air = copy(notch)
        arc(air, R, theta + a0, theta - a0, o.curve_tolerance, false)
        table.remove(air) -- the final point is already the first notch vertex
        slots[i] = {
            index = i,
            theta = theta,
            half = span,
            opening = a0,
            notch = notch,
            air = air,
            options = v,
        }
        sorted[#sorted + 1] = slots[i]
        anchors["slot " .. i] = point((v.depth + v.neck_depth) / 2, 0)
        anchors["slot " .. i .. " opening"] = polar(R, theta)
        anchors["slot " .. i .. " bottom"] = point(v.depth, 0)
        check(v.layers >= 1 and v.layers % 1 == 0 and v.layers <= 8, "layers must be 1..8")
        check(
            v.clearance > 0 and v.layer_gap >= 0,
            "positive conductor clearance and nonnegative layer gap required"
        )
        check(v.mesh_size >= 0 and v.conductor_mesh_size >= 0, "negative slot mesh size")
        local lo, hi = v.neck_depth + v.clearance, v.depth - v.clearance
        local thickness = (hi - lo - (v.layers - 1) * v.layer_gap) / v.layers
        check(thickness > 0, "no room for conductor layers")
        local slope = (w1 - w0) / (2 * (v.depth - v.neck_depth))
        local function half(d)
            return (w0 / 2 + slope * (d - v.neck_depth)) - v.clearance * sqrt(1 + slope * slope)
        end
        conductors[i] = {}
        for layer = 1, v.layers do
            local d0 = lo + (layer - 1) * (thickness + v.layer_gap)
            local d1 = d0 + thickness
            check(math.min(half(d0), half(d1)) > 0, "clearance consumes conductor")
            local p = {
                point(d0, -half(d0)),
                point(d1, -half(d1)),
                point(d1, half(d1)),
                point(d0, half(d0)),
            }
            conductors[i][layer] = {
                role = "conductor",
                material = v.conductor_material,
                points = p,
                style = v.conductor_style,
                current_density = v.current_density,
                mesh_size = v.conductor_mesh_size > 0 and v.conductor_mesh_size
                    or o.conductor_mesh_size,
            }
            anchors["slot " .. i .. " conductor " .. layer] = point((d0 + d1) / 2, 0)
        end
    end
    table.sort(sorted, function(a, b)
        return a.theta < b.theta
    end)
    for i, a in ipairs(sorted) do
        local b = sorted[i % #sorted + 1]
        check(
            (b.theta - a.theta) % 360 > a.half + b.half + 1e-8 or #sorted == 1,
            "overlapping slot angular envelopes"
        )
    end
    local boundary = {}
    if #sorted == 0 then
        boundary = circle(R, o.curve_tolerance)
    else
        for i, a in ipairs(sorted) do
            for _, p in ipairs(a.notch) do
                append(boundary, p)
            end
            local b = sorted[i % #sorted + 1]
            local finish = b.theta - b.opening + (i == #sorted and 360 or 0)
            arc(boundary, R, a.theta + a.opening, finish, o.curve_tolerance, false)
            table.remove(boundary) -- next notch supplies the shared endpoint
        end
    end
    local regions = {}
    local function add(role, material, contours, spacing)
        regions[#regions + 1] = {
            role = role,
            material = material,
            contours = contours,
            mesh_size = spacing or o.mesh_size,
        }
    end
    if rotor then
        local contours = { boundary }
        if o.shaft_radius > 0 then
            append(contours, circle(o.shaft_radius, o.curve_tolerance))
        end
        add("core", o.core_material, contours)
        if o.shaft_radius > 0 then
            add("shaft", o.shaft_material, { circle(o.shaft_radius, o.curve_tolerance) })
        end
    else
        local outside = o.outer_shape == "circle" and circle(outer, o.curve_tolerance)
            or square(outer)
        add("core", o.core_material, { outside, boundary })
        if o.housing_thickness > 0 then
            local h = outer + o.housing_thickness
            add(
                "housing",
                o.housing_material,
                { o.outer_shape == "circle" and circle(h, o.curve_tolerance) or square(h), outside }
            )
        end
    end
    for _, record in ipairs(o.conductors) do
        local v = options(
            { slot = 1, layer = 1, material = "", current_density = "", style = "" },
            record
        )
        local r = conductors[v.slot] and conductors[v.slot][v.layer]
        check(r, "invalid conductor slot/layer")
        if v.material ~= "" then
            r.material = v.material
        end
        if v.style ~= "" then
            r.style = v.style
        end
        if v.current_density ~= "" then
            r.current_density = number(v.current_density, "conductor current density")
        end
    end
    -- Coil records accumulate only on explicitly assigned conductor layers.
    -- A slot's raw current density and winding ampere-turns are exclusive.
    local function coil(record)
        local v = options(
            { slots = "", layer = 1, turns = 1, current = 0, ampere_turns = "", material = "" },
            record
        )
        local pair = list(v.slots)
        check(#pair == 2 and pair[1] ~= pair[2], "coil requires two distinct slots")
        local ni = v.ampere_turns ~= "" and number(v.ampere_turns, "ampere turns")
            or v.turns * v.current
        check(
            v.ampere_turns == ""
                or (record.turns == nil and record.current == nil)
                or (record.turns == "" and record.current == ""),
            "choose ampere turns or turns/current"
        )
        check(v.turns >= 0, "negative turn count; reverse current instead")
        for j, index in ipairs(pair) do
            local r = conductors[index] and conductors[index][v.layer]
            check(r, "invalid coil slot/layer")
            check(r.current_density == 0, "coil conflicts with current density")
            r.turns = (r.turns or 0) + (j == 1 and ni or -ni)
            if v.material ~= "" then
                r.material = v.material
            end
        end
    end
    for _, v in ipairs(o.coils) do
        coil(v)
    end
    local w =
        options({ distribution = "none", ampere_turns = 1000, pole_pairs = 1, axis = 0 }, o.winding)
    check(
        w.distribution == "none" or w.distribution == "uniform" or w.distribution == "sine",
        "unknown winding distribution"
    )
    if w.distribution ~= "none" then
        check(#o.coils == 0, "automatic and explicit winding cannot be combined")
        check(w.pole_pairs >= 1 and w.pole_pairs % 1 == 0, "positive integer pole pairs required")
        local weights, total, sum = {}, 0, 0
        for i, v in ipairs(slots) do
            local z = sin(math.rad(w.pole_pairs * (v.theta - w.axis)))
            z = math.abs(z) < 1e-12 and 0
                or (w.distribution == "uniform" and (z > 0 and 1 or -1) or z)
            weights[i] = z
            total = total + math.max(z, 0)
            sum = sum + z
        end
        check(
            total > 0 and math.abs(sum) < 1e-9 * total,
            "automatic winding requires balanced slot positions"
        )
        for i, z in ipairs(weights) do
            for _, r in ipairs(conductors[i]) do
                check(r.current_density == 0, "automatic winding conflicts with current density")
                r.turns = w.ampere_turns * z / total / #conductors[i]
            end
        end
    end
    for i, v in ipairs(slots) do
        add(
            "slot",
            v.options.material,
            { v.air },
            v.options.mesh_size > 0 and v.options.mesh_size or o.slot_mesh_size
        )
        regions[#regions].style = v.options.style
        for _, r in ipairs(conductors[i]) do
            append(regions, r)
        end
    end
    local extent = rotor and R or outer + o.housing_thickness
    for name, p in pairs({
        north = { 0, extent },
        south = { 0, -extent },
        east = { extent, 0 },
        west = { -extent, 0 },
        ["north east"] = { extent, extent },
        ["north west"] = { -extent, extent },
        ["south east"] = { extent, -extent },
        ["south west"] = { -extent, -extent },
        ["gap east"] = { R, 0 },
        ["gap north"] = { 0, R },
        ["gap west"] = { -R, 0 },
        ["gap south"] = { 0, -R },
    }) do
        anchors[name] = p
    end
    return {
        kind = kind,
        options = o,
        regions = regions,
        anchors = anchors,
        cycles = {},
        half_width = extent,
        half_height = extent,
        airgap_radius = R,
        machine = true,
    }
end
return M
