-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Passive conversion of static extensions to magnetic FEMM records.
-- No imported expression is evaluated as Lua code. FEMM's diagonal material
-- format cannot represent a rotated tensor with a nonzero off-diagonal entry.
-- @module luafemm-static-exchange
local X = {}
local top = require("luafemm-topology")
local function check(v, s)
    assert(v, "luafemm static exchange: " .. s)
end
local function copy(t)
    if type(t) ~= "table" then
        return t
    end
    local r = {}
    for k, v in pairs(t) do
        r[k] = copy(v)
    end
    return r
end

--- Emit the supported affine-in-polar-angle magnetisation family.
-- FEMM atan2 returns radians; label angles themselves are in degrees.
-- @tparam[opt] table r Source region.
-- @treturn[opt] string FEMM angular expression, or nil for uniform direction.
function X.expression(r)
    local p = r and r.magnetization
    if not p or p.type == "uniform" then
        return nil
    end
    local offset = r.angle + (p.type == "tangential" and 90 or 0)
    local order = p.type == "angular" and p.order or 1
    return string
        .format("%.17g+%.17g*atan2(y-(%.17g),x-(%.17g))*180/pi", offset, order, p.y, p.x)
        :gsub("/pi$", "/3.1415926535897931")
end

--- Recognise our passive grammar plus FEMM theta/radial/tangential presets.
-- Other user expressions remain unsupported; never use load/loadstring here.
-- @tparam[opt] string text Passive imported expression.
-- @tparam number scale File-unit to model-unit coordinate conversion.
-- @return Region angular options and angle offset, or nil for no expression.
function X.parse_expression(text, scale)
    if not text or text == "" then
        return nil
    end
    local s = text:gsub("%s", "")
    s = s:gsub("/3%.1415926535897931$", "/pi")
    local offset, order, y, x = s:match(
        "^([%d.eE+%-]+)%+([%d.eE+%-]+)%*atan2%(y%-%(([%d.eE+%-]+)%),x%-%(([%d.eE+%-]+)%)%)%*180/pi$"
    )
    if s == "theta" then
        offset, order, y, x = "0", "1", "0", "0"
    end
    if s == "theta+90" or s == "90+theta" then
        offset, order, y, x = "90", "1", "0", "0"
    end
    check(
        offset,
        "magnetisation expressions are unsupported outside the passive angular grammar (only passive angular laws are accepted)"
    )
    local values = { offset, order, y, x }
    for i, v in ipairs(values) do
        values[i] = tonumber(v)
        check(
            values[i] and values[i] == values[i] and math.abs(values[i]) < math.huge,
            "invalid angle expression number"
        )
    end
    return {
        magnetization = "angular",
        order = values[2],
        center_y = values[3] * scale,
        center_x = values[4] * scale,
    },
        values[1]
end

--- Relative permeabilities in file axes, rejecting unrepresentable cross terms.
-- @tparam table p Prepared static material.
-- @return File-axis relative permeabilities, or an error for a cross term.
function X.permeability(p)
    check(
        math.abs(p.tensor[2]) <= 1e-12 * math.max(p.tensor[1], p.tensor[3]),
        "rotated anisotropic tensors cannot be represented by FEMM mu_x/mu_y"
    )
    return 1 / p.tensor[1], 1 / p.tensor[3]
end
--- FEMM has no arbitrary sheet-angle field. Refuse a rotated nonlinear law.
function X.sheet_type(p)
    if not p.nonlinear_sheet then
        return 0
    end
    local angle = math.rad(p.material_angle + (p.lamination == "parallel y" and 90 or 0))
    if math.abs(math.sin(angle)) < 1e-12 then
        return 1
    end
    if math.abs(math.cos(angle)) < 1e-12 then
        return 2
    end
    check(false, "rotated nonlinear sheets cannot be represented by FEMM LamType")
end
local function transform(p, v, inverse)
    local c, s = math.cos(math.rad(p.angle)), math.sin(math.rad(p.angle))
    local scale = p.scale or 1
    if inverse then
        local x, y = v[1] - p.cx - p.x, v[2] - p.cy - p.y
        return { p.cx + (c * x + s * y) / scale, p.cy + (-s * x + c * y) / scale }
    end
    local x, y = v[1] - p.cx, v[2] - p.cy
    return { p.cx + scale * (c * x - s * y) + p.x, p.cy + scale * (s * x + c * y) + p.y }
end
local function path_segments(t, name)
    local a = {}
    for _, s in ipairs(t.segments) do
        if s.condition and s.condition.path == name then
            a[#a + 1] = s
        end
    end
    return a
end

--- Split exported periodic paths into pairs of equal straight segments.
-- FEMM attaches one periodic property to exactly two geometric primitives.
-- Native chains remain unrestricted and their original topology is unmodified.
-- @tparam table m Source model with optional periodic pairs.
-- @tparam table t Computational arrangement; never modified in place.
-- @treturn table Arrangement with paired export primitives.
function X.pairs(m, t)
    if not m.pairs or #m.pairs == 0 then
        return t
    end
    t = copy(t)
    local function vertex(p)
        for i, q in ipairs(t.points) do
            if (p[1] - q[1]) ^ 2 + (p[2] - q[2]) ^ 2 < t.tolerance ^ 2 then
                return i
            end
        end
        t.points[#t.points + 1] = p
        return #t.points
    end
    local serial = 0
    for _, p in ipairs(m.pairs) do
        local first, second = path_segments(t, p.first), path_segments(t, p.second)
        local cuts = {}
        for direction, list in ipairs({ first, second }) do
            local other = direction == 1 and second or first
            for _, s in ipairs(list) do
                for _, id in ipairs({ s[1], s[2] }) do
                    local q = transform(p, t.points[id], direction == 2)
                    local best, segment, u = math.huge
                    for _, s2 in ipairs(other) do
                        local d, v = top.distance(q, t.points[s2[1]], t.points[s2[2]])
                        if d < best then
                            best, segment, u = d, s2, v
                        end
                    end
                    check(
                        segment and best <= p.tolerance,
                        "periodic export requires matching traces"
                    )
                    u = math.max(0, math.min(1, u))
                    if u > 1e-8 and u < 1 - 1e-8 then
                        local a, b = t.points[segment[1]], t.points[segment[2]]
                        local list = cuts[segment] or {}
                        cuts[segment] = list
                        list[#list + 1] =
                            { u, vertex({ a[1] + u * (b[1] - a[1]), a[2] + u * (b[2] - a[2]) }) }
                    end
                end
            end
        end
        local split = {}
        for _, s in ipairs(t.segments) do
            local list = cuts[s] or {}
            list[#list + 1] = { 0, s[1] }
            list[#list + 1] = { 1, s[2] }
            table.sort(list, function(a, b)
                return a[1] < b[1]
            end)
            for i = 2, #list do
                if list[i][2] ~= list[i - 1][2] then
                    local own = copy(s)
                    own[1], own[2] = list[i - 1][2], list[i][2]
                    split[#split + 1] = own
                end
            end
        end
        t.segments = split
        first, second = path_segments(t, p.first), path_segments(t, p.second)
        check(#first == #second, "unmatched exported periodic pieces")
        local used = {}
        for _, s in ipairs(first) do
            local a, b = transform(p, t.points[s[1]]), transform(p, t.points[s[2]])
            local match
            for j, q in ipairs(second) do
                if
                    not used[j]
                    and top.distance(a, t.points[q[1]], t.points[q[2]]) <= p.tolerance
                    and top.distance(b, t.points[q[1]], t.points[q[2]]) <= p.tolerance
                then
                    match = j
                    break
                end
            end
            check(match, "unmatched exported periodic segment")
            used[match] = true
            serial = serial + 1
            local q = second[match]
            local aa, bb = t.points[s[1]], t.points[s[2]]
            local cc, dd = t.points[q[1]], t.points[q[2]]
            local length = math.sqrt((aa[1] - bb[1]) ^ 2 + (aa[2] - bb[2]) ^ 2)
            local length2 = math.sqrt((cc[1] - dd[1]) ^ 2 + (cc[2] - dd[2]) ^ 2)
            check(
                math.abs(length - length2) <= 1e-6,
                "FEMM requires paired segment lengths within 1e-6 file units"
            )
            s.condition.fem_pair, q.condition.fem_pair = serial, serial
            s.condition.sign, q.condition.sign = p.sign, p.sign
        end
    end
    return t
end

-- Orientation of an exterior segment with the active domain on its left.
local function oriented(t, s)
    local a, b = t.points[s[1]], t.points[s[2]]
    for _, f in ipairs(t.faces) do
        if not f.hole then
            for ci, poly in ipairs(f.contours) do
                for i, u in ipairs(poly) do
                    local v = poly[i % #poly + 1]
                    if
                        top.distance(a, u, v) <= t.tolerance * 30
                        and top.distance(b, u, v) <= t.tolerance * 30
                    then
                        local same = (b[1] - a[1]) * (v[1] - u[1]) + (b[2] - a[2]) * (v[2] - u[2])
                            > 0
                        if (ci == 1 and same) or (ci > 1 and not same) then
                            return a, b
                        end
                        return b, a
                    end
                end
            end
        end
    end
    error("luafemm static exchange: periodic segment does not bound an active face")
end

--- Infer the FEMM pairing of two exterior straight segments per property.
-- Reverse the second exterior orientation, as FEMM's mesher does. Circular
-- periodic arcs are deliberately excluded until their separate rules are ported.
-- @tparam table m Imported model, receiving sorted pairing declarations.
-- @tparam table t Compiled arrangement with active faces.
-- @tparam table groups Periodic property identifiers mapped to two segments.
function X.import_pairs(m, t, groups)
    m.pairs = {}
    for id, g in pairs(groups) do
        local paths, order = {}, {}
        for _, s in ipairs(g) do
            local name = s.condition.path
            if not paths[name] then
                paths[name] = {}
                order[#order + 1] = name
            end
            paths[name][#paths[name] + 1] = s
        end
        table.sort(order)
        check(#order == 2, "FEMM periodic property must have exactly two open primitives")
        local function chain(list)
            local following, incoming = {}, {}
            for _, s in ipairs(list) do
                local a, b = oriented(t, s)
                check(not following[a] and not incoming[b], "branched periodic primitive")
                following[a], incoming[b] = b, true
            end
            local first
            for a in pairs(following) do
                if not incoming[a] then
                    check(not first, "disconnected periodic primitive")
                    first = a
                end
            end
            check(first, "closed periodic primitive is unsupported")
            local points = { first }
            while following[points[#points]] do
                check(#points <= #list, "cyclic periodic primitive")
                points[#points + 1] = following[points[#points]]
            end
            check(#points == #list + 1, "disconnected periodic primitive")
            return points
        end
        local first, second = chain(paths[order[1]]), chain(paths[order[2]])
        local a, b = first[1], first[#first]
        local d, c = second[1], second[#second]
        local scale = 1
        local axial = m.options.problem_type == "axisymmetric"
        if axial then
            -- Native Kelvin arcs pair equal polar positions, not opposite
            -- tangential orientations. Flux psi=2*pi*r*A_phi is continuous.
            if a[2] > b[2] then
                a, b = b, a
            end
            c, d = second[1], second[#second]
            if c[2] > d[2] then
                c, d = d, c
            end
            check(
                math.abs(a[1]) + math.abs(b[1]) + math.abs(c[1]) + math.abs(d[1])
                    < t.tolerance * 120,
                "Kelvin pairs require meridional semicircles ending on the axis"
            )
            scale = (d[2] - c[2]) / (b[2] - a[2])
            check(scale > 0 and scale < math.huge, "invalid Kelvin interface scale")
        end
        -- Target starts at c and ends at d, opposite its exterior orientation.
        local angle = math.atan(d[2] - c[2], d[1] - c[1]) - math.atan(b[2] - a[2], b[1] - a[1])
        local co, si = math.cos(angle), math.sin(angle)
        m.pairs[#m.pairs + 1] = {
            first = order[1],
            second = order[2],
            angle = math.deg(angle),
            cx = 0,
            cy = 0,
            x = c[1] - co * a[1] + si * a[2],
            y = c[2] - si * a[1] - co * a[2],
            sign = g[1].condition.sign,
            scale = scale,
            kelvin = axial or nil,
            tolerance = math.max(t.tolerance * 30, 1e-8),
        }
        local pairing = m.pairs[#m.pairs]
        pairing.x = c[1] - scale * (co * a[1] - si * a[2])
        pairing.y = c[2] - scale * (si * a[1] + co * a[2])
        for _, q in ipairs(first) do
            local transformed = transform(pairing, q)
            local distance = math.huge
            for j = 2, #second do
                distance = math.min(distance, top.distance(transformed, second[j - 1], second[j]))
            end
            check(
                distance <= pairing.tolerance,
                "periodic primitives require rigidly matching traces"
            )
        end
    end
    table.sort(m.pairs, function(a, b)
        return a.first < b.first
    end)
end
return X
