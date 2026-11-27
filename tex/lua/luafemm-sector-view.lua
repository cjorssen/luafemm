-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Physical-frame continuation and clipping of concentric sector solutions.
-- Computation stays in body frames. This module never changes the mesh or its
-- solution: scalar A changes sign on antiperiodic copies, vectors also rotate.
local V = {}
local pi = math.pi
local function rotate(x, y, cx, cy, a)
    local c, s = math.cos(a), math.sin(a)
    x, y = x - cx, y - cy
    return cx + c * x - s * y, cy + s * x + c * y
end
--- Find a bulk element from a physical SI point, returning its reference point,
-- body rotation (radians) and continuation sign. Gap points use the gap operator.
function V.locate(m, x, y)
    local o, u = m.sliding.options, m.options.unit
    local cx, cy = o.center_x * u, o.center_y * u
    local rotor = (x - cx) ^ 2 + (y - cy) ^ 2 < (o.inner_radius * u) ^ 2
    local span, start = 2 * pi / o.sectors, math.rad(o.start_angle or 0)
    local angle = rotor and math.rad(m.sliding.angle) or 0
    local theta = math.atan(y - cy, x - cx)
    local sector = math.floor((theta - start - angle) / span)
    local turn = angle + sector * span
    local xx, yy = rotate(x, y, cx, cy, -turn)
    local sign = o.periodicity == "antiperiodic" and (-1) ^ sector or 1
    for i = 1, m.sliding.bulk_count do
        local e = m.triangles[i]
        local p = m.nodes[e.ids[1]]
        local b = e.gx[2] * (xx - p[1]) + e.gy[2] * (yy - p[2])
        local c = e.gx[3] * (xx - p[1]) + e.gy[3] * (yy - p[2])
        if b >= -1e-9 and c >= -1e-9 and b + c <= 1 + 1e-9 then
            return e, xx, yy, turn, sign
        end
    end
end
function V.vector(x, y, angle, sign)
    local c, s = math.cos(angle), math.sin(angle)
    return sign * (c * x - s * y), sign * (s * x + c * y)
end
-- Convex wedge clipping, interpolating A at intersections. Sectors have span
-- at most pi, so the two oriented half-planes define the physical display.
local function clip(poly, cx, cy, start, span)
    for _, side in ipairs({ { start, 1 }, { start + span, -1 } }) do
        local c, s = math.cos(side[1]), math.sin(side[1])
        local function distance(p)
            local v = side[2] * (c * (p[2] - cy) - s * (p[1] - cx))
            -- Trigonometric roundoff must not erase a segment lying on a
            -- radial cut (cos(pi/2) is not represented as exact zero).
            local scale = math.abs(p[1] - cx) + math.abs(p[2] - cy)
            return math.abs(v) <= 1e-13 * scale and 0 or v
        end
        local out = {}
        for i, p in ipairs(poly) do
            local q = poly[i % #poly + 1]
            local a, b = distance(p), distance(q)
            if a >= 0 then
                out[#out + 1] = p
            end
            if (a >= 0) ~= (b >= 0) then
                local t = a / (a - b)
                out[#out + 1] =
                    { p[1] + t * (q[1] - p[1]), p[2] + t * (q[2] - p[2]), p[3] + t * (q[3] - p[3]) }
            end
        end
        poly = out
    end
    return poly
end
--- Build a rendering-only P1 view. Bulk interpolation is exact; the analytical
-- gap is sampled at a resolution controlled by contour subdivisions. Traces of
-- separate operators deliberately keep separate vertices (no false field joins).
function V.mesh(m)
    local o, u = m.sliding.options, m.options.unit
    local cx, cy = o.center_x * u, o.center_y * u
    local start, span = math.rad(o.start_angle or 0), 2 * pi / o.sectors
    local view = { options = m.options, stats = m.stats, nodes = {}, triangles = {}, A = {} }
    local map = {}
    local scale = o.outer_radius * u
    local function triangle(poly, part)
        local ids = {}
        for _, p in ipairs(poly) do
            local key = part
                .. ":"
                .. string.format("%.11f:%.11f", (p[1] - cx) / scale, (p[2] - cy) / scale)
            local id = map[key]
            if not id then
                id = #view.nodes + 1
                map[key] = id
                view.nodes[id] = { p[1], p[2] }
                view.A[id] = p[3]
            end
            ids[#ids + 1] = id
        end
        for i = 2, #ids - 1 do
            local a, b, c = ids[1], ids[i], ids[i + 1]
            local p, q, r = view.nodes[a], view.nodes[b], view.nodes[c]
            local twice = (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])
            if twice > 1e-22 * scale ^ 2 then
                local aa, bb, cc = view.A[a], view.A[b], view.A[c]
                view.triangles[#view.triangles + 1] = {
                    ids = { a, b, c },
                    bx = (aa * (r[1] - q[1]) + bb * (p[1] - r[1]) + cc * (q[1] - p[1])) / twice,
                    by = -(aa * (q[2] - r[2]) + bb * (r[2] - p[2]) + cc * (p[2] - q[2])) / twice,
                }
            end
        end
    end
    for i = 1, m.sliding.bulk_count do
        local e = m.triangles[i]
        local x, y = 0, 0
        for _, id in ipairs(e.ids) do
            x = x + m.nodes[id][1] / 3
            y = y + m.nodes[id][2] / 3
        end
        local rotor = (x - cx) ^ 2 + (y - cy) ^ 2 < (o.inner_radius * u) ^ 2
        for k = 0, (rotor and o.sectors - 1 or 0) do
            local a = rotor and math.rad(m.sliding.angle) + k * span or 0
            local sign = o.periodicity == "antiperiodic" and (-1) ^ k or 1
            local poly = {}
            for _, id in ipairs(e.ids) do
                local p = m.nodes[id]
                local xx, yy = rotate(p[1], p[2], cx, cy, a)
                poly[#poly + 1] = { xx, yy, sign * m.A[id] }
            end
            triangle(clip(poly, cx, cy, start, span), rotor and "rotor" or "stator")
        end
    end
    local G = require("luafemm-airgap")
    local angular = math.max(16, math.ceil(4 * G.prepare(m).modes / o.sectors))
        * (m.options.contour_subdivisions or 1)
    local radial = math.max(2, math.ceil((o.outer_radius - o.inner_radius) / m.options.h))
        * (m.options.contour_subdivisions or 1)
    local function point(i, j)
        local r = (o.inner_radius + (o.outer_radius - o.inner_radius) * j / radial) * u
        local a = start + span * i / angular
        local x, y = cx + r * math.cos(a), cy + r * math.sin(a)
        local _, _, v = G.field(m, m.A, x, y)
        return { x, y, v }
    end
    for i = 0, angular - 1 do
        for j = 0, radial - 1 do
            local a, b, c, d = point(i, j), point(i + 1, j), point(i + 1, j + 1), point(i, j + 1)
            triangle({ a, d, c }, "gap")
            triangle({ a, c, b }, "gap")
        end
    end
    return view
end
--- Transform and clip interface segments for the geometry pic (model units).
function V.geometry(m, t)
    local o = m.sliding.options
    local start, span = math.rad(o.start_angle or 0), 2 * pi / o.sectors
    local out = {}
    for _, s in ipairs(t.segments) do
        local p, q = t.points[s[1]], t.points[s[2]]
        local x, y = (p[1] + q[1]) / 2 - o.center_x, (p[2] + q[2]) / 2 - o.center_y
        local rotor = x * x + y * y < o.inner_radius ^ 2 * (1 + 1e-9)
        if s.condition and s.condition.type == "periodic" then
            rotor = false
        end
        for k = 0, (rotor and o.sectors - 1 or 0) do
            local a = rotor and math.rad(m.sliding.angle) + k * span or 0
            local px, py = rotate(p[1], p[2], o.center_x, o.center_y, a)
            local qx, qy = rotate(q[1], q[2], o.center_x, o.center_y, a)
            -- A segment is a two-point degenerate polygon; clipping retains
            -- intersections twice. Its extreme points define the visible line.
            local poly = clip({ { px, py, 0 }, { qx, qy, 0 } }, o.center_x, o.center_y, start, span)
            local first = poly[1]
            if first then
                local last, best = first, 0
                for _, v in ipairs(poly) do
                    local d = (v[1] - first[1]) ^ 2 + (v[2] - first[2]) ^ 2
                    if d > best then
                        last, best = v, d
                    end
                end
                if best > 1e-20 then
                    out[#out + 1] = { first, last, group = s.group }
                end
            end
        end
    end
    return out
end
return V
