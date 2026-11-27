-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Physical 3D integration followed by an orthographic camera projection.
-- Clipping boxes stop traces in physical space. A camera never changes a source,
-- field sample or streamline. No hidden-surface or perspective claim is made.
local V = {}
local G = require("luafemm-biot-geometry")
local B = require("luafemm-biot")
local check, n = G.check, G.number
local add, sub, mul, norm = G.add, G.sub, G.mul, G.norm
function V.options(o)
    o = o or {}
    local r = {}
    for k, d in pairs({
        azimuth = -45,
        elevation = 25,
        step = 0.08,
        tolerance = 0.0002,
        max_length = 30,
        max_steps = 3000,
        null_field = 1e-14,
        geometry_samples = 64,
        xmin = -3,
        xmax = 3,
        ymin = -3,
        ymax = 3,
        zmin = -3,
        zmax = 3,
        spacing = 0.5,
        vector_length = 0.25,
        offset = 0,
    }) do
        r[k] = n(o[k], k, d)
    end
    check(
        r.step > 0 and r.tolerance > 0 and r.max_length > 0 and r.null_field >= 0,
        "invalid tracing controls"
    )
    check(
        r.max_steps % 1 == 0 and r.max_steps >= 1 and r.max_steps <= 100000,
        "invalid maximum steps"
    )
    check(
        r.geometry_samples % 1 == 0 and r.geometry_samples >= 8 and r.geometry_samples <= 1000,
        "geometry samples must be 8 to 1000"
    )
    check(r.spacing > 0 and r.vector_length > 0, "invalid vector size or spacing")
    for _, a in ipairs({ "x", "y", "z" }) do
        check(r[a .. "max"] > r[a .. "min"], "reversed view bounds")
    end
    r.plane = o.plane or "xz"
    check(r.plane == "xy" or r.plane == "xz" or r.plane == "yz", "unknown sampling plane")
    r.direction = o.direction or "both"
    check(
        r.direction == "both" or r.direction == "forward" or r.direction == "backward",
        "unknown tracing direction"
    )
    return r
end
function V.project(p, o)
    local a, e = o.azimuth * math.pi / 180, o.elevation * math.pi / 180
    return {
        math.cos(a) * p[1] + math.sin(a) * p[2],
        -math.sin(a) * math.sin(e) * p[1] + math.cos(a) * math.sin(e) * p[2] + math.cos(e) * p[3],
    }
end
function V.lines(m, seeds, o)
    o = V.options(o)
    return require("luafemm-stream3d").lines(function(p)
        local f = B.sample(m, p)
        return { valid = f.valid, vector = f.vector, magnitude = f.B }
    end, seeds, o)
end
function V.geometry(m, o)
    o = V.options(o)
    local paths = {}
    for source_id, curves in ipairs(m.compiled) do
        for _, c in ipairs(curves) do
            local points = {}
            if c.infinite then
                local u = mul(sub(c.b, c.a), 1 / c.speed)
                local lo, hi = -math.huge, math.huge
                for k, a in ipairs({ "x", "y", "z" }) do
                    if math.abs(u[k]) < 1e-14 then
                        if c.a[k] < o[a .. "min"] or c.a[k] > o[a .. "max"] then
                            lo, hi = 1, 0
                        end
                    else
                        local x, y =
                            (o[a .. "min"] - c.a[k]) / u[k], (o[a .. "max"] - c.a[k]) / u[k]
                        lo, hi = math.max(lo, math.min(x, y)), math.min(hi, math.max(x, y))
                    end
                end
                if hi > lo then
                    points = { add(c.a, mul(u, lo)), add(c.a, mul(u, hi)) }
                end
            else
                local samples = c.line and 1 or math.ceil(o.geometry_samples * c.panels / 16)
                for i = 0, samples do
                    points[#points + 1] = c.evaluate(i / samples)
                end
            end
            if #points >= 2 then
                paths[#paths + 1] = {
                    points = points,
                    part = c.part or "winding",
                    radius = m.sources[source_id].wire_radius,
                    source = source_id,
                }
            end
        end
    end
    return paths
end
function V.vectors(m, o)
    o = V.options(o)
    local paths = {}
    local a, b, k
    if o.plane == "xy" then
        a, b, k = 1, 2, 3
    elseif o.plane == "xz" then
        a, b, k = 1, 3, 2
    else
        a, b, k = 2, 3, 1
    end
    local axis = { "x", "y", "z" }
    local na = math.floor((o[axis[a] .. "max"] - o[axis[a] .. "min"]) / o.spacing)
    local nb = math.floor((o[axis[b] .. "max"] - o[axis[b] .. "min"]) / o.spacing)
    check((na + 1) * (nb + 1) <= 10000, "vector sample budget exceeded")
    for i = 0, na do
        for j = 0, nb do
            local p = { 0, 0, 0 }
            p[a] = o[axis[a] .. "min"] + i * o.spacing
            p[b] = o[axis[b] .. "min"] + j * o.spacing
            p[k] = o.offset
            local f = B.sample(m, p)
            if f.valid and f.B > o.null_field then
                paths[#paths + 1] = { points = { p, add(p, mul(f.vector, o.vector_length / f.B)) } }
            end
        end
    end
    return paths
end
return V
