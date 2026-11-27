-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Electric capture adapters for provider-independent spatial scenes.
-- Slice mode embeds an existing section without changing its physics. Spatial
-- mode evaluates Coulomb fields or the exact symmetry lift of an existing model.
-- Neither mode introduces a three-dimensional finite-element discretization.
local C = {}
local E = require("luafemm-electrostatic")
local R = require("luafemm-electric-view")
local S = require("luafemm-scene3d")
local T = require("luafemm-stream3d")
local check, n = S.check, S.number
local sqrt, pi = math.sqrt, math.pi
local function norm(p)
    return sqrt(p[1] ^ 2 + p[2] ^ 2 + p[3] ^ 2)
end
local function lift(p, plane, offset)
    if plane == "xy" then
        return { p[1], p[2], offset }
    elseif plane == "xz" then
        return { p[1], offset, p[2] }
    else
        return { offset, p[1], p[2] }
    end
end
--- Samples are SI; p and exclusion radii are expressed in model length units.
-- Point-source exclusions are display masks, not finite charged bodies.
function C.sample(m, p)
    local x, y, z = n(p[1], "x"), n(p[2], "y"), n(p[3], "z")
    E.freeze(m)
    local f, v
    if m.fe and m.fe.axisymmetric then
        local r = sqrt(x * x + y * y)
        f = E.sample(m, r, z)
        if f then
            v = { r > 0 and f.Ex * x / r or 0, r > 0 and f.Ex * y / r or 0, f.Ey }
        end
    elseif m.fe or (m.classic and (m.classic.planar or m.classic.cylindrical)) then
        f = E.sample(m, x, y)
        if f then
            v = { f.Ex, f.Ey, 0 }
        end
    elseif m.classic then
        local c = m.classic
        local d = { x - c.x, y - c.y, z }
        local r = norm(d)
        f = E.sample(m, c.x + r, c.y)
        if f then
            v = r > 0 and { f.Ex * d[1] / r, f.Ex * d[2] / r, f.Ex * d[3] / r } or { 0, 0, 0 }
        end
    else
        local potential = 0
        v = { 0, 0, 0 }
        local k = 1 / (4 * pi * m.epsilon)
        for _, q in ipairs(m.sources) do
            local d = { (x - q.x) * m.unit, (y - q.y) * m.unit, z * m.unit }
            local r = norm(d)
            if r == 0 then
                return { valid = false, reason = "point-charge singularity" }
            end
            local a = k * q.q / r
            potential = potential + a
            for j = 1, 3 do
                v[j] = v[j] + a * d[j] / (r * r)
            end
        end
        f = { V = potential }
    end
    if not f then
        return { valid = false, reason = "outside provider domain or on interface" }
    end
    local mag = norm(v)
    if not (mag < math.huge and math.abs(f.V) < math.huge) then
        return { valid = false, reason = "electric evaluation overflow" }
    end
    return {
        valid = true,
        V = f.V,
        vector = v,
        magnitude = mag,
        E = mag,
        Ex = v[1],
        Ey = v[2],
        Ez = v[3],
    }
end
function C.clearance(m, p)
    -- FE source-correction exclusions retain their native translational/revolved
    -- symmetry; the analytic coplanar point-charge provider uses genuine balls.
    if m.fe then
        if m.fe.axisymmetric then
            return E.clearance(m, sqrt(p[1] ^ 2 + p[2] ^ 2), p[3])
        end
        return E.clearance(m, p[1], p[2])
    end
    local d = math.huge
    for _, q in ipairs(m.sources) do
        d = math.min(d, norm({ p[1] - q.x, p[2] - q.y, p[3] }) - q.radius)
    end
    return d
end
local function map(paths, fn)
    local out = {}
    for _, path in ipairs(paths) do
        if #path.points >= 2 then
            local q = { points = {} }
            for _, p in ipairs(path.points) do
                q.points[#q.points + 1] = fn(p)
            end
            if path.closed then
                q.points[#q.points + 1] = q.points[1]
            end
            out[#out + 1] = q
        end
    end
    return out
end
local function circles(center, r, spatial, samples)
    local out = {}
    for a = 1, spatial and 3 or 1 do
        local points = {}
        for i = 0, samples do
            local u, v = r * math.cos(2 * pi * i / samples), r * math.sin(2 * pi * i / samples)
            local p = lift({ u, v }, ({ "xy", "xz", "yz" })[a], 0)
            for j = 1, 3 do
                p[j] = p[j] + center[j]
            end
            points[#points + 1] = p
        end
        out[#out + 1] = { points = points }
    end
    return out
end
local function geometry(m, spatial, plane, offset, view, o)
    local paths = {}
    local function append(a)
        for _, p in ipairs(a) do
            paths[#paths + 1] = p
        end
    end
    if m.fe or m.classic then
        local base = m.fe and require("luafemm-electric-fe").geometry(m)
            or require("luafemm-electric-classics").geometry(m, view.plate_span)
        if spatial and m.classic and not m.classic.planar then
            local c = m.classic
            for _, r in ipairs(c.interfaces) do
                if not c.cylindrical then
                    append(circles({ c.x, c.y, 0 }, r, true, o.geometry_samples))
                else
                    append(circles({ c.x, c.y, o.zmin }, r, false, o.geometry_samples))
                    append(circles({ c.x, c.y, o.zmax }, r, false, o.geometry_samples))
                    for i = 0, 3 do
                        local x, y = c.x + r * math.cos(i * pi / 2), c.y + r * math.sin(i * pi / 2)
                        append({ { points = { { x, y, o.zmin }, { x, y, o.zmax } } } })
                    end
                end
            end
        elseif spatial and m.fe and m.fe.axisymmetric then
            -- Four meridians are annotations, not opaque surfaces or a new volume mesh.
            for i = 0, 3 do
                local a = i * pi / 2
                append(map(base, function(p)
                    return { p[1] * math.cos(a), p[1] * math.sin(a), p[2] }
                end))
            end
        else
            append(map(base, function(p)
                return lift(p, spatial and "xy" or plane, spatial and 0 or offset)
            end))
        end
    else
        for _, q in ipairs(m.sources) do
            if spatial then
                append(circles({ q.x, q.y, 0 }, q.radius, true, o.geometry_samples))
            else
                append(
                    map(circles({ q.x, q.y, 0 }, q.radius, false, o.geometry_samples), function(p)
                        return lift(p, plane, offset)
                    end)
                )
            end
        end
    end
    return paths
end
local function seeds(m, text, view)
    if type(text) == "table" then
        return text
    end
    local out = {}
    if text and text ~= "" then
        for entry in text:gmatch("[^,]+") do
            local x, y, z = entry:match("^%s*([^/]+)/([^/]+)/([^/]+)%s*$")
            check(x, "spatial seeds must be comma-separated x/y/z triples")
            out[#out + 1] = { n(x, "seed x"), n(y, "seed y"), n(z, "seed z") }
        end
    else
        check(not m.fe, "spatial finite-element fields require explicit seeds")
        local opts = R.options(m, view)
        local count = opts.lines
        local centers = {}
        if m.classic then
            check(not m.classic.planar, "spatial plate fields require explicit seeds")
            local c = m.classic
            centers = { { x = c.x, y = c.y, radius = (c.radius or c.inner) } }
        else
            local sign = opts.seed_sign
            if sign == "auto" then
                sign = "negative"
                for _, q in ipairs(m.sources) do
                    if q.q > 0 then
                        sign = "positive"
                    end
                end
            end
            for _, q in ipairs(m.sources) do
                if
                    sign == "all"
                    or (sign == "positive" and q.q > 0)
                    or (sign == "negative" and q.q < 0)
                then
                    centers[#centers + 1] = q
                end
            end
        end
        check(#centers * count <= 500, "spatial seed budget exceeded (500)")
        for _, q in ipairs(centers) do
            for i = 0, count - 1 do
                local z = m.classic and m.classic.cylindrical and 0 or 1 - 2 * (i + 0.5) / count
                local r = sqrt(1 - z * z)
                local a = i * pi * (3 - sqrt(5))
                local d = q.radius * 1.05
                out[#out + 1] = { q.x + d * r * math.cos(a), q.y + d * r * math.sin(a), d * z }
            end
        end
    end
    return out
end
--- Capture transactionally into owned physical-coordinate primitives.
-- In spatial mode, contours are intersections with a world coordinate plane;
-- in slice mode, the plane only places the original two-dimensional section.
function C.capture(m, options, view)
    options = options or {}
    view = view or {}
    E.freeze(m)
    local mode = options.mode or "slice"
    local plane = options.plane or "xy"
    check(mode == "slice" or mode == "spatial", "electric scene mode must be slice or spatial")
    check(plane == "xy" or plane == "xz" or plane == "yz", "unknown electric scene plane")
    local offset = n(options.offset, "slice offset", 0)
    local o = R.options(m, view)
    o.geometry_samples = n(options.geometry_samples, "geometry samples", 96)
    check(
        o.geometry_samples % 1 == 0 and o.geometry_samples >= 8 and o.geometry_samples <= 1000,
        "geometry samples must be 8 to 1000"
    )
    o.zmin = n(options.zmin, "zmin", -3)
    o.zmax = n(options.zmax, "zmax", 3)
    check(o.zmin < o.zmax, "reversed spatial z bounds")
    local scene = S.new(m.descriptor.length_unit)
    local spatial = mode == "spatial"
    local function enabled(k, default)
        local v = options[k] or default
        check(v == "true" or v == "false", "electric scene " .. k .. " must be true or false")
        return v == "true"
    end
    local field, contours, vectors =
        enabled("field", "true"), enabled("equipotentials", "true"), enabled("vectors", "false")
    if enabled("geometry", "true") then
        S.add_paths(scene, "source", geometry(m, spatial, plane, offset, view, o))
    end
    local function embed(p)
        return lift(p, plane, offset)
    end
    if not spatial then
        for _, entry in ipairs({
            { "field", field, R.lines },
            { "equipotential", contours, R.contours },
            { "vectors", vectors, R.vectors },
        }) do
            if entry[2] then
                S.add_paths(scene, entry[1], map(entry[3](m, view), embed))
            end
        end
    else
        -- Only native providers with stated translation/rotation symmetry are lifted.
        local tracing = S.copy(view)
        for _, key in ipairs({
            "step",
            "tolerance",
            "max_length",
            "max_steps",
            "null_field",
            "direction",
        }) do
            tracing[key] = o[key]
        end
        tracing.xmin, tracing.ymin, tracing.xmax, tracing.ymax = table.unpack(o.bounds)
        tracing.zmin, tracing.zmax = o.zmin, o.zmax
        local points = (field or vectors) and seeds(m, options.seeds, view) or {}
        check(#points <= 500, "spatial seed budget exceeded (500)")
        for _, p in ipairs(points) do
            check(#p == 3, "seed needs three coordinates")
            for j, axis in ipairs({ "x", "y", "z" }) do
                p[j] = n(p[j], "seed " .. axis)
                check(
                    p[j] >= tracing[axis .. "min"] and p[j] <= tracing[axis .. "max"],
                    "seed outside the tracing box"
                )
            end
        end
        local function sample(p)
            if C.clearance(m, p) <= 0 then
                return { valid = false }
            end
            local f = C.sample(m, p)
            f.clearance = C.clearance(m, p)
            return f
        end
        if field then
            local lines = T.lines(sample, points, tracing)
            S.add_paths(
                scene,
                "field",
                map(lines, function(p)
                    return p
                end)
            )
        end
        if vectors then
            local paths = {}
            for _, p in ipairs(points) do
                local f = sample(p)
                if f.valid and f.E > o.null_field then
                    local l = o.vector_mode == "scaled" and f.E * o.vector_scale or o.vector_length
                    local q = {}
                    for j = 1, 3 do
                        q[j] = p[j] + l * f.vector[j] / f.E
                    end
                    if l < C.clearance(m, p) and sample(q).valid then
                        paths[#paths + 1] = { points = { p, q } }
                    end
                end
            end
            S.add_paths(scene, "vectors", paths)
        end
        if contours then
            -- Sample a world-coordinate slice. Discard a cell if its center is closer
            -- to an exclusion than its half-diagonal (conservative for distance masks).
            local sampler = {
                sample = function(x, y)
                    local f = C.sample(m, embed({ x, y }))
                    return f.valid and f or nil
                end,
                clearance = function(x, y)
                    return C.clearance(m, embed({ x, y }))
                end,
                mask = function(x0, y0, x1, y1)
                    return C.clearance(m, embed({ (x0 + x1) / 2, (y0 + y1) / 2 }))
                        <= sqrt((x1 - x0) ^ 2 + (y1 - y0) ^ 2) / 2
                end,
            }
            S.add_paths(scene, "equipotential", map(R.contours(m, view, sampler), embed))
        end
    end
    return scene
end
return C
