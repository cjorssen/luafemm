-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Scoped generic TeX adapter for prescribed-current models and owned profiles.
-- Coordinates and options are passive numeric strings, never evaluated Lua.
local T = {}
local B = require("luafemm-biot")
local G = require("luafemm-biot-geometry")
local V = require("luafemm-biot-view")
local P = require("luafemm-path")
local M = require("luafemm")
local models, names, previous, profiles = {}, {}, {}, {}
local check, n = G.check, G.number
local function model(id)
    return assert(models[tonumber(id)], "luafemm: no Biot-Savart picture")
end
local function identifier(s)
    check(type(s) == "string" and s:match("^[%w][%w:_%-]*$"), "invalid result name")
    return s
end
function T.new(o, transform, basis)
    local m
    if o.native_file ~= "" then
        check(o.explicit ~= "1", "native file owns the physical and numerical settings")
        m = require("luafemm-tex").native_read(o.native_file)
        check(m.biot, "native file is not a prescribed-current model")
    else
        m = B.new(o)
    end
    m.frame = P.frame(transform, basis)
    local id = #models + 1
    models[id] = m
    previous[id] = M.current
    M.current = nil
    tex.sprint("\\expandafter\\def\\csname luafemm@biotid\\endcsname{" .. id .. "}")
end
function T.finish(id, name)
    if name ~= "" then
        names[identifier(name)] = B.snapshot(model(id))
    end
    M.current = previous[tonumber(id)]
end
function T.native_model(id, name)
    if name and name ~= "" then
        return assert(names[identifier(name)], "unknown Biot-Savart model")
    end
    return model(id)
end
function T.has_model(name)
    return names[name] ~= nil
end
function T.source(id, o)
    B.add(model(id), o)
end
local function lift(p, plane, offset)
    if plane == "xy" then
        return { p[1], p[2], offset }
    elseif plane == "xz" then
        return { p[1], offset, p[2] }
    elseif plane == "yz" then
        return { offset, p[1], p[2] }
    end
    error("luafemm: unknown source plane")
end
function T.path(id, path, o)
    local m = model(id)
    local points = P.read(path, m.frame, n(o.curve_tolerance, "curve tolerance", 0.002), true)
    o.points = {}
    o.kind = "wire"
    for _, p in ipairs(points) do
        o.points[#o.points + 1] = lift(p, o.plane, n(o.offset, "offset", 0))
    end
    B.add(m, o)
end
local function vectors(s)
    local result = {}
    for point in tostring(s):gmatch("[^;]+") do
        result[#result + 1] = G.vector(point, "point list")
    end
    return result
end
local function draw(m, paths, o, style)
    local f = m.frame
    for _, path in ipairs(paths) do
        if #path.points >= 2 then
            local out = { "\\draw[" .. style .. "] " }
            for i, p in ipairs(path.points) do
                local q = V.project(p, o)
                out[#out + 1] = string.format(
                    "%s(%.7fpt,%.7fpt)",
                    i == 1 and "" or "--",
                    f[1] * q[1] + f[3] * q[2] + f[5],
                    f[2] * q[1] + f[4] * q[2] + f[6]
                )
            end
            out[#out + 1] = ";"
            tex.sprint(table.concat(out))
        end
    end
end
-- Aesthetic seeds describe the first source only; they carry no flux weight.
local function automatic_seeds(m)
    local r = m.sources[1]
    check(r, "declare a source before drawing field lines")
    local seeds = {}
    if r.kind == "loop" or r.kind == "solenoid" then
        for _, x in ipairs({ -0.65, -0.3, 0, 0.3, 0.65 }) do
            seeds[#seeds + 1] = G.place(r, { x * r.radius, 0, 0 })
        end
    elseif r.kind == "toroid" then
        for _, z in ipairs({ -0.45, 0, 0.45 }) do
            seeds[#seeds + 1] = G.place(r, { r.major_radius, 0, z * r.radius })
        end
    elseif #r.points == 2 then
        local u = G.sub(r.points[2], r.points[1])
        u = G.mul(u, 1 / G.norm(u))
        local normal = G.cross(u, math.abs(u[3]) < 0.9 and { 0, 0, 1 } or { 1, 0, 0 })
        normal = G.mul(normal, 1 / G.norm(normal))
        local mid = G.mul(G.add(r.points[1], r.points[2]), 0.5)
        for _, a in ipairs({ 0.6, 1, 1.4 }) do
            seeds[#seeds + 1] =
                G.place(r, G.add(mid, G.mul(normal, a * math.max(1, 3 * r.wire_radius))))
        end
    else
        local center = { 0, 0, 0 }
        local count = #r.points - 1
        for i = 1, count do
            center = G.add(center, r.points[i])
        end
        seeds[1] = G.place(r, G.mul(center, 1 / count))
    end
    return seeds
end
local function view_paths(m, kind, o)
    local paths, style
    if kind == "geometry" then
        paths, style = V.geometry(m, o), "femm biot conductor"
    elseif kind == "vectors" then
        paths, style = V.vectors(m, o), "femm biot vector"
    else
        paths, style =
            V.lines(m, o.seeds == "auto" and automatic_seeds(m) or vectors(o.seeds), o),
            "femm biot field"
        for _, p in ipairs(paths) do
            if p.start == "steps" or p.finish == "steps" then
                texio.write_nl(
                    "term and log",
                    "luafemm Biot-Savart warning: streamline reached maximum steps"
                )
                break
            end
        end
    end
    return paths, style
end
function T.draw(id, kind, o)
    local m = model(id)
    m.frozen = true
    local paths, style = view_paths(m, kind, o)
    draw(m, paths, V.options(o), style)
end
--- Capture owned spatial primitives once. Subsequent cameras do not evaluate B.
function T.scene(id, options, view)
    local m = T.native_model(id, options.model)
    local S = require("luafemm-scene3d")
    local scene = S.new(m.options.length_unit)
    V.options(view)
    for _, entry in ipairs({
        { "geometry", "geometry" },
        { "field", "lines" },
        { "vectors", "vectors" },
    }) do
        local enabled = options[entry[1]]
        check(enabled == "true" or enabled == "false", "scene selection must be true or false")
        if enabled == "true" then
            S.add_paths(scene, entry[1], view_paths(m, entry[2], view))
        end
    end
    require("luafemm-scene3d-tex").store(identifier(options.name), scene)
    m.frozen = true
end
--- A profile owns evaluated SI samples, not a reference to a later builder.
function T.profile(id, o)
    local m = T.native_model(id, o.model)
    m.frozen = true
    local points = o.points ~= "" and vectors(o.points)
        or { G.vector(o.from, "from"), G.vector(o.to, "to") }
    check(#points >= 2, "profile needs at least two points")
    local samples = n(o.samples, "samples", 101)
    check(samples % 1 == 0 and samples >= 2 and samples <= 10001, "samples must be 2 to 10001")
    local cumulative = { 0 }
    for i = 2, #points do
        cumulative[i] = cumulative[i - 1] + G.norm(G.sub(points[i], points[i - 1]))
    end
    local length = cumulative[#points]
    check(length > 0, "zero-length profile")
    local rows = {}
    local segment = 2
    for i = 0, samples - 1 do
        local s = i * length / (samples - 1)
        while segment < #points and cumulative[segment] < s do
            segment = segment + 1
        end
        local a, b = points[segment - 1], points[segment]
        local size = cumulative[segment] - cumulative[segment - 1]
        check(size > 0, "profile has repeated points")
        local p = G.add(a, G.mul(G.sub(b, a), (s - cumulative[segment - 1]) / size))
        local f = B.sample(m, p)
        f.x, f.y, f.z, f.s = p[1], p[2], p[3], s
        rows[#rows + 1] = f
    end
    profiles[identifier(o.name)] = { rows = rows, unit = m.unit }
end
local components =
    { B = true, Bx = true, By = true, Bz = true, H = true, Hx = true, Hy = true, Hz = true }
function T.plot(o)
    local p = assert(profiles[identifier(o.profile)], "unknown Biot-Savart profile")
    check(components[o.component], "unknown field component")
    check(
        o.abscissa == "s" or o.abscissa == "x" or o.abscissa == "y" or o.abscissa == "z",
        "unknown abscissa"
    )
    local out = {}
    local scale = n(o.scale, "abscissa scale", 1)
    check(scale > 0, "abscissa scale must be positive")
    for _, f in ipairs(p.rows) do
        out[#out + 1] = string.format(
            "(%.12g,%s)",
            f[o.abscissa] * scale,
            f.valid and string.format("%.12g", f[o.component]) or "nan"
        )
    end
    tex.sprint(table.concat(out, " "))
end
function T.value(id, o)
    check(components[o.component], "unknown field component")
    local m = T.native_model(id, o.model)
    m.frozen = true
    local f = B.sample(m, G.vector(o.at, "at"))
    check(f.valid, "field value lies in a filament exclusion")
    tex.sprint(string.format("%.12g", f[o.component]))
end
return T
