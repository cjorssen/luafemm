-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- LuaTeX boundary: validated key values, path capture and generated TeX.
-- Numerical modules do not access TeX globals. This module deliberately does.
-- @module luafemm-tex
local T = {}
local M = require("luafemm")
local abs, max = math.abs, math.max
local function check(value, message)
    assert(value, "luafemm: " .. message)
end

--- Parse a finite numeric key without evaluating user text as Lua code.
-- @tparam string value Numeric literal supplied through luaescapestring.
-- @tparam string key Name included in diagnostics.
-- @param default Value returned for an empty string, if supplied.
function T.number(value, key, default)
    if value == "" and default ~= nil then
        return default
    end
    local number = tonumber(value)
    check(
        number and number == number and abs(number) < math.huge,
        "expected a finite number for " .. key .. ", got " .. tostring(value)
    )
    return number
end

--- Apply one key-based exterior condition to the current model.
function T.boundary(side, options)
    for key, value in pairs(options) do
        if key ~= "type" then
            options[key] = T.number(value, "boundary " .. key)
        end
    end
    M.boundary(M.current, side, options)
end

--- Initialise the current picture from validated strings and the PGF transform.
function T.new_problem(options, transform, basis, cache)
    for key, value in pairs(options) do
        if key ~= "mesher" then
            options[key] = T.number(value, key)
        end
    end
    options.log = function(message)
        texio.write_nl("luafemm: " .. message)
    end
    -- Follow TeX's output directory for relative cache names. The standalone
    -- Lua API deliberately retains ordinary working-directory semantics.
    if
        cache
        and cache.file
        and status
        and status.output_directory
        and status.output_directory ~= ""
        and not cache.file:match("^[/\\]")
        and not cache.file:match("^%a:")
    then
        cache.file = status.output_directory .. "/" .. cache.file
    end
    options.cache = cache
    M.current = M.new(options)
    M.current.frame = require("luafemm-path").frame(transform, basis)
    -- Shape dimensions use the original coordinate basis, independently of
    -- text metrics and later local x/y vector changes.
    M.current.component_basis = require("luafemm-path").frame("{1}{0}{0}{1}{0pt}{0pt}", basis)
end

--- Print the current cache outcome; this command never starts a calculation.
function T.cache_status()
    check(M.current, "no current model")
    tex.sprint(M.current.cache_info.status)
end

--- Resolve one PGF material style; omitted overrides inherit library properties.
function T.material(name, library, mur, bh, remanence, coercivity)
    local function optional(value, key)
        if value ~= "" then
            return T.number(value, key)
        end
    end
    local points = M.parse_pairs(bh)
    M.material(M.current, name, {
        library = library ~= "" and library or nil,
        mur = optional(mur, "mu r"),
        bh = next(points) and points or nil,
        remanence = optional(remanence, "remanence"),
        coercivity = optional(coercivity, "coercivity"),
    })
end

--- Print a solved-model statistic as plain numeric TeX input.
function T.statistic(key, mesh)
    local m = M.current
    check(m, "no current model")
    local source
    if mesh then
        source = m.mesh_stats
    else
        source = m.stats
    end
    local value = source and source[key]
    check(type(value) == "number", "statistic unavailable: " .. key)
    tex.sprint(string.format(mesh and "%.6g" or "%.8g", value))
end

--- Print the field magnitude in tesla at a point in model units.
function T.field_value(x, y)
    local bx, by = M.sample(M.current, T.number(x, "x"), T.number(y, "y"))
    tex.sprint(string.format("%.3f", math.sqrt(bx * bx + by * by)))
end

--- Emit field contours for the legacy coordinate-based TeX API.
-- @tparam table m Solved model.
-- @tparam number count Number of potential levels.
function T.tikz_lines(m, count)
    local paths = M.contours(m, count)
    local out = {}
    for _, path in ipairs(paths) do
        local parts = { "\\draw[femm field] " }
        for i, p in ipairs(path) do
            parts[#parts + 1] = string.format("%s(%.7f,%.7f)", i > 1 and "--" or "", p.x, p.y)
        end
        parts[#parts + 1] = ";"
        out[#out + 1] = table.concat(parts)
    end
    tex.sprint(table.concat(out, " "))
    return paths
end

--- Capture a closed material path and apply the PGF coordinate-precision convention.
-- @tparam table m Unmeshed model with a frame.
-- @tparam string s Evaluated soft-path tokens.
-- @tparam string material Model material name.
-- @tparam number current Current density in A/m^2.
-- @tparam number tolerance Polygonalisation tolerance in model units.
-- @tparam number angle Magnetisation direction in degrees.
-- @tparam number mesh_size Local spacing, zero for inheritance.
-- @tparam[opt] string turns Component winding ampere-turns (single contour only).
-- @tparam[opt="nonzero"] string fill_rule Evaluated native PGF fill rule.
function T.capture_path(m, s, material, current, tolerance, angle, mesh_size, turns, fill_rule)
    check(m and m.frame, "set femm/problem on the tikzpicture first")
    local contours = require("luafemm-path").read_contours(s, m.frame, tolerance)
    -- Round below the useful precision of PGF fixed-point coordinates.
    -- Snap domain endpoints to suppress TeX dimension-rounding artifacts.
    for _, points in ipairs(contours) do
        for _, p in ipairs(points) do
            for axis = 1, 2 do
                p[axis] = math.floor(p[axis] * 10000 + 0.5) / 10000
                local lo = axis == 1 and m.options.xmin or m.options.ymin
                local hi = axis == 1 and m.options.xmax or m.options.ymax
                if abs(p[axis] - lo) < 0.001 then
                    p[axis] = lo
                end
                if abs(p[axis] - hi) < 0.001 then
                    p[axis] = hi
                end
            end
        end
    end
    if turns and turns ~= "" then
        check(#contours == 1, "component winding needs one contour")
        current = require("luafemm-components").current(
            contours[1],
            m.options.unit,
            T.number(turns, "ampere turns")
        )
    end
    M.region_contours(m, material, contours, current, angle, mesh_size, fill_rule)
end

-- Saved shape IDs keep anchor geometry independent of the current model.
-- These records hold no solver/model reference and never enter disk caches.
local components = {}

--- Create immutable node geometry; PGF stores the returned numeric ID.
-- @tparam string kind Component type.
-- @tparam table options Escaped numeric strings from component keys.
function T.component_new(kind, options)
    check(M.current and M.current.frame, "set femm/problem before declaring components")
    check(not M.current.nodes, "declare components before meshing")
    for key, value in pairs(options) do
        options[key] = T.number(value, key)
    end
    local c = require("luafemm-components").make(kind, options)
    c.basis = M.current.component_basis
    components[#components + 1] = c
    tex.sprint(tostring(#components))
end

local function component_point(c, p)
    local f = c.basis
    return f[1] * p[1] + f[3] * p[2], f[2] * p[1] + f[4] * p[2]
end

--- Emit one saved anchor in the shape's local PGF coordinate system.
function T.component_anchor(id, name)
    local c = assert(components[id])
    local p = c.anchors[name]
    check(p, "anchor '" .. name .. "' is unavailable on " .. c.kind)
    local x, y = component_point(c, p)
    tex.sprint(string.format("\\pgfpoint{%.9fpt}{%.9fpt}", x, y))
end

--- Intersect a connection ray with the component envelope, not its iron.
function T.component_border(id, x, y)
    local c = assert(components[id])
    local f = c.basis
    local function dimension(s)
        return assert(tonumber(s:match("([%d%.%+%-eE]+)pt")))
    end
    x, y = dimension(x), dimension(y)
    local det = f[1] * f[4] - f[2] * f[3]
    local u, v = (f[4] * x - f[3] * y) / det, (-f[2] * x + f[1] * y) / det
    local scale = max(abs(u) / c.half_width, abs(v) / c.half_height)
    if scale == 0 then
        tex.sprint("\\pgfpointorigin")
    else
        local px, py = component_point(c, { u / scale, v / scale })
        tex.sprint(string.format("\\pgfpoint{%.9fpt}{%.9fpt}", px, py))
    end
end

--- Draw/capture each part once, at the final node transform, before its text.
-- A component-local magnetization direction is mapped into the model frame.
-- Coil ampere-turns are converted after capture, using the rounded polygon.
function T.component_paths(id, transform)
    local c = assert(components[id])
    local n = require("luafemm-path").frame(transform, { "1pt", "0pt", "0pt", "1pt" })
    local f = M.current.frame
    local theta = math.rad(c.options.magnetization_angle)
    local vx, vy = component_point(c, { math.cos(theta), math.sin(theta) })
    vx, vy = n[1] * vx + n[3] * vy, n[2] * vx + n[4] * vy
    local det = f[1] * f[4] - f[2] * f[3]
    local mx, my = (f[4] * vx - f[3] * vy) / det, (-f[2] * vx + f[1] * vy) / det
    check(mx * mx + my * my > 0, "singular component transform")
    local angle = math.deg(math.atan(my, mx))
    for _, r in ipairs(c.regions) do
        local points = {}
        for i, p in ipairs(r.points) do
            local x, y = component_point(c, p)
            points[i] = string.format("(%.9fpt,%.9fpt)", x, y)
        end
        tex.sprint(
            string.format(
                "\\csname luafemm@componentpart\\endcsname{%s}{%s}{%.17g}{%.17g}{%s--cycle}",
                r.role,
                r.turns and string.format("%.17g", r.turns) or "",
                r.mesh_size,
                angle,
                table.concat(points, "--")
            )
        )
    end
end

local function canvas(m, x, y)
    local f = m.frame
    return string.format("(%.7fpt,%.7fpt)", f[1] * x + f[3] * y + f[5], f[2] * x + f[4] * y + f[6])
end

--- Solve if needed and emit contours in the original picture frame.
-- @tparam table m Model with frame.
-- @tparam number count Number of contour levels.
function T.tikz_field_picture(m, count)
    if not m.stats then
        M.solve(m)
    end
    for _, path in ipairs(M.contours(m, count)) do
        local out = { "\\draw[femm field,femm current field] " }
        for i, p in ipairs(path) do
            out[#out + 1] = (i > 1 and "--" or "") .. canvas(m, p.x, p.y)
        end
        out[#out + 1] = ";"
        tex.sprint(table.concat(out))
    end
end

--- Mesh if needed and emit each edge exactly once in the picture frame.
-- @tparam table m Model with frame.
function T.tikz_mesh_picture(m)
    if not m.nodes then
        M.mesh(m)
    end
    local seen = {}
    local out = { "\\draw[femm mesh,femm current mesh] " }
    for _, e in ipairs(m.triangles) do
        for i = 1, 3 do
            local a, b = e.ids[i], e.ids[i % 3 + 1]
            local k = math.min(a, b) .. ":" .. max(a, b)
            if not seen[k] then
                seen[k] = true
                local p, q = m.nodes[a], m.nodes[b]
                local u = m.options.unit
                out[#out + 1] = canvas(m, p[1] / u, p[2] / u)
                    .. "--"
                    .. canvas(m, q[1] / u, q[2] / u)
                    .. " "
            end
        end
    end
    out[#out + 1] = ";"
    tex.sprint(table.concat(out))
end

-- Identifiers contain only ASCII letters, digits and hyphens: safe PGF keys.
--- Install safe ASCII identifiers as PGF choice handlers at library load time.
-- The generic loader calls this with @ a letter; no document data are emitted.
function T.register_material_choices()
    for _, r in ipairs(require("luafemm-materials").list()) do
        tex.sprint(
            "\\pgfkeys{/femm/region/material/"
                .. r.id
                .. "/.code={\\luafemm@selectlibrary{"
                .. r.id
                .. "}}}"
        )
        tex.sprint(
            "\\pgfkeys{/tikz/femm/materials/library/"
                .. r.id
                .. "/.code={\\def\\luafemm@library{"
                .. r.id
                .. "}}}"
        )
    end
end

return T
