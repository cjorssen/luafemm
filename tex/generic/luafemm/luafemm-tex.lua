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

--- Initialise the current picture from validated strings and the PGF transform.
function T.new_problem(options, transform, basis)
    for key, value in pairs(options) do
        if key ~= "mesher" then
            options[key] = T.number(value, key)
        end
    end
    options.log = function(message)
        texio.write_nl("luafemm: " .. message)
    end
    M.current = M.new(options)
    M.current.frame = require("luafemm-path").frame(transform, basis)
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
function T.capture_path(m, s, material, current, tolerance, angle, mesh_size)
    check(m and m.frame, "set femm/problem on the tikzpicture first")
    local points = require("luafemm-path").read(s, m.frame, tolerance)
    -- Round below the useful precision of PGF fixed-point coordinates.
    -- Snap domain endpoints to suppress TeX dimension-rounding artifacts.
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
    M.region(m, material, points, current, angle, mesh_size)
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
