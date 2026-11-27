-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Owned dielectric catalogue records and safe generic PGF choice installation.
-- Typical FEMM 5 properties, not certified design data. The data remain MIT;
-- this original adapter is LPPL. No files are parsed at document compilation.
-- @module luafemm-dielectrics
local D = {}
local copy = require("luafemm-physical").copy
local records = require("luafemm-dielectrics-data")
local by_id = {}
for _, r in ipairs(records) do
    by_id[r.id] = r
end
function D.list()
    return copy(records)
end
function D.get(id)
    return copy(assert(by_id[id], "luafemm dielectric: unknown material: " .. tostring(id)))
end
--- Return the scalar-compatible record, or nil for the user's numeric inputs.
function D.scalar(id)
    if id == nil or id == "" or id == "custom" then
        return nil
    end
    local r = D.get(id)
    assert(r.ex == r.ey, "luafemm dielectric: anisotropic permittivity is not supported: " .. id)
    return r
end
--- Symmetric positive relative tensor in the physical x/y (or r/z) frame.
-- Principal values are dimensionless; angle is counterclockwise in degrees.
-- A scalar override supplies both principal values; explicit axes override it.
function D.tensor(o, fallback)
    local E = require("luafemm-electrostatic")
    local function present(x)
        return x ~= nil and x ~= ""
    end
    local r = present(o.material) and o.material ~= "custom" and D.get(o.material) or nil
    local scalar = present(o.epsilon_r) and E.positive(o.epsilon_r, "relative permittivity") or nil
    fallback = fallback or { xx = 1, xy = 0, yy = 1 }
    if
        not r
        and not scalar
        and not present(o.epsilon_x)
        and not present(o.epsilon_y)
        and not present(o.material_angle)
    then
        return copy(fallback)
    end
    local x = E.positive(
        present(o.epsilon_x) and o.epsilon_x
            or scalar
            or (r and r.ex)
            or fallback.principal_x
            or fallback.xx,
        "relative permittivity x"
    )
    local y = E.positive(
        present(o.epsilon_y) and o.epsilon_y
            or scalar
            or (r and r.ey)
            or fallback.principal_y
            or fallback.yy,
        "relative permittivity y"
    )
    local angle = E.number(present(o.material_angle) and o.material_angle or 0, "material angle")
        * math.pi
        / 180
    local c, s = math.cos(angle), math.sin(angle)
    return {
        xx = x * c * c + y * s * s,
        xy = (x - y) * c * s,
        yy = x * s * s + y * c * c,
        principal_x = x,
        principal_y = y,
    }
end
--- Generate only vetted ASCII identifiers and finite numeric literals as TeX.
-- Material selection resets numeric defaults; subsequent explicit keys override
-- them according to ordinary PGF order and grouping. Custom only clears the id.
function D.register()
    for _, family in ipairs({ "electrostatic", "electric region" }) do
        for _, r in ipairs(records) do
            assert(r.id:match("^[a-z0-9%-]+$"), "unsafe dielectric identifier")
            local path = "/femm/" .. family .. "/"
            local body = "\\pgfkeyssetvalue{"
                .. path
                .. "material}{"
                .. r.id
                .. "}"
                .. "\\pgfkeyssetvalue{"
                .. path
                .. "relative permittivity}{"
                .. r.ex
                .. "}"
            body = body
                .. "\\pgfkeyssetvalue{"
                .. path
                .. "relative permittivity x}{"
                .. r.ex
                .. "}"
                .. "\\pgfkeyssetvalue{"
                .. path
                .. "relative permittivity y}{"
                .. r.ey
                .. "}"
                .. "\\pgfkeyssetvalue{"
                .. path
                .. "material angle}{0}"
            if family == "electric region" then
                body = body
                    .. "\\pgfkeyssetvalue{"
                    .. path
                    .. "volume charge density}{"
                    .. r.qv
                    .. "}"
            else
                body = body .. "\\def\\luafemm@classicexplicit{1}"
            end
            tex.sprint("\\pgfkeys{" .. path .. "material/" .. r.id .. "/.code={" .. body .. "}}")
        end
    end
end
return D
