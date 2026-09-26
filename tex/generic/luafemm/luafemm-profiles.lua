-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- luafemm-profiles implementation.
-- @module luafemm-profiles
-- Arc-length sampling of connected TikZ paths and in-memory pgfplots output.
local P = {}
local profiles = {}
local abs, sqrt = math.abs, math.sqrt
local function finite(x)
    return type(x) == "number" and x == x and abs(x) < math.huge
end
local function check(v, s)
    assert(v, "luafemm profile: " .. s)
end
--- Register an arc-length profile, retaining its model independently of TeX groups.
-- Reusing a name replaces its geometry and invalidates its sample cache.
-- @tparam table m Source model.
-- @tparam string name Document-wide identifier.
-- @tparam table points Connected polyline in model units.
-- @tparam[opt] table o samples, offset and outside (error or nan).
-- @treturn table Registered profile.
function P.register(m, name, points, o)
    o = o or {}
    local n = o.samples or 201
    check(
        type(name) == "string" and name:match("^[%w_-]+$"),
        "name must contain only letters, digits, underscores or hyphens"
    )
    check(
        finite(n) and n >= 2 and n % 1 == 0 and n <= 100000,
        "samples must be an integer from 2 to 100000"
    )
    check(finite(o.offset or 0), "invalid normal offset")
    check(
        o.outside == nil or o.outside == "error" or o.outside == "nan",
        "outside must be error or nan"
    )
    local p, lengths, total = {}, {}, 0
    for _, v in ipairs(points) do
        check(finite(v[1]) and finite(v[2]), "invalid coordinate")
        local last = p[#p]
        local l = last and sqrt((v[1] - last[1]) ^ 2 + (v[2] - last[2]) ^ 2) or 0
        if not last or l > 0 then
            p[#p + 1] = { v[1], v[2] }
            total = total + l
            lengths[#p] = total
        end
    end
    check(total > 0, "path must have positive length")
    -- Named profiles retain their own model, even after another picture starts.
    profiles[name] = {
        model = m,
        points = p,
        lengths = lengths,
        length = total,
        samples = n,
        offset = o.offset or 0,
        outside = o.outside or "error",
    }
    return profiles[name]
end
--- Register a profile from an evaluated PGF soft path.
-- @tparam table m Model with a picture frame.
-- @tparam string name Profile identifier.
-- @tparam string softpath Evaluated PGF soft-path tokens.
-- @tparam table o Profile options, including curve tolerance.
function P.capture(m, name, softpath, o)
    check(m and m.frame, "set femm/problem on the source picture")
    local points = require("luafemm-path").read(softpath, m.frame, o.tolerance or 0.03, true)
    -- Same capture precision as physical regions, including domain endpoints.
    for _, p in ipairs(points) do
        for a = 1, 2 do
            p[a] = math.floor(p[a] * 10000 + 0.5) / 10000
            local lo = a == 1 and m.options.xmin or m.options.ymin
            local hi = a == 1 and m.options.xmax or m.options.ymax
            if abs(p[a] - lo) < 0.001 then
                p[a] = lo
            end
            if abs(p[a] - hi) < 0.001 then
                p[a] = hi
            end
        end
    end
    return P.register(m, name, points, o)
end
--- Look up a registered profile; unknown names raise an error.
-- @tparam string name Profile identifier.
-- @treturn table Internal profile record (treat as read-only).
function P.get(name)
    return assert(profiles[name], "luafemm profile: unknown profile " .. tostring(name))
end
--- Return all components at uniformly spaced arc-length positions.
-- Solves on demand and caches values. Normals point to the left of the path.
-- @tparam string name Registered profile name.
-- @treturn table Rows with s,t,x,y,B,Bx,By,Bt,Bn,H,Hx,Hy,Ht,Hn,Az.
function P.sample(name)
    local p = P.get(name)
    if p.rows then
        return p.rows
    end
    local f = require("luafemm")
    local m = p.model
    if not m.stats then
        f.solve(m)
    end
    local rows = {}
    local j = 2
    for i = 0, p.samples - 1 do
        local s = p.length * i / (p.samples - 1)
        while j < #p.points and s >= p.lengths[j] do
            j = j + 1
        end
        local a, b = p.points[j - 1], p.points[j]
        local l = p.lengths[j] - p.lengths[j - 1]
        local tx, ty = (b[1] - a[1]) / l, (b[2] - a[2]) / l
        local u = (s - p.lengths[j - 1]) / l
        local x, y =
            a[1] + u * (b[1] - a[1]) - p.offset * ty, a[2] + u * (b[2] - a[2]) + p.offset * tx
        local o = m.options
        local row = { s = s, t = i / (p.samples - 1), x = x, y = y }
        if x < o.xmin or x > o.xmax or y < o.ymin or y > o.ymax then
            check(
                p.outside == "nan",
                string.format("%s: sample %d at (%g,%g) outside domain", name, i + 1, x, y)
            )
            for _, k in ipairs({ "B", "Bx", "By", "Bt", "Bn", "H", "Hx", "Hy", "Ht", "Hn", "Az" }) do
                row[k] = 0 / 0
            end
        else
            local bx, by, az, hx, hy = f.sample(m, x, y)
            row.B, row.Bx, row.By = sqrt(bx * bx + by * by), bx, by
            row.H, row.Hx, row.Hy = sqrt(hx * hx + hy * hy), hx, hy
            row.Bt, row.Bn = bx * tx + by * ty, -bx * ty + by * tx
            row.Ht, row.Hn = hx * tx + hy * ty, -hx * ty + hy * tx
            row.Az = az
        end
        rows[#rows + 1] = row
    end
    p.rows = rows
    return rows
end
local components = {
    B = true,
    Bx = true,
    By = true,
    Bt = true,
    Bn = true,
    H = true,
    Hx = true,
    Hy = true,
    Ht = true,
    Hn = true,
    Az = true,
}
local abscissas = { s = true, t = true, x = true, y = true }
--- Format a chosen component for the pgfplots coordinates input handler.
-- @tparam string name Profile identifier.
-- @tparam[opt=B] string component Magnetic component or Az.
-- @tparam[opt=s] string abscissa s, t, x or y.
-- @treturn string Parenthesised numeric coordinate pairs; nan preserves gaps.
function P.coordinates(name, component, abscissa)
    component = component or "B"
    abscissa = abscissa or "s"
    check(components[component], "unknown component " .. component)
    check(abscissas[abscissa], "unknown abscissa " .. abscissa)
    local out = {}
    for _, r in ipairs(P.sample(name)) do
        out[#out + 1] = string.format(
            "(%.12g,%s)",
            r[abscissa],
            r[component] == r[component] and string.format("%.12g", r[component]) or "nan"
        )
    end
    return table.concat(out, " ")
end
return P
