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
-- @tparam[opt] table o samples, offset, outside and polar centre/origin/pole count.
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
        if not last or (l > 0 and total + l > total) then
            p[#p + 1] = { v[1], v[2] }
            total = total + l
            lengths[#p] = total
        end
    end
    check(total > 0, "path must have positive length")
    -- Named profiles retain their own model, even after another picture starts.
    for _, k in ipairs({ "center_x", "center_y", "angle_origin", "pole_pairs" }) do
        check(o[k] == nil or finite(o[k]), "invalid polar option " .. k)
    end
    check(
        (o.pole_pairs or 1) >= 1 and (o.pole_pairs or 1) % 1 == 0,
        "pole pairs must be a positive integer"
    )
    profiles[name] = {
        model = m,
        points = p,
        lengths = lengths,
        length = total,
        samples = n,
        offset = o.offset or 0,
        center_x = o.center_x or 0,
        center_y = o.center_y or 0,
        angle_origin = o.angle_origin or 0,
        pole_pairs = o.pole_pairs or 1,
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
        if p.circle then
            local a = math.rad(p.angle_origin + 360 * i / (p.samples - 1))
            tx, ty = -math.sin(a), math.cos(a)
            x, y = p.center_x + p.circle * math.cos(a), p.center_y + p.circle * math.sin(a)
            s = 2 * math.pi * p.circle * i / (p.samples - 1)
        end
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
            local ok, bx, by, az, hx, hy = pcall(f.sample, m, x, y)
            if not ok then
                check(
                    p.outside == "nan" and tostring(bx):find("sample outside domain", 1, true),
                    tostring(bx)
                )
                bx, by, az, hx, hy = 0 / 0, 0 / 0, 0 / 0, 0 / 0, 0 / 0
            end
            row.B, row.Bx, row.By = sqrt(bx * bx + by * by), bx, by
            row.H, row.Hx, row.Hy = sqrt(hx * hx + hy * hy), hx, hy
            row.Bt, row.Bn = bx * tx + by * ty, -bx * ty + by * tx
            row.Ht, row.Hn = hx * tx + hy * ty, -hx * ty + hy * tx
            row.Az = az
        end
        local dx, dy = x - p.center_x, y - p.center_y
        local radius = math.sqrt(dx * dx + dy * dy)
        if radius == 0 then
            radius = 0 / 0
        end
        row.Br, row.Btheta =
            (row.Bx * dx + row.By * dy) / radius, (-row.Bx * dy + row.By * dx) / radius
        row.Hr, row.Htheta =
            (row.Hx * dx + row.Hy * dy) / radius, (-row.Hx * dy + row.Hy * dx) / radius
        local angle = math.deg(math.atan(dy, dx)) - p.angle_origin
        if #rows > 0 then
            local previous = rows[#rows].angle
            angle = previous + (angle - previous + 180) % 360 - 180
        end
        row.angle = p.circle and 360 * i / (p.samples - 1) or angle
        row.electrical_angle = p.pole_pairs * row.angle
        rows[#rows + 1] = row
    end
    p.rows = rows
    return rows
end
local components = {
    Br = true,
    Btheta = true,
    Hr = true,
    Htheta = true,
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
    circulation = true,
}
local abscissas = { s = true, t = true, x = true, y = true, angle = true, electrical_angle = true }
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
    local rows = P.sample(name)
    if component == "circulation" then
        local p = P.get(name)
        check(not p.circle, "use an ordinary captured path for circulation")
        check(p.offset == 0, "circulation requires offset=0; declare the displaced path explicitly")
        if not p.integrals then
            local integrate = require("luafemm-integrals").segment
            local base, j = 0, 2
            for _, row in ipairs(rows) do
                while j < #p.points and row.s >= p.lengths[j] do
                    base = base + integrate(p.model, p.points[j - 1], p.points[j])
                    j = j + 1
                end
                row.circulation = base + integrate(p.model, p.points[j - 1], { row.x, row.y })
            end
            p.integrals = true
        end
    end
    local out = {}
    for _, r in ipairs(rows) do
        out[#out + 1] = string.format(
            "(%.12g,%s)",
            r[abscissa],
            r[component] == r[component] and string.format("%.12g", r[component]) or "nan"
        )
    end
    return table.concat(out, " ")
end

--- Register an exact circular scan, sampled uniformly in mechanical angle.
-- The closing endpoint is retained for plotting, excluded from Fourier sums.
-- @tparam table m Source model.
-- @tparam string name Document-wide profile identifier.
-- @tparam number x Centre x in model units.
-- @tparam number y Centre y in model units.
-- @tparam number radius Positive radius in model units.
-- @tparam[opt] table o Sampling and polar options; offset must be zero.
-- @treturn table Registered profile retaining its source model.
function P.circle(m, name, x, y, radius, o)
    o = o or {}
    check(not o.offset or o.offset == 0, "circular profiles require zero offset")
    check(finite(radius) and radius > 0, "positive circle radius required")
    local settings = {}
    for k, v in pairs(o) do
        settings[k] = v
    end
    settings.center_x, settings.center_y = x, y
    local points = {}
    for i = 0, 72 do
        local a = math.rad((o.angle_origin or 0) + 5 * i)
        points[#points + 1] = { x + radius * math.cos(a), y + radius * math.sin(a) }
    end
    local p = P.register(m, name, points, settings)
    p.circle = radius
    p.length = 2 * math.pi * radius
    return p
end

--- Fourier coefficients of a full circular profile in the mechanical angle.
-- Phase is for A*cos(n*theta+phase), degrees relative to angle_origin.
-- THD includes orders 1..maximum except the requested fundamental; DC is separate.
-- @tparam string name Full-circle profile name.
-- @tparam[opt=Br] string component Field component to analyse.
-- @tparam[opt=15] integer maximum Largest mechanical harmonic below Nyquist.
-- @tparam[opt=1] integer fundamental Reference order for distortion.
-- @treturn table Indexed coefficients plus mean and thd (NaN for zero fundamental).
function P.harmonics(name, component, maximum, fundamental)
    component = component or "Br"
    maximum = maximum or 15
    fundamental = fundamental or 1
    check(components[component] and component ~= "circulation", "invalid harmonic component")
    local p = P.get(name)
    check(p.circle, "harmonics require a circular gap profile")
    local rows = P.sample(name)
    local n = #rows - 1
    check(
        maximum >= 1 and maximum % 1 == 0 and 2 * maximum < n,
        "harmonic order must be below Nyquist"
    )
    check(
        fundamental >= 1 and fundamental % 1 == 0 and fundamental <= maximum,
        "invalid fundamental order"
    )
    local mean = 0
    for i = 1, n do
        check(finite(rows[i][component]), "nonfinite harmonic sample")
        mean = mean + rows[i][component] / n
    end
    local out = { mean = mean }
    for k = 1, maximum do
        local a, b = 0, 0
        for i = 1, n do
            local theta = 2 * math.pi * (i - 1) / n
            a = a + 2 * rows[i][component] * math.cos(k * theta) / n
            b = b + 2 * rows[i][component] * math.sin(k * theta) / n
        end
        out[k] = {
            cosine = a,
            sine = b,
            amplitude = math.sqrt(a * a + b * b),
            phase = math.deg(math.atan(-b, a)),
        }
    end
    local energy = 0
    for k = 1, maximum do
        if k ~= fundamental then
            energy = energy + out[k].amplitude ^ 2
        end
    end
    out.thd = out[fundamental].amplitude > 1e-30 and math.sqrt(energy) / out[fundamental].amplitude
        or 0 / 0
    return out
end
return P
