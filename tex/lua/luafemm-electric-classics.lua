-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Classical electrostatic symmetry models, evaluated in a physical section.
-- Cylinders and plates are infinite. Spheres are genuinely three dimensional.
-- A model contains one distribution or one ideal capacitor: adding prescribed
-- point charges would generally violate the electrode conditions and is rejected.
-- Coordinates/radii use model units; densities, voltages and results use SI.
-- @module luafemm-electric-classics
local S = {}
local abs, sqrt, min, max = math.abs, math.sqrt, math.min, math.max
local function check(v, s)
    assert(v, "luafemm electric symmetry: " .. s)
end
local function number(v, name, default)
    if v == nil or v == "" then
        v = default
    end
    v = tonumber(v)
    check(v and v == v and abs(v) < math.huge, "finite value required for " .. name)
    return v
end
local function positive(v, name, default)
    v = number(v, name, default)
    check(v > 0, name .. " must be positive")
    return v
end
local kinds = {
    ["charged sphere"] = true,
    ["charged cylinder"] = true,
    ["parallel plate capacitor"] = true,
    ["cylindrical capacitor"] = true,
    ["spherical capacitor"] = true,
}
--- Validate the complete physical declaration before drawing anything.
function S.configure(m, o)
    local kind = o.model
    check(kinds[kind], "unknown analytical model: " .. tostring(kind))
    local c = {
        kind = kind,
        x = number(o.center_x, "center x", 0),
        y = number(o.center_y, "center y", 0),
        reference = number(o.reference_potential, "reference potential", 0),
        surface_value = o.interface_value or "upper",
    }
    check(c.surface_value == "lower" or c.surface_value == "upper", "unknown interface value")
    c.capacitor = kind:find("capacitor", 1, true) ~= nil
    c.planar = kind == "parallel plate capacitor"
    c.cylindrical = kind == "charged cylinder" or kind == "cylindrical capacitor"
    if c.capacitor then
        check(not o.volume_density or o.volume_density == "", "density is not a capacitor input")
        check(not o.surface_density or o.surface_density == "", "density is not a capacitor input")
        c.voltage = number(o.voltage, "voltage")
        if c.planar then
            c.separation = positive(o.separation, "separation", 1)
            c.angle = number(o.normal_angle, "normal angle", 0)
            c.nx, c.ny = math.cos(math.rad(c.angle)), math.sin(math.rad(c.angle))
            c.interfaces = { -c.separation / 2, c.separation / 2 }
        else
            c.inner = positive(o.inner_radius, "inner radius", 0.7)
            c.outer = positive(o.outer_radius, "outer radius", 1.4)
            check(c.inner < c.outer, "inner radius must be smaller than outer radius")
            c.interfaces = { c.inner, c.outer }
        end
    else
        check(not o.voltage or o.voltage == "", "voltage requires a capacitor")
        c.radius = positive(o.radius, "radius", 1)
        c.distribution = o.distribution or "volume"
        check(
            c.distribution == "volume" or c.distribution == "surface",
            "unknown charge distribution"
        )
        local selected = c.distribution == "volume" and "volume_density" or "surface_density"
        local unused = c.distribution == "volume" and "surface_density" or "volume_density"
        check(not o[unused] or o[unused] == "", "supply only the selected charge density")
        c.density = number(o[selected], (selected:gsub("_", " ")))
        if c.cylindrical then
            c.reference_radius = positive(o.reference_radius, "reference radius")
        else
            check(
                not o.reference_radius or o.reference_radius == "",
                "reference radius requires a charged cylinder"
            )
        end
        c.interfaces = { c.radius }
    end
    if c.capacitor then
        check(
            not o.reference_radius or o.reference_radius == "",
            "reference radius requires a charged cylinder"
        )
    end
    m.classic = c
    local d = m.descriptor
    d.geometry = kind
    d.symmetry = c.planar and "translation in two directions"
        or (c.cylindrical and "infinite cylinder" or "spherical")
    d.reference = c.capacitor and (c.planar and "positive-normal electrode" or "outer electrode")
        or (c.cylindrical and "finite reference radius" or "infinity")
    d.reference_potential = c.reference
    d.reference_radius = c.reference_radius
    d.validity = "all finite section points; explicit one-sided field values at interfaces"
    d.units.volume_charge_density = "C/m^3"
    d.units.surface_charge_density = "C/m^2"
end
local function coordinate(c, x, y)
    x, y = x - c.x, y - c.y
    if c.planar then
        return x * c.nx + y * c.ny, x, y
    end
    return sqrt(x * x + y * y), x, y
end
-- Exact equality is interpreted with a relative geometric tolerance, rather
-- than allowing roundoff in an intersected circle to select an arbitrary side.
local function below(q, r, side, gap)
    local tolerance = 1e-11 * max(abs(r), 1e-100)
    if gap then
        tolerance = min(tolerance, gap / 8)
    end
    if abs(q - r) <= tolerance then
        return side == "lower"
    end
    return q < r
end
-- Stable logarithmic ratio for thin coaxial gaps and nearby radii.
local function logratio(a, b)
    local u = (a - b) / b
    if abs(u) < 1e-4 then
        return u * (1 - u * (0.5 - u * (1 / 3 - u * (0.25 - u / 5))))
    end
    return math.log(a / b)
end
local function radial_body(m, c, r, side)
    local R = c.radius * m.unit
    local inside = below(r, R, side)
    if c.cylindrical then
        local a = c.distribution == "volume" and c.density * R * R / (2 * m.epsilon)
            or c.density * R / m.epsilon
        local function potential(t)
            if t < R then
                return c.distribution == "volume" and a * (1 - t * t / (R * R)) / 2 or 0
            end
            return -a * logratio(t, R)
        end
        local er = inside and (c.distribution == "volume" and a * r / (R * R) or 0) or a / r
        return c.reference + potential(r) - potential(c.reference_radius * m.unit), er
    end
    local a = c.distribution == "volume" and c.density * R ^ 3 / (3 * m.epsilon)
        or c.density * R * R / m.epsilon
    if inside then
        if c.distribution == "volume" then
            return c.reference + a * (3 - r * r / (R * R)) / (2 * R), a * r / R ^ 3
        end
        return c.reference + a / R, 0
    end
    return c.reference + a / r, a / (r * r)
end
--- One-sided field convention: lower/upper means smaller/larger radial or
-- signed normal coordinate. V is continuous and independent of this choice.
function S.sample(m, x, y, side)
    local c = m.classic
    side = side or c.surface_value
    check(side == "lower" or side == "upper", "unknown interface value")
    local q, dx, dy = coordinate(c, x, y)
    local v, er
    if not c.capacitor then
        v, er = radial_body(m, c, q * m.unit, side)
    elseif c.planar then
        local d = c.separation
        if below(q, -d / 2, side) then
            v, er = c.reference + c.voltage, 0
        elseif below(q, d / 2, side) then
            v, er = c.reference + c.voltage * (0.5 - q / d), c.voltage / (d * m.unit)
        else
            v, er = c.reference, 0
        end
    else
        local a, b = c.inner, c.outer
        if below(q, a, side, b - a) then
            v, er = c.reference + c.voltage, 0
        elseif below(q, b, side, b - a) then
            if c.cylindrical then
                local factor = c.voltage / logratio(b, a)
                v, er = c.reference + factor * logratio(b, q), factor / (q * m.unit)
            else
                local factor = c.voltage * a * b / (b - a)
                v, er = c.reference + factor * (b - q) / (q * b), factor / (q * q * m.unit)
            end
        else
            v, er = c.reference, 0
        end
    end
    local nx, ny = 0, 0
    if c.planar then
        nx, ny = c.nx, c.ny
    elseif q > 0 then
        nx, ny = dx / q, dy / q
    end
    local ex, ey = er * nx, er * ny
    if not (abs(v) < math.huge and abs(er) < math.huge and abs(m.epsilon * er) < math.huge) then
        return nil, "analytical evaluation overflow"
    end
    return {
        V = v,
        Ex = ex,
        Ey = ey,
        E = abs(er),
        Dx = m.epsilon * ex,
        Dy = m.epsilon * ey,
        D = m.epsilon * abs(er),
    }
end
--- SI integral quantities, with their normalization encoded in the key.
function S.properties(m)
    local c = m.classic
    local eps, u = m.epsilon, m.unit
    if c.capacitor then
        local cap, suffix
        if c.planar then
            cap, suffix = eps / (c.separation * u), " per area"
        elseif c.cylindrical then
            cap, suffix = 2 * math.pi * eps / logratio(c.outer, c.inner), " per length"
        else
            cap, suffix = 4 * math.pi * eps * u * c.inner * c.outer / (c.outer - c.inner), ""
        end
        return {
            ["capacitance" .. suffix] = cap,
            ["charge" .. suffix] = cap * c.voltage,
            ["energy" .. suffix] = cap * c.voltage ^ 2 / 2,
            voltage = c.voltage,
        }
    end
    local r = c.radius * u
    if c.cylindrical then
        return {
            ["charge per length"] = c.distribution == "volume" and math.pi * r * r * c.density
                or 2 * math.pi * r * c.density,
        }
    end
    return {
        charge = c.distribution == "volume" and 4 * math.pi * r ^ 3 * c.density / 3
            or 4 * math.pi * r * r * c.density,
    }
end
--- Intersections of a segment with physical interfaces; zero-length segments
-- and lines lying in a plane interface need no artificial crossing event.
function S.intersections(m, a, b)
    local c = m.classic
    local dx, dy = b[1] - a[1], b[2] - a[2]
    local out = {}
    local aa = dx * dx + dy * dy
    if aa == 0 then
        return out
    end
    for _, r in ipairs(c.interfaces) do
        local roots = {}
        if c.planar then
            local slope = dx * c.nx + dy * c.ny
            if abs(slope) > 1e-15 * sqrt(aa) then
                roots[1] = (r - coordinate(c, a[1], a[2])) / slope
            end
        else
            local x, y = a[1] - c.x, a[2] - c.y
            local bb = x * dx + y * dy
            local cc = x * x + y * y - r * r
            local discriminant = bb * bb - aa * cc
            if discriminant >= 0 then
                roots = { (-bb - sqrt(discriminant)) / aa, (-bb + sqrt(discriminant)) / aa }
            end
        end
        for _, t in ipairs(roots) do
            if t >= -1e-12 and t <= 1 + 1e-12 then
                out[#out + 1] = { t = max(0, min(1, t)), surface = r }
            end
        end
    end
    return out
end
--- Insert exact path-interface positions and one-sided values in profiles.
-- This prevents a finite sample interval from smoothing a field jump, even
-- when the complete dielectric gap fits between two requested sample points.
function S.profile_events(m, points, lengths)
    local total = lengths[#lengths]
    local hits = {}
    for i = 2, #points do
        for _, event in ipairs(S.intersections(m, points[i - 1], points[i])) do
            hits[#hits + 1] = {
                s = lengths[i - 1] + event.t * (lengths[i] - lengths[i - 1]),
                surface = event.surface,
            }
        end
    end
    table.sort(hits, function(a, b)
        return a.s < b.s
    end)
    local unique = {}
    for _, hit in ipairs(hits) do
        if
            #unique == 0
            or abs(hit.s - unique[#unique].s) > 1e-12 * total
            or hit.surface ~= unique[#unique].surface
        then
            unique[#unique + 1] = hit
        end
    end
    local function point(s)
        local i = 2
        while i < #points and lengths[i] < s do
            i = i + 1
        end
        local len = lengths[i] - lengths[i - 1]
        local t = len > 0 and (s - lengths[i - 1]) / len or 0
        return points[i - 1][1] + t * (points[i][1] - points[i - 1][1]),
            points[i - 1][2] + t * (points[i][2] - points[i - 1][2])
    end
    for i, hit in ipairs(unique) do
        local previous = i > 1 and unique[i - 1].s or 0
        local following = i < #unique and unique[i + 1].s or total
        -- At an endpoint only the limit along the existing contour is used.
        local dl = min(total * 1e-7, (hit.s - previous) / 4)
        local dr = min(total * 1e-7, (following - hit.s) / 4)
        local function side(s)
            local x, y = point(s)
            return coordinate(m.classic, x, y) < hit.surface and "lower" or "upper"
        end
        hit.left = dl > 0 and side(hit.s - dl) or nil
        hit.right = dr > 0 and side(hit.s + dr) or nil
        hit.left = hit.left or hit.right or m.classic.surface_value
        hit.right = hit.right or hit.left
        hit.jump = (
            m.classic.capacitor and m.classic.voltage ~= 0
            or m.classic.distribution == "surface" and m.classic.density ~= 0
        ) and hit.left ~= hit.right
    end
    return unique
end
--- True discontinuities only: omit arrows crossing them rather than shortening.
function S.crosses_jump(m, a, b)
    if not (m.classic.capacitor or m.classic.distribution == "surface") then
        return false
    end
    for _, e in ipairs(S.intersections(m, a, b)) do
        if e.t > 1e-12 and e.t < 1 - 1e-12 then
            return true
        end
    end
    return false
end
local function line_box(b, p, n)
    local low, high = -math.huge, math.huge
    for i = 1, 2 do
        if abs(n[i]) < 1e-15 then
            if p[i] < b[i] or p[i] > b[i + 2] then
                return nil
            end
        else
            local a, z = (b[i] - p[i]) / n[i], (b[i + 2] - p[i]) / n[i]
            if a > z then
                a, z = z, a
            end
            low, high = max(low, a), min(high, z)
        end
    end
    if low > high then
        return nil
    end
    return low, high
end
--- Exact straight radial/normal field lines. Zero-field interiors are not
-- traced. RK controls remain available for point-charge superpositions.
function S.lines(m, o)
    local c = m.classic
    local seeds = o.seeds or {}
    local result = {}
    local strength = c.capacitor and c.voltage or c.density
    if strength == 0 then
        return result
    end
    local active_low, active_high
    if c.planar then
        if abs(c.voltage) / (c.separation * m.unit) <= o.null_field then
            return result
        end
        active_low, active_high = -c.separation / 2, c.separation / 2
    else
        local radius = c.capacitor and c.inner or c.radius
        local f =
            assert(S.sample(m, c.x + radius, c.y, "upper"), "luafemm: analytical field overflow")
        if f.E <= o.null_field then
            return result
        end
        active_low = c.capacitor and c.inner
            or (c.distribution == "surface" and c.radius or radius * o.null_field / f.E)
        active_high = radius * (f.E / o.null_field) ^ (c.cylindrical and 1 or 0.5)
        if c.capacitor then
            active_high = min(active_high, c.outer)
        end
    end
    local sign = strength > 0 and 1 or -1
    if not o.seeds then
        if (o.seed_sign == "positive" and sign < 0) or (o.seed_sign == "negative" and sign > 0) then
            return result
        end
        if c.planar then
            local b = o.bounds
            local lo, hi = math.huge, -math.huge
            for _, p in ipairs({ { b[1], b[2] }, { b[1], b[4] }, { b[3], b[2] }, { b[3], b[4] } }) do
                local t = -(p[1] - c.x) * c.ny + (p[2] - c.y) * c.nx
                lo, hi = min(lo, t), max(hi, t)
            end
            for i = 1, o.lines do
                local t = lo + (hi - lo) * (i - 0.5) / o.lines
                seeds[#seeds + 1] = { c.x - t * c.ny, c.y + t * c.nx, automatic = true }
            end
        else
            local r = c.capacitor and (c.inner + c.outer) / 2 or c.radius * 1.05
            for i = 1, o.lines do
                local a = 2 * math.pi * (i - 0.5) / o.lines
                seeds[#seeds + 1] =
                    { c.x + r * math.cos(a), c.y + r * math.sin(a), automatic = true }
            end
        end
    end
    for _, seed in ipairs(seeds) do
        local q, dx, dy = coordinate(c, seed[1], seed[2])
        local p, n, low, high
        if c.planar then
            p = { seed[1] - q * c.nx, seed[2] - q * c.ny }
            n = { c.nx, c.ny }
            low, high = active_low, active_high
        elseif q > 0 then
            p = { c.x, c.y }
            n = { dx / q, dy / q }
            low, high = active_low, active_high
        end
        if p then
            local a, b = line_box(o.bounds, p, n)
            local explicit_ok = seed.automatic
                or (
                    seed[1] >= o.bounds[1]
                    and seed[1] <= o.bounds[3]
                    and seed[2] >= o.bounds[2]
                    and seed[2] <= o.bounds[4]
                    and q > low
                    and q < high
                )
            if a and explicit_ok then
                low, high = max(low, a), min(high, b)
                local start = max(low, min(high, q))
                if o.direction == "both" then
                    low, high = max(low, start - o.max_length), min(high, start + o.max_length)
                elseif (o.direction == "forward") == (sign > 0) then
                    low, high = start, min(high, start + o.max_length)
                else
                    low, high = max(low, start - o.max_length), start
                end
                if low < high then
                    local points = {
                        { p[1] + low * n[1], p[2] + low * n[2] },
                        {
                            p[1] + high * n[1],
                            p[2] + high * n[2],
                        },
                    }
                    if sign < 0 then
                        points = { points[2], points[1] }
                    end
                    result[#result + 1] = {
                        points = points,
                        start = "analytic endpoint",
                        finish = "analytic endpoint",
                    }
                end
            end
        end
    end
    return result
end
--- Model-derived geometry; plate span is a display length, never a finite plate.
function S.geometry(m, span)
    local c = m.classic
    local out = {}
    local function circle(r, style)
        local points = {}
        for i = 0, 180 do
            local a = 2 * math.pi * i / 180
            points[#points + 1] = { c.x + r * math.cos(a), c.y + r * math.sin(a) }
        end
        out[#out + 1] = { points = points, style = style, closed = true }
    end
    if c.planar then
        span = positive(span, "plate span", 4)
        for _, x in ipairs(c.interfaces) do
            out[#out + 1] = {
                points = {
                    { c.x + x * c.nx + span * c.ny / 2, c.y + x * c.ny - span * c.nx / 2 },
                    { c.x + x * c.nx - span * c.ny / 2, c.y + x * c.ny + span * c.nx / 2 },
                },
                style = "femm electric electrode",
            }
        end
    elseif c.capacitor then
        circle(c.outer, "femm electric electrode")
        circle(c.inner, "femm electric conductor")
    else
        circle(
            c.radius,
            c.distribution == "volume" and "femm charged body" or "femm charged surface"
        )
    end
    return out
end
return S
