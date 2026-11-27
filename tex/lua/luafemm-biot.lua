-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Prescribed DC currents in a homogeneous, isotropic 3D medium.
-- Geometry uses model units; all returned B/H components use SI. No FEM mesh,
-- material response, eddy currents or current redistribution is implied.
-- Thick curves are divergence-preserving radial convolutions, NOT an arbitrary
-- replacement of distance by sqrt(distance^2+a^2). Their straight-wire section
-- has exactly uniform current density. See kernel() and the manual derivation.
local B = {}
local G = require("luafemm-biot-geometry")
local copy = require("luafemm-physical").copy
local check, n = G.check, G.number
local add, sub, mul, dot, cross, norm = G.add, G.sub, G.mul, G.dot, G.cross, G.norm
local units = { m = 1, cm = 0.01, mm = 0.001, um = 1e-6 }
local pi = math.pi
local function positive(v, name, default)
    v = n(v, name, default)
    check(v > 0, name .. " must be positive")
    return v
end
local function integer(v, name, default, lo, hi)
    v = n(v, name, default)
    check(v % 1 == 0 and v >= lo and v <= hi, name .. " outside its integer budget")
    return v
end
--- New independent model. Relative/absolute tolerances control integration,
-- not geometry approximation or streamline accuracy.
function B.new(o)
    o = o or {}
    local u = o.length_unit or "cm"
    check(units[u], "unknown length unit")
    local options = {
        length_unit = u,
        mu_r = positive(o.mu_r, "relative permeability", 1),
        tolerance = positive(o.tolerance, "quadrature tolerance", 1e-5),
        absolute_tolerance = positive(o.absolute_tolerance, "absolute tolerance", 1e-12),
        max_depth = integer(o.max_depth, "maximum quadrature depth", 20, 1, 30),
        max_evaluations = integer(o.max_evaluations, "maximum evaluations", 200000, 24, 10000000),
        exclusion = positive(o.exclusion, "exclusion radius", 0.001),
    }
    check(options.tolerance < 1, "quadrature tolerance must be below one")
    return {
        biot = true,
        options = options,
        unit = units[u],
        sources = {},
        compiled = {},
        descriptor = {
            schema = 1,
            physics = "magnetism",
            physical_dimension = 3,
            evaluation_dimension = 3,
            geometry = "three dimensional",
            backend = "prescribed currents",
            potential = "none",
            length_unit = u,
            length_scale = units[u],
            reference = "infinity",
            units = { length = "m", B = "T", H = "A/m" },
        },
    }
end
function B.add(m, o)
    check(not m.frozen, "declare sources before capturing or rendering the model")
    check(#m.sources < 2000, "maximum source count exceeded")
    local r = G.record(o)
    local curves = G.compile(r)
    local panels = 0
    for _, group in ipairs(m.compiled) do
        for _, c in ipairs(group) do
            panels = panels + c.panels
        end
    end
    for _, c in ipairs(curves) do
        check(c.speed > 0 and c.speed < math.huge, "degenerate or overflowing source geometry")
        panels = panels + c.panels
    end
    check(panels <= 100000, "source geometry panel budget exceeded")
    m.sources[#m.sources + 1] = r
    m.compiled[#m.compiled + 1] = curves
    return r
end
-- Enclosed mass of K_a(r)=1/(pi^2 a^2 sqrt(a^2-r^2)), r<a.
-- Its transverse projection is 1/(pi a^2) on a disk. Convolving the
-- filament current with K_a commutes with divergence and curl. Newton's
-- radial-kernel identity gives the enclosed fraction multiplying R/|R|^3.
local function kernel(distance, a)
    if a == 0 or distance >= a then
        return 1 / distance ^ 3
    end
    local q = distance / a
    if q < 1e-3 then
        return (4 / (3 * pi) + 2 * q * q / (5 * pi) + 3 * q ^ 4 / (14 * pi)) / a ^ 3
    end
    return 2 / pi * (math.asin(q) - q * math.sqrt(1 - q * q)) / distance ^ 3
end
B.kernel = kernel
local gx4 = { 0.3399810435848563, 0.8611363115940526 }
local gw4 = { 0.6521451548625461, 0.3478548451374538 }
local gx8 = { 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 }
local gw8 = { 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 }
local function quadrature(fn, a, b, x, w, state)
    local c, h = (a + b) / 2, (b - a) / 2
    local sum = { 0, 0, 0 }
    for i, t in ipairs(x) do
        state.evaluations = state.evaluations + 2
        check(state.evaluations <= state.limit, "quadrature evaluation budget exhausted")
        local p, q = fn(c - h * t), fn(c + h * t)
        sum = add(sum, mul(add(p, q), h * w[i]))
    end
    return sum
end
-- Refine the largest error globally within one curve. Assigning vanishing
-- absolute tolerances to individual children is unreliable at a finite-core
-- surface, whose kernel has a square-root endpoint derivative.
local function integrate(fn, knots, abs_tol, o, state)
    local cells = {}
    local sum = { 0, 0, 0 }
    local error = 0
    local function cell(a, b, depth)
        local q4 = quadrature(fn, a, b, gx4, gw4, state)
        local q8 = quadrature(fn, a, b, gx8, gw8, state)
        return { a = a, b = b, depth = depth, value = q8, error = norm(sub(q8, q4)) }
    end
    for i = 2, #knots do
        if knots[i] > knots[i - 1] then
            local c = cell(knots[i - 1], knots[i], 0)
            cells[#cells + 1] = c
            sum = add(sum, c.value)
            error = error + c.error
        end
    end
    while error > abs_tol + o.tolerance * norm(sum) do
        local index = 1
        for i = 2, #cells do
            if cells[i].error > cells[index].error then
                index = i
            end
        end
        local c = cells[index]
        check(c.depth < o.max_depth, "quadrature depth exhausted; refine tolerances or geometry")
        local mid = (c.a + c.b) / 2
        local l, r = cell(c.a, mid, c.depth + 1), cell(mid, c.b, c.depth + 1)
        sum = add(sub(sum, c.value), add(l.value, r.value))
        error = math.max(0, error - c.error + l.error + r.error)
        cells[index] = l
        cells[#cells + 1] = r
    end
    return sum, error
end
local function distance_line(c, p)
    local d = sub(c.b, c.a)
    local t = dot(sub(p, c.a), d) / dot(d, d)
    if not c.infinite then
        t = math.max(0, math.min(1, t))
    end
    return norm(sub(p, add(c.a, mul(d, t)))), t
end
-- A speed bound encloses each analytic curve panel. Only nearby panels need
-- minimisation. Small base panels prevent a many-turn winding aliasing one
-- quadrature panel or concealing a source intersection between quadrature nodes.
local function excluded(c, p, r)
    if c.line then
        return distance_line(c, p) <= r
    end
    for i = 1, c.panels do
        local a, b = (i - 1) / c.panels, i / c.panels
        local q = c.evaluate((a + b) / 2)
        if norm(sub(p, q)) <= c.speed * (b - a) / 2 + r then
            local ratio = (math.sqrt(5) - 1) / 2
            local x, y = b - ratio * (b - a), a + ratio * (b - a)
            local function d(t)
                return norm(sub(p, c.evaluate(t)))
            end
            local dx, dy = d(x), d(y)
            for _ = 1, 48 do
                if dx < dy then
                    b, y, dy = y, x, dx
                    x = b - ratio * (b - a)
                    dx = d(x)
                else
                    a, x, dx = x, y, dy
                    y = a + ratio * (b - a)
                    dy = d(y)
                end
            end
            if math.min(dx, dy, d(a), d(b)) <= r then
                return true
            end
        end
    end
    return false
end
local function straight(c, p, a)
    local d = sub(c.b, c.a)
    local length = norm(d)
    local u = mul(d, 1 / length)
    local r = sub(p, c.a)
    local z = dot(r, u)
    local perpendicular = sub(r, mul(u, z))
    local rho2 = dot(perpendicular, perpendicular)
    if c.infinite then
        if rho2 == 0 then
            return { 0, 0, 0 }
        end
        return mul(cross(u, perpendicular), 2 / math.max(rho2, a * a))
    end
    if rho2 == 0 then
        return { 0, 0, 0 }
    end
    local r1, r2 = math.sqrt(rho2 + z * z), math.sqrt(rho2 + (z - length) ^ 2)
    if z < 0 or z > length then
        -- Stable endpoint expression outside the segment: do not subtract two
        -- nearly identical direction cosines in a distant axial observation.
        local h = length * (r1 + r2) / (r1 * r2 * (r1 * r2 + rho2 + z * (z - length)))
        return mul(cross(u, perpendicular), h)
    end
    local h = z / r1 - (z - length) / r2
    return mul(cross(u, perpendicular), h / rho2)
end
--- Evaluate one point. Excluded filament positions return a validity flag;
-- numerical budget failures are errors and are never turned into plotted zeros.
function B.sample(m, p)
    p = G.vector(p, "sample position")
    local total = { 0, 0, 0 }
    local error = 0
    local state = { evaluations = 0, limit = m.options.max_evaluations }
    local count = 0
    for _, curves in ipairs(m.compiled) do
        count = count + #curves
    end
    for j, r in ipairs(m.sources) do
        local factor = 1e-7 * m.options.mu_r * r.current / m.unit
        if factor ~= 0 then
            for _, c in ipairs(m.compiled[j]) do
                if r.wire_radius == 0 and excluded(c, p, m.options.exclusion) then
                    return { valid = false, reason = "filament exclusion", source = j }
                end
                local value, err
                if
                    c.infinite
                    or (c.line and (r.wire_radius == 0 or distance_line(c, p) >= r.wire_radius))
                then
                    value, err = straight(c, p, r.wire_radius), 0
                else
                    local function fn(t)
                        local q, d = c.evaluate(t)
                        local delta = sub(p, q)
                        local distance = norm(delta)
                        if distance == 0 then
                            return { 0, 0, 0 }
                        end
                        return mul(cross(d, delta), kernel(distance, r.wire_radius))
                    end
                    local knots = { 0, 1 }
                    if c.line then
                        local _, t = distance_line(c, p)
                        -- Resolve a narrow core before testing nested rules.
                        for _, dt in ipairs({ -2, -1, 0, 1, 2 }) do
                            local x = t + dt * r.wire_radius / c.speed
                            if x > 0 and x < 1 then
                                knots[#knots + 1] = x
                            end
                        end
                    else
                        for i = 1, c.panels - 1 do
                            knots[#knots + 1] = i / c.panels
                        end
                    end
                    table.sort(knots)
                    value, err = integrate(
                        fn,
                        knots,
                        m.options.absolute_tolerance / math.abs(factor) / math.max(1, count),
                        m.options,
                        state
                    )
                end
                total = add(total, mul(value, factor))
                error = error + math.abs(factor) * err
            end
        end
    end
    local mu = 4 * pi * 1e-7 * m.options.mu_r
    local result = {
        valid = true,
        Bx = total[1],
        By = total[2],
        Bz = total[3],
        B = norm(total),
        Hx = total[1] / mu,
        Hy = total[2] / mu,
        Hz = total[3] / mu,
        H = norm(total) / mu,
        vector = total,
        estimated_error = error,
        evaluations = state.evaluations,
    }
    -- Overflow is a failed evaluation, not a drawable field or an exclusion gap.
    for _, key in ipairs({ "Bx", "By", "Bz", "B", "Hx", "Hy", "Hz", "H", "estimated_error" }) do
        local value = result[key]
        check(
            value == value and math.abs(value) < math.huge,
            "non-finite field; rescale geometry or current"
        )
    end
    return result
end
function B.snapshot(m)
    local p = { options = copy(m.options), sources = copy(m.sources) }
    return B.restore(p)
end
--- Revalidate passive source records and reconstruct executable evaluators.
function B.restore(p)
    check(
        type(p) == "table" and type(p.options) == "table" and type(p.sources) == "table",
        "invalid prescribed-current records"
    )
    local m = B.new(p.options)
    for key in pairs(p) do
        check(key == "options" or key == "sources", "unknown prescribed-current problem field")
    end
    for key in pairs(p.options) do
        check(m.options[key] ~= nil, "unknown prescribed-current option: " .. tostring(key))
    end
    local count = 0
    for k in pairs(p.sources) do
        count = count + 1
        check(
            type(k) == "number" and k % 1 == 0 and k >= 1 and k <= #p.sources,
            "invalid source index"
        )
    end
    check(count == #p.sources, "sparse source records")
    for _, r in ipairs(p.sources) do
        check(type(r) == "table" and type(r.rotation) == "table", "invalid source record")
        local o = copy(r)
        o.rotate_x, o.rotate_y, o.rotate_z = table.unpack(G.vector(r.rotation, "rotation"))
        o.phase = r.phase and n(r.phase, "phase") * 180 / pi or 0
        local canonical = G.record(o)
        for key in pairs(r) do
            check(
                canonical[key] ~= nil,
                "unknown prescribed-current source field: " .. tostring(key)
            )
        end
        B.add(m, o)
    end
    m.frozen = true
    return m
end
return B
