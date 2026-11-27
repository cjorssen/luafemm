-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Provider-neutral adaptive spatial streamline integration.
-- The callback returns valid, vector and magnitude in consistent SI field units.
-- Coordinates and all geometric tolerances use the model's length unit.
local T = {}
local S = require("luafemm-scene3d")
local check, n = S.check, S.number
local G = {}
function G.vector(p)
    check(type(p) == "table" and #p == 3, "seed needs three coordinates")
    return { n(p[1], "seed"), n(p[2], "seed"), n(p[3], "seed") }
end
local function add(a, b)
    return { a[1] + b[1], a[2] + b[2], a[3] + b[3] }
end
local function sub(a, b)
    return { a[1] - b[1], a[2] - b[2], a[3] - b[3] }
end
local function mul(a, s)
    return { a[1] * s, a[2] * s, a[3] * s }
end
local function norm(a)
    return math.sqrt(a[1] ^ 2 + a[2] ^ 2 + a[3] ^ 2)
end
function T.options(o)
    o = o or {}
    local r = {}
    for k, v in pairs({
        step = 0.08,
        tolerance = 0.0002,
        max_length = 30,
        max_steps = 3000,
        null_field = 1e-14,
        xmin = -3,
        xmax = 3,
        ymin = -3,
        ymax = 3,
        zmin = -3,
        zmax = 3,
    }) do
        r[k] = n(o[k], k, v)
    end
    check(
        r.step > 0 and r.tolerance > 0 and r.max_length > 0 and r.null_field >= 0,
        "invalid tracing controls"
    )
    check(
        r.max_steps % 1 == 0 and r.max_steps >= 1 and r.max_steps <= 100000,
        "invalid maximum steps"
    )
    for _, a in ipairs({ "x", "y", "z" }) do
        check(r[a .. "max"] > r[a .. "min"], "reversed trace bounds")
    end
    r.direction = o.direction or "both"
    check(
        r.direction == "both" or r.direction == "forward" or r.direction == "backward",
        "unknown tracing direction"
    )
    return r
end
local function inside(p, o)
    return p[1] >= o.xmin
        and p[1] <= o.xmax
        and p[2] >= o.ymin
        and p[2] <= o.ymax
        and p[3] >= o.zmin
        and p[3] <= o.zmax
end
local function trace(sample, seed, sign, o)
    local points = { seed }
    local length, h = 0, o.step
    local reason
    local function direction(p)
        local f = sample(p)
        if not f.valid then
            reason = "source"
            return
        end
        if f.magnitude <= o.null_field then
            reason = "null"
            return
        end
        return mul(f.vector, sign / f.magnitude)
    end
    local function rk(p, s)
        local a = direction(p)
        if not a then
            return
        end
        local b = direction(add(p, mul(a, s / 2)))
        if not b then
            return
        end
        local c = direction(add(p, mul(b, s / 2)))
        if not c then
            return
        end
        local d = direction(add(p, mul(c, s)))
        if not d then
            return
        end
        return add(p, mul(add(add(a, mul(b, 2)), add(mul(c, 2), d)), s / 6))
    end
    for _ = 1, o.max_steps do
        local p = points[#points]
        if not inside(p, o) then
            return points, "window"
        end
        local first = sample(p)
        if first.clearance then
            if first.clearance <= o.tolerance then
                return points, "source"
            end
            h = math.min(h, first.clearance * 0.45)
        end
        local q
        local accepted = false
        for _ = 1, 30 do
            local full = rk(p, h)
            if not full then
                return points, reason
            end
            local half = rk(p, h / 2)
            if not half then
                return points, reason
            end
            q = rk(half, h / 2)
            if not q then
                return points, reason
            end
            local error = norm(sub(q, full)) / 15
            if error <= o.tolerance then
                accepted = true
                if error < o.tolerance / 32 then
                    h = math.min(o.step, h * 2)
                end
                break
            end
            h = h / 2
        end
        check(accepted, "streamline tolerance cannot be reached")
        local dl = norm(sub(q, p))
        if dl < 1e-12 then
            return points, "stationary"
        end
        if length + dl > o.max_length then
            return points, "length"
        end
        points[#points + 1] = q
        length = length + dl
        if length > 8 * o.step and norm(sub(q, seed)) < 0.5 * o.step then
            points[#points + 1] = seed
            return points, "closed"
        end
    end
    return points, "steps"
end
--- Return oriented physical 3D polylines. Backward branches are reversed so
-- that every emitted arrow follows the supplied field, including direction=backward requests.
function T.lines(sample, seeds, o)
    local result = {}
    o = T.options(o)
    check(type(seeds) == "table" and #seeds > 0 and #seeds <= 500, "provide 1 to 500 seeds")
    for _, seed in ipairs(seeds) do
        seed = G.vector(seed, "seed")
        check(inside(seed, o), "seed outside the tracing box")
        local back, forward, rs, rf = { seed }, { seed }, "seed", "seed"
        if o.direction ~= "forward" then
            back, rs = trace(sample, seed, -1, o)
        end
        if o.direction ~= "backward" and rs ~= "closed" then
            forward, rf = trace(sample, seed, 1, o)
        end
        local points = {}
        for i = #back, 1, -1 do
            points[#points + 1] = back[i]
        end
        for i = 2, #forward do
            points[#points + 1] = forward[i]
        end
        result[#result + 1] = { points = points, start = rs, finish = rf }
    end
    return result
end
return T
