-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- luafemm-predicates implementation.
-- @module luafemm-predicates
-- Filtered geometric predicates with exact floating-point expansions.
-- Original Lua implementation of error-free sum/product arithmetic.
-- References and numerical range assumptions: docs/meshing.md.
local P = { stats = { orientation = 0, incircle = 0 } }
local abs = math.abs
local function two_sum(a, b)
    local s = a + b
    local v = s - a
    return s, (a - (s - v)) + (b - v)
end
local function two_product(a, b)
    local p = a * b
    local c = 134217729.0 * a
    local ah = c - (c - a)
    local al = a - ah
    c = 134217729.0 * b
    local bh = c - (c - b)
    local bl = b - bh
    return p, ((ah * bh - p) + ah * bl + al * bh) + al * bl
end
local function grow(e, b)
    local out = {}
    local q = b
    for _, v in ipairs(e) do
        local s, r = two_sum(q, v)
        if r ~= 0 then
            out[#out + 1] = r
        end
        q = s
    end
    if q ~= 0 or #out == 0 then
        out[#out + 1] = q
    end
    return out
end
local function sum(e, f, sign)
    local out = e
    sign = sign or 1
    for _, v in ipairs(f) do
        out = grow(out, sign * v)
    end
    return out
end
local function product(e, f)
    local out = { 0 }
    for _, a in ipairs(e) do
        for _, b in ipairs(f) do
            local p, r = two_product(a, b)
            out = grow(out, r)
            out = grow(out, p)
        end
    end
    return out
end
local function diff(a, b)
    local h, l = two_sum(a + 0.0, -b)
    return l ~= 0 and { l, h } or { h }
end
local function determinant(ax, ay, bx, by)
    return sum(product(ax, by), product(ay, bx), -1)
end
local function leading(e)
    return e[#e]
end
--- Return a determinant estimate with the exact orientation sign of binary inputs.
-- @tparam table a First {x,y} point.
-- @tparam table b Second point.
-- @tparam table c Third point.
-- @treturn number Positive for counterclockwise orientation, zero if collinear.
function P.orient(a, b, c)
    local ax, ay = a[1] + 0.0 - c[1], a[2] + 0.0 - c[2]
    local bx, by = b[1] + 0.0 - c[1], b[2] + 0.0 - c[2]
    local left, right = ax * by, ay * bx
    local det = left - right
    if abs(det) > 8.881784197001252e-16 * (abs(left) + abs(right)) then
        return det
    end
    P.stats.orientation = P.stats.orientation + 1
    return leading(
        determinant(diff(a[1], c[1]), diff(a[2], c[2]), diff(b[1], c[1]), diff(b[2], c[2]))
    )
end
--- Return the exact incircle sign, subject to the documented exponent assumptions.
-- @tparam table a First point of a counterclockwise triangle.
-- @tparam table b Second triangle point.
-- @tparam table c Third triangle point.
-- @tparam table d Query point.
-- @treturn number Positive inside the circumcircle; negative outside; zero on it.
function P.incircle(a, b, c, d)
    local ax, ay = a[1] + 0.0 - d[1], a[2] + 0.0 - d[2]
    local bx, by = b[1] + 0.0 - d[1], b[2] + 0.0 - d[2]
    local cx, cy = c[1] + 0.0 - d[1], c[2] + 0.0 - d[2]
    local al, bl, cl = ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy
    local det = al * (bx * cy - by * cx) + bl * (cx * ay - cy * ax) + cl * (ax * by - ay * bx)
    local bound = al * (abs(bx * cy) + abs(by * cx))
        + bl * (abs(cx * ay) + abs(cy * ax))
        + cl * (abs(ax * by) + abs(ay * bx))
    if abs(det) > 7.105427357601002e-15 * bound then
        return det
    end
    P.stats.incircle = P.stats.incircle + 1
    ax, ay = diff(a[1], d[1]), diff(a[2], d[2])
    bx, by = diff(b[1], d[1]), diff(b[2], d[2])
    cx, cy = diff(c[1], d[1]), diff(c[2], d[2])
    al = sum(product(ax, ax), product(ay, ay))
    bl = sum(product(bx, bx), product(by, by))
    cl = sum(product(cx, cx), product(cy, cy))
    return leading(
        sum(
            sum(product(al, determinant(bx, by, cx, cy)), product(bl, determinant(cx, cy, ax, ay))),
            product(cl, determinant(ax, ay, bx, by))
        )
    )
end
--- Test exact collinearity and inclusive coordinate bounds.
-- @tparam table p Query point.
-- @tparam table a Segment start.
-- @tparam table b Segment end.
-- @treturn boolean Whether p lies on the closed segment.
function P.on_segment(p, a, b)
    return P.orient(a, b, p) == 0
        and p[1] >= math.min(a[1], b[1])
        and p[1] <= math.max(a[1], b[1])
        and p[2] >= math.min(a[2], b[2])
        and p[2] <= math.max(a[2], b[2])
end
return P
