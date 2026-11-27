-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Owned prescribed-current geometry in model units, independent of projection.
-- Curves return position and derivative with respect to a dimensionless parameter.
-- Only registered records are persisted; evaluator closures are rebuilt locally.
local G = {}
function G.check(v, s)
    assert(v, "luafemm Biot-Savart: " .. s)
end
function G.number(v, name, default)
    if v == nil or v == "" then
        v = default
    end
    v = tonumber(v)
    G.check(v and v == v and math.abs(v) < math.huge, "finite " .. name .. " required")
    return v
end
function G.vector(v, name)
    if type(v) == "string" then
        local a, b, c = v:match("^%s*([^/]+)/([^/]+)/([^/]+)%s*$")
        G.check(a, name .. " must be x/y/z")
        v = { a, b, c }
    end
    G.check(type(v) == "table" and #v == 3, name .. " must have three components")
    local count = 0
    for k in pairs(v) do
        count = count + 1
        G.check(k == 1 or k == 2 or k == 3, "invalid " .. name .. " component index")
    end
    G.check(count == 3, "sparse " .. name)
    return { G.number(v[1], name), G.number(v[2], name), G.number(v[3], name) }
end
function G.add(a, b)
    return { a[1] + b[1], a[2] + b[2], a[3] + b[3] }
end
function G.sub(a, b)
    return { a[1] - b[1], a[2] - b[2], a[3] - b[3] }
end
function G.mul(a, s)
    return { a[1] * s, a[2] * s, a[3] * s }
end
function G.dot(a, b)
    return a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
end
function G.cross(a, b)
    return { a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1] }
end
function G.norm(a)
    return math.sqrt(G.dot(a, a))
end
local add, sub, mul = G.add, G.sub, G.mul
local pi = math.pi
local function rotate(p, o)
    local x, y, z = p[1], p[2], p[3]
    for k, angle in ipairs(o.rotation) do
        angle = angle % 360
        local c, s = math.cos(angle * pi / 180), math.sin(angle * pi / 180)
        if k == 1 then
            y, z = c * y - s * z, s * y + c * z
        elseif k == 2 then
            x, z = c * x + s * z, -s * x + c * z
        else
            x, y = c * x - s * y, s * x + c * y
        end
    end
    return { x, y, z }
end
function G.place(r, p)
    return add(rotate(p, r), r.center)
end
--- Validate and normalize one declaration, without retaining caller-owned tables.
function G.record(o)
    local n = G.number
    local r = {
        kind = o.kind or "loop",
        current = n(o.current, "current"),
        wire_radius = n(o.wire_radius, "wire radius", 0),
        center = G.vector(o.center or "0/0/0", "center"),
        rotation = {
            n(o.rotate_x, "rotate x", 0),
            n(o.rotate_y, "rotate y", 0),
            n(o.rotate_z, "rotate z", 0),
        },
        closure = o.closure or "closed",
        winding = o.winding or "turns",
    }
    G.check(r.wire_radius >= 0, "wire radius must be nonnegative")
    G.check(
        r.closure == "closed" or r.closure == "contribution" or r.closure == "axis",
        "unknown closure"
    )
    G.check(r.winding == "turns" or r.winding == "helix", "unknown winding")
    if r.kind == "wire" then
        r.extent = o.extent or "finite"
        G.check(r.extent == "finite" or r.extent == "infinite", "unknown wire extent")
        r.points = {}
        local points = o.points
            or { G.vector(o.from or "0/0/-1", "from"), G.vector(o.to or "0/0/1", "to") }
        G.check(
            type(points) == "table" and #points >= 2 and #points <= 100000,
            "wire needs 2 to 100000 points"
        )
        local count = 0
        for k in pairs(points) do
            count = count + 1
            G.check(
                type(k) == "number" and k % 1 == 0 and k >= 1 and k <= #points,
                "invalid wire point index"
            )
        end
        G.check(count == #points, "sparse wire points")
        for _, p in ipairs(points) do
            r.points[#r.points + 1] = G.vector(p, "wire point")
        end
        for i = 2, #r.points do
            G.check(G.norm(sub(r.points[i], r.points[i - 1])) > 1e-12, "zero-length wire segment")
        end
        if r.extent == "infinite" then
            G.check(#r.points == 2, "infinite wire needs two points")
        else
            G.check(r.closure ~= "axis", "axis return is only available for a solenoid helix")
            if r.closure == "closed" then
                G.check(
                    G.norm(sub(r.points[1], r.points[#r.points])) < 1e-9,
                    "close the current path or declare closure=contribution"
                )
                r.points[#r.points] = G.vector(r.points[1], "closure")
            end
        end
    else
        G.check(
            r.kind == "loop" or r.kind == "solenoid" or r.kind == "toroid",
            "unknown source kind"
        )
        G.check(not o.extent or o.extent == "finite", "infinite extent requires kind=wire")
        G.check(
            r.closure ~= "axis" or r.kind == "solenoid" and r.winding == "helix",
            "axis return requires a solenoid helix"
        )
        r.radius = n(o.radius, "radius", 1)
        r.turns = n(o.turns, "turns", 10)
        G.check(r.radius > 0, "radius must be positive")
        G.check(
            r.turns % 1 == 0 and r.turns >= 1 and r.turns <= 1000,
            "turns must be an integer from 1 to 1000"
        )
        r.phase = (n(o.phase, "phase", 0) % 360) * pi / 180
        r.handedness = o.handedness or "right"
        G.check(r.handedness == "right" or r.handedness == "left", "unknown handedness")
        if r.kind == "solenoid" then
            r.length = n(o.length, "length", 3)
            G.check(r.length > 0, "length must be positive")
            if r.winding == "helix" then
                G.check(r.closure ~= "closed", "solenoid helix needs closure=axis or contribution")
            end
        elseif r.kind == "toroid" then
            r.major_radius = n(o.major_radius, "major radius", 2)
            G.check(r.major_radius > r.radius, "major radius must exceed the winding radius")
        end
        G.check(r.wire_radius < r.radius, "wire radius must be smaller than the bend radius")
    end
    return r
end
--- Compile validated records into parametric curves and exact straight segments.
function G.compile(r)
    local curves = {}
    local function transform(p)
        return add(rotate(p, r), r.center)
    end
    local function line(a, b, infinite, part)
        a, b = transform(a), transform(b)
        local d = sub(b, a)
        curves[#curves + 1] = {
            line = true,
            part = part or "winding",
            infinite = infinite,
            a = a,
            b = b,
            speed = G.norm(d),
            panels = 1,
            evaluate = function(t)
                return add(a, mul(d, t)), d
            end,
        }
    end
    local function curve(fn, panels, speed)
        curves[#curves + 1] = {
            panels = panels,
            speed = speed,
            evaluate = function(t)
                local p, d = fn(t)
                return transform(p), rotate(d, r)
            end,
        }
    end
    local function circle(center, u, v, phase)
        curve(function(t)
            local a = 2 * pi * t + (phase or 0)
            local c, s = math.cos(a), math.sin(a)
            return add(center, add(mul(u, c), mul(v, s))), mul(add(mul(u, -s), mul(v, c)), 2 * pi)
        end, 16, 2 * pi * r.radius)
    end
    if r.kind == "wire" then
        for i = 2, #r.points do
            line(r.points[i - 1], r.points[i], r.extent == "infinite")
        end
    elseif r.kind == "loop" then
        circle({ 0, 0, 0 }, { r.radius, 0, 0 }, { 0, r.radius, 0 }, r.phase)
    elseif r.kind == "solenoid" and r.winding == "turns" then
        for i = 1, r.turns do
            local z = r.turns == 1 and 0 or r.length * ((i - 1) / (r.turns - 1) - 0.5)
            circle({ 0, 0, z }, { r.radius, 0, 0 }, { 0, r.radius, 0 }, r.phase)
        end
    elseif r.kind == "toroid" and r.winding == "turns" then
        for i = 1, r.turns do
            local a = 2 * pi * (i - 1) / r.turns + r.phase
            local c, s = math.cos(a), math.sin(a)
            circle(
                { r.major_radius * c, r.major_radius * s, 0 },
                { r.radius * c, r.radius * s, 0 },
                { 0, 0, r.radius }
            )
        end
    elseif r.kind == "solenoid" then
        local w = 2 * pi * r.turns * (r.handedness == "right" and 1 or -1)
        local function helix(t)
            local a = w * t + r.phase
            local c, s = math.cos(a), math.sin(a)
            return { r.radius * c, r.radius * s, r.length * (t - 0.5) }, {
                -w * r.radius * s,
                w * r.radius * c,
                r.length,
            }
        end
        curve(helix, 16 * r.turns, math.sqrt((w * r.radius) ^ 2 + r.length ^ 2))
        if r.closure == "axis" then
            local first, last = helix(0), helix(1)
            line(last, { 0, 0, r.length / 2 }, nil, "leads")
            line({ 0, 0, r.length / 2 }, { 0, 0, -r.length / 2 }, nil, "return")
            line({ 0, 0, -r.length / 2 }, first, nil, "leads")
        end
    else
        local w = 2 * pi * r.turns * (r.handedness == "right" and 1 or -1)
        curve(function(t)
            local a, b = 2 * pi * t + r.phase, w * t
            local c, s, cb, sb = math.cos(a), math.sin(a), math.cos(b), math.sin(b)
            local h = r.major_radius + r.radius * cb
            return { h * c, h * s, r.radius * sb }, {
                -w * r.radius * sb * c - 2 * pi * h * s,
                -w * r.radius * sb * s + 2 * pi * h * c,
                w * r.radius * cb,
            }
        end, 16 * r.turns, 2 * pi * (r.major_radius + r.radius) + math.abs(w) * r.radius)
    end
    return curves
end
return G
