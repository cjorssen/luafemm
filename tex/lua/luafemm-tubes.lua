-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Physics-neutral faceted tubes around owned polygonal centre lines.
-- Radius and vertices use scene units. Cross sections use rotation-minimising
-- frames; closed curves distribute residual twist to close the seam. This is
-- display geometry, not a swept-volume field quadrature or solid modeller.
local M = {}
local function sub(a, b)
    return { a[1] - b[1], a[2] - b[2], a[3] - b[3] }
end
local function dot(a, b)
    return a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
end
local function cross(a, b)
    return { a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1] }
end
local function unit(a)
    local n = math.sqrt(dot(a, a))
    assert(n > 1e-14, "luafemm 3D tubes: duplicate points or reversing tangent")
    return { a[1] / n, a[2] / n, a[3] / n }
end
local function rotate(a, axis, c, s)
    local b = cross(axis, a)
    local d = dot(axis, a) * (1 - c)
    return {
        a[1] * c + b[1] * s + axis[1] * d,
        a[2] * c + b[2] * s + axis[2] * d,
        a[3] * c + b[3] * s + axis[3] * d,
    }
end
--- Return planar triangles and end caps; call reserve before allocating faces.
function M.faces(points, radius, sides, reserve)
    local n = #points
    local closed = dot(sub(points[1], points[n]), sub(points[1], points[n])) < 1e-24
    local tangents, frames, rings = {}, {}, {}
    for i, p in ipairs(points) do
        local previous = points[i - 1] or (closed and points[n - 1]) or p
        local following = points[i + 1] or (closed and points[2]) or p
        tangents[i] = unit(sub(following, previous))
    end
    local t = tangents[1]
    frames[1] = unit(cross(t, math.abs(t[3]) < 0.9 and { 0, 0, 1 } or { 0, 1, 0 }))
    for i = 2, n do
        local a, b = tangents[i - 1], tangents[i]
        local axis = cross(a, b)
        local s = math.sqrt(dot(axis, axis))
        local c = math.max(-1, math.min(1, dot(a, b)))
        assert(c > -0.999999, "luafemm 3D tubes: reversing tangent; split the path")
        frames[i] = s < 1e-12 and frames[i - 1] or rotate(frames[i - 1], unit(axis), c, s)
    end
    if closed then
        local angle = math.atan(dot(t, cross(frames[n], frames[1])), dot(frames[n], frames[1]))
        for i = 2, n do
            local a = angle * (i - 1) / (n - 1)
            frames[i] = rotate(frames[i], tangents[i], math.cos(a), math.sin(a))
        end
    end
    for i, p in ipairs(points) do
        local u, v = frames[i], cross(tangents[i], frames[i])
        local ring = {}
        for j = 1, sides do
            local a = 2 * math.pi * (j - 1) / sides
            local c, s = radius * math.cos(a), radius * math.sin(a)
            ring[j] = {
                p[1] + c * u[1] + s * v[1],
                p[2] + c * u[2] + s * v[2],
                p[3] + c * u[3] + s * v[3],
            }
        end
        rings[i] = ring
    end
    if closed then
        rings[n] = rings[1]
    end
    local faces = {}
    local function add(p)
        reserve()
        faces[#faces + 1] = p
    end
    for i = 2, n do
        for j = 1, sides do
            local k = j % sides + 1
            add({ rings[i - 1][j], rings[i - 1][k], rings[i][j] })
            add({ rings[i - 1][k], rings[i][k], rings[i][j] })
        end
    end
    if not closed then
        local cap = {}
        for j = sides, 1, -1 do
            cap[#cap + 1] = rings[1][j]
        end
        add(cap)
        add(rings[n])
    end
    return faces
end
return M
