-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- luafemm-path implementation.
-- @module luafemm-path
-- Read evaluated PGF soft paths, including cubic Beziers and rectangles.
local P = {}
local function check(v, s)
    assert(v, "luafemm path: " .. s)
end
local function dim(s)
    return assert(tonumber(s:match("([%d%.%+%-eE]+)pt")), "invalid PGF dimension " .. s)
end
--- Compose the PGF canvas transform with the model coordinate basis.
-- @tparam string transform Six braced values returned by pgfgettransform.
-- @tparam table basis Four PGF dimensions (xx, xy, yx, yy).
-- @treturn table Affine map from model coordinates to TeX points.
function P.frame(transform, basis)
    local values = {}
    for v in transform:gmatch("{([^}]+)}") do
        values[#values + 1] = v
    end
    local a, b, c, d =
        tonumber(values[1]), tonumber(values[2]), tonumber(values[3]), tonumber(values[4])
    local xx, xy, yx, yy = dim(basis[1]), dim(basis[2]), dim(basis[3]), dim(basis[4])
    local f = {
        a * xx + c * xy,
        b * xx + d * xy,
        a * yx + c * yy,
        b * yx + d * yy,
        dim(values[5]),
        dim(values[6]),
    }
    check(math.abs(f[1] * f[4] - f[2] * f[3]) > 1e-14, "singular picture transform")
    return f
end
-- Internal parser shared by physical compound paths and single profiles.
local function read_contours(s, frame, tolerance, allow_open)
    check(tolerance > 0, "curve tolerance must be positive")
    local f = frame
    local det = f[1] * f[4] - f[2] * f[3]
    local function point(x, y, vector)
        x, y = dim(x), dim(y)
        if not vector then
            x, y = x - f[5], y - f[6]
        end
        return { (f[4] * x - f[3] * y) / det, (-f[2] * x + f[1] * y) / det }
    end
    local points = {}
    local closed = false
    local drawn = false
    local ca, cb, corner
    local contours, closures = {}, {}
    local function append(p)
        if #points > 0 then
            local a = points[#points]
            if (p[1] - a[1]) ^ 2 + (p[2] - a[2]) ^ 2 < 1e-16 then
                return
            end
        end
        points[#points + 1] = p
    end
    local function finish()
        -- PGF emits isolated moves after circles/rectangles. They are not
        -- additional contours. Drawing commands, however, require closure.
        if drawn or closed then
            check(
                closed or allow_open,
                "physical regions require a closed TikZ path (use -- cycle)"
            )
            if closed and #points > 1 then
                local a, b = points[1], points[#points]
                if (a[1] - b[1]) ^ 2 + (a[2] - b[2]) ^ 2 < 1e-12 then
                    table.remove(points)
                end
            end
            contours[#contours + 1] = points
            closures[#contours] = closed
        end
        points, closed, ca, cb, corner = {}, false, nil, nil, nil
        drawn = false
    end
    local function mid(a, b)
        return { (a[1] + b[1]) / 2, (a[2] + b[2]) / 2 }
    end
    local function bezier(a, b, c, d, depth)
        local dx, dy = d[1] - a[1], d[2] - a[2]
        local l = dx * dx + dy * dy
        local function distance(p)
            local t = l > 0
                    and math.max(0, math.min(1, ((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / l))
                or 0
            return (p[1] - a[1] - t * dx) ^ 2 + (p[2] - a[2] - t * dy) ^ 2
        end
        if math.max(distance(b), distance(c)) <= tolerance * tolerance then
            append(d)
            return
        end
        check(depth < 20, "curve subdivision limit reached")
        local ab, bc, cd = mid(a, b), mid(b, c), mid(c, d)
        local abc, bcd = mid(ab, bc), mid(bc, cd)
        local center = mid(abc, bcd)
        bezier(a, ab, abc, center, depth + 1)
        bezier(center, bcd, cd, d, depth + 1)
    end
    local tokens = 0
    for kind, x, y in s:gmatch("\\pgfsyssoftpath@(%a+)token%s*{([^}]+)}%s*{([^}]+)}") do
        tokens = tokens + 1
        if kind == "moveto" then
            check(
                not allow_open or closed or not drawn,
                "one connected contour per path; use separate profiles"
            )
            finish()
            points = { point(x, y) }
        elseif kind == "lineto" then
            check(not closed and #points > 0, "start a new contour with a move after cycle")
            drawn = true
            append(point(x, y))
        elseif kind == "curvetosupporta" then
            ca = point(x, y)
        elseif kind == "curvetosupportb" then
            cb = point(x, y)
        elseif kind == "curveto" then
            check(not closed and #points > 0 and ca and cb, "invalid cubic path")
            drawn = true
            bezier(points[#points], ca, cb, point(x, y), 0)
            ca, cb = nil, nil
        elseif kind == "rectcorner" then
            finish()
            corner = point(x, y)
        elseif kind == "rectsize" then
            check(corner, "rectangle size without corner")
            -- PGF uses rect tokens only for axis-aligned rectangles in canvas
            -- coordinates; invert the two separate vectors for the model frame.
            local u, v = point(x, "0pt", true), point("0pt", y, true)
            append(corner)
            append({ corner[1] + u[1], corner[2] + u[2] })
            append({ corner[1] + u[1] + v[1], corner[2] + u[2] + v[2] })
            append({ corner[1] + v[1], corner[2] + v[2] })
            closed = true
        elseif kind == "closepath" then
            check(#points > 0, "closing an empty path")
            closed = true
        else
            error("luafemm path: unsupported PGF token " .. kind)
        end
    end
    finish()
    check(tokens > 0 and #contours > 0, "path contains no drawable contour")
    return contours, closures
end

--- Flatten all closed subpaths of a physical region into model coordinates.
-- Each cubic uses de Casteljau subdivision with a chord-distance bound.
-- @tparam string s Evaluated PGF soft-path token list.
-- @tparam table frame Model-to-canvas affine map.
-- @tparam number tolerance Positive chord tolerance in model units.
-- @treturn table Contours, each an array of vertices without closure repeat.
function P.read_contours(s, frame, tolerance)
    local contours = read_contours(s, frame, tolerance, false)
    return contours
end

--- Flatten one connected contour, retaining the single-path profile contract.
-- @tparam string s Evaluated PGF soft-path token list.
-- @tparam table frame Model-to-canvas affine map.
-- @tparam number tolerance Positive chord tolerance in model units.
-- @tparam[opt=false] boolean allow_open Accept open measurement paths.
-- @return Vertices and closure flag. Closed profiles repeat their first vertex.
function P.read(s, frame, tolerance, allow_open)
    local contours, closures = read_contours(s, frame, tolerance, allow_open)
    check(#contours == 1, "one connected contour per path; use separate profiles")
    local points, closed = contours[1], closures[1]
    if allow_open and closed and #points > 1 then
        points[#points + 1] = { points[1][1], points[1][2] }
    end
    return points, closed
end
return P
