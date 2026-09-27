-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Exact integration of the piecewise-constant P1 field along a polyline.
-- @module luafemm-integrals
local L = {}

--- Integrate H.dl over a segment, splitting at every crossed triangle edge.
-- Coordinates are model units; result is amperes. At an edge, the first
-- containing triangle determines H, matching sample(). Missing coverage errors.
function L.segment(m, a, b)
    local unit = m.options.unit
    local x, y = a[1] * unit, a[2] * unit
    local dx, dy = (b[1] - a[1]) * unit, (b[2] - a[2]) * unit
    if dx == 0 and dy == 0 then
        return 0
    end
    local intervals, cuts = {}, { 0, 1 }
    for index, e in ipairs(m.triangles) do
        local p = m.nodes[e.ids[1]]
        local v = {
            0,
            e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2]),
            e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2]),
        }
        v[1] = 1 - v[2] - v[3]
        local lo, hi = 0, 1
        for j = 1, 3 do
            local slope = e.gx[j] * dx + e.gy[j] * dy
            -- Match the dimensionless barycentric tolerance in sample(). A
            -- nearly parallel segment on a mesh edge must not acquire a large
            -- artificial parameter cut from dividing two roundoff-sized values.
            if math.abs(slope) < 1e-10 then
                if v[j] < -1e-10 then
                    hi = -1
                end
            elseif slope > 0 then
                lo = math.max(lo, -v[j] / slope)
            else
                hi = math.min(hi, -v[j] / slope)
            end
        end
        if hi - lo > 1e-12 then
            intervals[#intervals + 1] = { lo, hi, index }
            cuts[#cuts + 1], cuts[#cuts + 2] = lo, hi
        end
    end
    table.sort(cuts)
    local total = 0
    for i = 2, #cuts do
        local lo, hi = cuts[i - 1], cuts[i]
        if hi - lo > 1e-11 then
            local mid, element = (lo + hi) / 2
            for _, interval in ipairs(intervals) do
                if mid >= interval[1] - 1e-12 and mid <= interval[2] + 1e-12 then
                    element = m.triangles[interval[3]]
                    break
                end
            end
            assert(
                element,
                string.format(
                    "luafemm integral: path leaves the computational domain near (%g,%g), interval %.12g:%.12g",
                    (x + mid * dx) / unit,
                    (y + mid * dy) / unit,
                    lo,
                    hi
                )
            )
            total = total + (hi - lo) * (element.hx * dx + element.hy * dy)
        end
    end
    return total
end
return L
