-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Shared compound-region membership for material classification and meshing.
-- Contours are individually simple. Distinct contours may overlap: the
-- mesher splits their intersections and this predicate follows the PGF fill rule.
-- @module luafemm-geometry
local G = {}
local orient = require("luafemm-predicates").orient

--- Test membership using half-open ray crossings and robust orientation signs.
-- Boundary points follow the half-open convention; element centroids are
-- classified after constraint recovery, so none should straddle an interface.
-- @tparam table region contours and fill_rule.
-- @tparam number x Coordinate in the same units as the contours.
-- @tparam number y Coordinate in the same units as the contours.
-- @treturn boolean Whether the selected fill rule includes this point.
function G.contains(region, x, y)
    local winding, point = 0, { x, y }
    for _, contour in ipairs(region.contours) do
        for i, a in ipairs(contour) do
            local b = contour[i % #contour + 1]
            if a[2] <= y and b[2] > y and orient(a, b, point) > 0 then
                winding = winding + 1
            elseif a[2] > y and b[2] <= y and orient(a, b, point) < 0 then
                winding = winding - 1
            end
        end
    end
    if region.fill_rule == "even odd" then
        return winding % 2 ~= 0
    end
    return winding ~= 0
end

return G
