-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Exact measures for planar and axisymmetric scalar Poisson P1 elements.
-- Coordinates are model units. We integrate r*N_i exactly for constant volume
-- and sheet densities; replacing this with centroid lumping changes the load.
-- The radial factor has no singularity at the regular scalar-potential axis.
-- @module luafemm-electric-metric
local M = {}
function M.triangle(m, t)
    local sum = 0
    for _, id in ipairs(t.ids) do
        sum = sum + m.fe.nodes[id][1]
    end
    if not m.fe.axisymmetric then
        return 1, { 1 / 3, 1 / 3, 1 / 3 }
    end
    local fractions = {}
    for i, id in ipairs(t.ids) do
        fractions[i] = (sum + m.fe.nodes[id][1]) / (4 * sum)
    end
    return 2 * math.pi * m.unit * sum / 3, fractions
end
function M.edge(m, a, b)
    local length = math.sqrt((a[1] - b[1]) ^ 2 + (a[2] - b[2]) ^ 2) * m.unit
    if not m.fe.axisymmetric then
        return length, { 0.5, 0.5 }
    end
    local sum = a[1] + b[1]
    if sum == 0 then
        return 0, { 0.5, 0.5 }
    end
    return length * math.pi * m.unit * sum,
        { (2 * a[1] + b[1]) / (3 * sum), (a[1] + 2 * b[1]) / (3 * sum) }
end
--- Exact integral of N_i*N_j, including the linear radial weight.
function M.edge_mass(m, a, b)
    local measure = M.edge(m, a, b)
    if not m.fe.axisymmetric then
        return { { measure / 3, measure / 6 }, { measure / 6, measure / 3 } }
    end
    local sum = a[1] + b[1]
    if sum == 0 then
        return { { 0, 0 }, { 0, 0 } }
    end
    return {
        { measure * (3 * a[1] + b[1]) / (6 * sum), measure / 6 },
        { measure / 6, measure * (a[1] + 3 * b[1]) / (6 * sum) },
    }
end
return M
