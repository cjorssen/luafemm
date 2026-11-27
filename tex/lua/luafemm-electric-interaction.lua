-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Finite source charging work relative to removed homogeneous self-kernels.
-- The linear boundary/material response is charged with factor 1/2; an existing
-- applied field contributes with factor 1. This is an interaction potential at
-- fixed external boundary data, not the divergent total field energy. Planar
-- values depend on the declared logarithmic reference by a position-independent
-- constant. No drawing mask enters the definition.
local M = {}
function M.calculate(m)
    assert(m.hybrid, "luafemm electric interaction: source correction model required")
    local H = require("luafemm-electric-hybrid")
    local F = require("luafemm-electric-fe")
    H.solve(m)
    local h = m.hybrid
    if h.interaction_energy ~= nil then
        return h.interaction_energy
    end
    local baseline = require("luafemm-physical").copy(m)
    baseline.sources = {}
    baseline.frozen = false
    baseline.hybrid.solution = nil
    H.solve(baseline)
    local energy = 0
    for i, s in ipairs(m.sources) do
        local id, w = F.locate(m, s.x, s.y)
        assert(id, "luafemm electric interaction: source outside dielectric")
        local u, external = 0, 0
        for j, node in ipairs(m.fe.triangles[id].ids) do
            u = u + w[j] * h.solution[node]
            external = external + w[j] * baseline.hybrid.solution[node]
        end
        energy = energy + 0.5 * s.q * (u + external)
        -- Count each direct source pair once; never evaluate its self-kernel.
        for j = 1, i - 1 do
            local t = m.sources[j]
            local d = math.sqrt((s.x - t.x) ^ 2 + (s.y - t.y) ^ 2) * m.unit
            local epsilon = require("luafemm-electrostatic").epsilon0 * h.epsilon_r
            local v = m.fe.axisymmetric and t.q / (4 * math.pi * epsilon * d)
                or -t.q * math.log(d / (h.reference_radius * m.unit)) / (2 * math.pi * epsilon)
            energy = energy + s.q * v
        end
    end
    assert(
        energy == energy and math.abs(energy) < math.huge,
        "luafemm electric interaction: energy overflow"
    )
    h.interaction_energy = energy
    return energy
end
return M
