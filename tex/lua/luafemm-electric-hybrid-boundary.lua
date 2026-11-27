-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Analytic reference flux, natural boundary loads and total-field reactions.
-- The active dielectric's outward normal points into a solid conductor. Thus
-- Ku-b is conductor charge, including the reference-source flux, not merely
-- the correction's charge. Internal thin electrodes have two cancelling
-- reference traces. All assembly quantities are per extrusion length in 2D,
-- or full-revolution quantities in axisymmetry, before the reporting depth.
local B = {}
local E = require("luafemm-electrostatic")
local metric = require("luafemm-electric-metric")
local Q = require("luafemm-electric-quadrature")
local function geometry(m, e)
    local f = m.fe
    local a, c = f.nodes[e[1]], f.nodes[e[2]]
    local dx, dy = c[1] - a[1], c[2] - a[2]
    local length = math.sqrt(dx * dx + dy * dy)
    local nx, ny = dy / length, -dx / length
    local active, n = nil, 0
    for _, id in ipairs(e.triangles) do
        if not f.triangles[id].electrode then
            active, n = id, n + 1
        end
    end
    if n == 1 then
        for _, id in ipairs(f.triangles[active].ids) do
            if id ~= e[1] and id ~= e[2] then
                local p = f.nodes[id]
                if nx * (p[1] - a[1]) + ny * (p[2] - a[2]) > 0 then
                    nx, ny = -nx, -ny
                end
            end
        end
    end
    return function(t)
        local x, y = a[1] + t * dx, a[2] + t * dy
        local weight = length * m.unit * (f.axisymmetric and 2 * math.pi * x * m.unit or 1)
        return x, y, weight
    end,
        n == 1,
        nx,
        ny
end
function B.assemble(m, K, b, base, budget)
    local f, h = m.fe, m.hybrid
    local integrate, count =
        Q.new(h.quadrature_tolerance, h.quadrature_depth, budget, "source boundary")
    for _, e in ipairs(f.edges) do
        local a, c = f.nodes[e[1]], f.nodes[e[2]]
        local edge, exterior, nx, ny = geometry(m, e)
        local coefficient = e.condition == "mixed" and e.coefficient or 0
        local outward = e.outer
                and (e.condition == "displacement" or e.condition == "mixed")
                and e.displacement
            or 0
        if coefficient > 0 then
            local mass = metric.edge_mass(m, a, c)
            for i = 1, 2 do
                for j = 1, 2 do
                    K[e[i]][e[j]] = (K[e[i]][e[j]] or 0) + coefficient * mass[i][j] / E.epsilon0
                end
            end
        end
        local measure, fractions = metric.edge(m, a, c)
        for i = 1, 2 do
            b[e[i]] = b[e[i]] - outward * measure * fractions[i] / E.epsilon0
        end
        if #m.sources > 0 and not e.axis and (exterior or coefficient > 0) then
            local load = integrate(function(t)
                local x, y, w = edge(t)
                local s = assert(base(x, y))
                local flux = (exterior and h.epsilon_r * (s.Ex * nx + s.Ey * ny) or 0)
                    - coefficient * s.V / E.epsilon0
                return { w * flux * (1 - t), w * flux * t }
            end)
            for i = 1, 2 do
                b[e[i]] = b[e[i]] + load[i]
            end
        end
    end
    return count()
end
function B.reactions(m, K, b, u, base, budget)
    local f = m.fe
    local charges, voltages = {}, {}
    charges["outer boundary"] = 0
    for _, name in ipairs(f.electrodes) do
        charges[name] = 0
    end
    local all, periodic = 0, 0
    local variable = {}
    local function reaction(i)
        local r = -b[i]
        for j, k in pairs(K[i]) do
            r = r + k * u[j]
        end
        return r * E.epsilon0
    end
    for i, name in pairs(f.owners) do
        local p = f.nodes[i]
        local v = u[i] + assert(base(p[1], p[2])).V
        if voltages[name] and math.abs(v - voltages[name]) > 1e-9 * math.max(1, math.abs(v)) then
            variable[name] = true
        end
        voltages[name] = v
        local r = reaction(i)
        charges[name], all = charges[name] + r, all + r
    end
    for name in pairs(variable) do
        voltages[name] = nil
    end
    for i in pairs(f.periodic_nodes) do
        if not f.owners[i] then
            periodic = periodic + reaction(i)
        end
    end
    local flux = 0
    -- Reuse the source boundary tolerance; the same nonlinear Vs trace enters
    -- the mixed load and the reported exterior flux. No total P1 interpolation.
    local h = m.hybrid
    local integrate, count =
        Q.new(h.quadrature_tolerance, h.quadrature_depth, budget, "source flux")
    for _, e in ipairs(f.edges) do
        if e.outer and (e.condition == "displacement" or e.condition == "mixed") then
            local edge = geometry(m, e)
            local coefficient = e.condition == "mixed" and e.coefficient or 0
            local value = integrate(function(t)
                local x, y, w = edge(t)
                local v = assert(base(x, y)).V + (1 - t) * u[e[1]] + t * u[e[2]]
                return { w * (e.displacement + coefficient * v) }
            end)
            flux = flux + value[1]
        end
    end
    local source = 0
    for _, s in ipairs(m.sources) do
        source = source + s.q
    end
    f.charges, f.voltages, f.source_charge = charges, voltages, source
    f.periodic_flux = -periodic
    f.balance = source + all + periodic - flux
    return count()
end
return B
