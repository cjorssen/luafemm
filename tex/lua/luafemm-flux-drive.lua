-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Static prescribed-linkage constraints with unknown winding currents.
-- Newton eliminates the SPD magnetic tangent before solving the small Schur
-- system. The augmented indefinite problem is never passed to conjugate gradient.
-- Linkage is a linear, quadrature-consistent functional of A; nonlinear material
-- tangents and permanent-magnet offsets remain in the ordinary field assembler.
local D = {}
local function check(v, s)
    assert(v, "luafemm flux drive: " .. s)
end
local function dot(a, b)
    local s = 0
    for k, v in pairs(a) do
        s = s + v * b[k]
    end
    return s
end
local function norm(a)
    local scale, s = 0, 1
    for _, v in pairs(a) do
        v = math.abs(v)
        if v > 0 then
            if v > scale then
                s, scale = 1 + s * (scale / v) ^ 2, v
            else
                s = s + (v / scale) ^ 2
            end
        end
    end
    return scale * math.sqrt(s)
end
function D.active(m)
    if m.options.frequency > 0 then
        return false
    end
    for _, c in pairs(m.circuits or {}) do
        if c.drive == "linkage" or c.drive == "flux" then
            return true
        end
    end
    return false
end
local function prepare(m)
    local ports = {}
    for id, c in pairs(m.circuits or {}) do
        if c.drive == "linkage" or c.drive == "flux" then
            check(c.type == "series", "use series windings for a linkage constraint")
            ports[#ports + 1] = {
                name = id,
                target = c.drive == "flux" and c.flux * c.flux_turns or c.linkage,
                load = {},
                link = {},
                density = {},
            }
        end
    end
    table.sort(ports, function(a, b)
        return a.name < b.name
    end)
    check(#ports <= 32, "at most 32 controlled windings")
    local axi = m.options.problem_type == "axisymmetric"
    check(axi or m.options.problem_type == "planar", "planar or axisymmetric model required")
    check(m.options.field_model == "open", "flux drive requires field model=open")
    for _, p in ipairs(ports) do
        local c = m.static_circuits[p.name]
        for _, g in pairs(c.windings) do
            check(axi or math.abs(g.turns) < 1e-10, "closed signed-turn windings are required")
        end
        for k, e in ipairs(m.triangles) do
            local id = e.face and m.topology.faces[e.face].region
            local r = id and m.regions[id]
            if r and r.circuit == p.name then
                local density = r.turns / c.areas[id]
                p.density[k] = density
                for _, q in ipairs(require("luafemm-harmonic").quadrature(m, e)) do
                    for j, n in ipairs(e.ids) do
                        local w = q.weight * q.n[j] * density
                        p.link[n] = (p.link[n] or 0)
                            + w * (axi and 2 * math.pi or m.options.depth * m.options.unit)
                        local node = m.nodes[n]
                        if node.dof then
                            p.load[node.dof] = (p.load[node.dof] or 0)
                                + 4 * math.pi * 1e-7 * w * (node.dof_sign or 1)
                        end
                    end
                end
            end
        end
    end
    return ports
end
--- Restore physical source data from persisted unknown currents.
function D.recover(m)
    if not D.active(m) then
        return
    end
    local ports = prepare(m)
    local state = m.circuit_state and m.circuit_state.dc
    check(state and #state == #ports, "missing static flux-driven current state")
    for j, p in ipairs(ports) do
        local c = m.static_circuits[p.name]
        c.current = state[j]
        c.voltage = c.resistance and c.resistance * state[j]
        for _, g in pairs(c.windings) do
            g.current = state[j]
        end
        for k, density in pairs(p.density) do
            local e = m.triangles[k]
            e.current = density * state[j]
            e.meridional_current = e.current * e.area
            for _, q in ipairs(e.quadrature or {}) do
                q.current = e.current
            end
        end
    end
end
function D.solve(m, assemble)
    local ports = prepare(m)
    local M, L = require("luafemm"), require("luafemm-linear")
    local A, I = m.A, {}
    for j = 1, #ports do
        I[j] = 0
    end
    local started, steps, stats = os.clock(), 0, nil
    local field_scale = norm(assemble(m, A, false))
    for _, p in ipairs(ports) do
        p.load_norm = norm(p.load)
    end
    local function residual(a, current, tangent)
        local r, k = assemble(m, a, tangent)
        for j, p in ipairs(ports) do
            for row, v in pairs(p.load) do
                r[row] = r[row] - v * current[j]
            end
        end
        local g, cr = {}, 0
        for j, p in ipairs(ports) do
            g[j] = dot(p.link, a) - p.target
            local scale = p.target ~= 0 and math.abs(p.target) or norm(p.link) * norm(a)
            cr = math.max(cr, scale == 0 and 0 or math.abs(g[j]) / scale)
        end
        local scale = field_scale
        for j, p in ipairs(ports) do
            scale = scale + p.load_norm * math.abs(current[j])
        end
        local residual_norm = norm(r)
        return r,
            k,
            g,
            residual_norm == 0 and 0 or residual_norm / math.max(scale, residual_norm),
            cr
    end
    local function free_link(p, x)
        local v = 0
        for n, w in pairs(p.link) do
            local node = m.nodes[n]
            if node.dof then
                v = v + w * (node.dof_sign or 1) * x[node.dof]
            end
        end
        return v
    end
    for iteration = 0, m.options.max_newton do
        local r, k, g, fr, cr = residual(A, I, true)
        if fr <= m.options.tolerance and cr <= m.options.tolerance then
            m.A = A
            m.circuit_state = { dc = I }
            m.stats = {
                newton = iteration,
                cg = steps,
                gmres = 0,
                linear = stats,
                residual = fr,
                constraint_residual = cr,
                seconds = os.clock() - started,
                nodes = #m.nodes,
                triangles = #m.triangles,
            }
            D.recover(m)
            M.update_fields(m)
            require("luafemm-cache").store(m)
            return m
        end
        check(iteration < m.options.max_newton, "constraint Newton iteration limit")
        local function solve(rhs)
            local x, it, t = L.real(m, k, rhs, 1e-11)
            steps = steps + it
            stats = L.accumulate(stats, t)
            return x
        end
        for j = 1, #r do
            r[j] = -r[j]
        end
        local u, Y = solve(r), {}
        for j, p in ipairs(ports) do
            local b = {}
            for row = 1, #m.free do
                b[row] = p.load[row] or 0
            end
            Y[j] = solve(b)
        end
        local sr, si, br, bi = {}, {}, {}, {}
        for j, p in ipairs(ports) do
            sr[j], si[j], br[j], bi[j] = {}, {}, -g[j] - free_link(p, u), 0
            for l = 1, #ports do
                sr[j][l], si[j][l] = free_link(p, Y[l]), 0
            end
        end
        local di = require("luafemm-ac-circuits").dense(sr, si, br, bi)
        for j = 1, #u do
            for l = 1, #ports do
                u[j] = u[j] + Y[l][j] * di[l]
            end
        end
        local accepted = false
        for back = 0, 24 do
            local alpha = 2 ^ -back
            local a, current = {}, {}
            for n, node in ipairs(m.nodes) do
                a[n] = A[n] + (node.dof and alpha * (node.dof_sign or 1) * u[node.dof] or 0)
            end
            for j = 1, #I do
                current[j] = I[j] + alpha * di[j]
            end
            local _, _, _, rf, rc = residual(a, current, false)
            if math.max(rf, rc) < math.max(fr, cr) * (1 - 1e-4 * alpha) then
                A, I, accepted = a, current, true
                break
            end
        end
        check(accepted, "constraint Newton line search failed")
    end
end
return D
