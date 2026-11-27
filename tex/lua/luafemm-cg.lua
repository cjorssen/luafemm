-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Private SPD solve, also used for force-weight potentials.
-- @module luafemm-cg
local L = require("luafemm-linear")
local C = {}
--- Solve sparse real SPD rows without modifying them; return x, iterations, stats.
-- Relative original-system tolerance; no hidden absolute floor. The model layer
-- owns gauge removal. A known indefinite operator must use a different solver.
function C.solve(K, b, tol, opts)
    local empty, zero = {}, {}
    for i = 1, #b do
        empty[i], zero[i] = {}, 0
    end
    local p = L.prepare(K, empty, b, zero, tol, opts, true)
    p.stats.solver = "cg"
    local x = L.zeros(p.n)
    local started = os.clock()
    local function finish(it)
        local a, second, third, _, stats = L.finish(p, x, it, started)
        if a then
            return a, third, stats
        end
        x = second
        return nil, third
    end
    if p.zero then
        return finish(0)
    end
    local pre = require("luafemm-factor").build(p)
    started = os.clock()
    local r = {}
    for i = 1, #p.b do
        r[i] = p.b[i]
    end
    local z = pre(r)
    local direction = {}
    for i = 1, #z do
        direction[i] = z[i]
    end
    local rz = L.dot(r, z)
    local last_relative, last_backward
    for it = 1, p.max_iterations do
        local q = L.multiply(p, direction)
        local pq = L.dot(direction, q)
        L.check(
            pq > 0 and rz > 0,
            "non-positive CG curvature; operator is not SPD or iteration broke down"
        )
        local alpha = rz / pq
        for i = 1, #x do
            x[i], r[i] = x[i] + alpha * direction[i], r[i] - alpha * q[i]
        end
        -- Inspect the true residual periodically, but retain conjugacy unless
        -- recurrence drift is significant. Unconditional CG(50) restarts can
        -- stall high-permeability machine operators that ordinary CG solves.
        local replace = false
        if it % 50 == 0 or L.norm(r) <= p.tolerance * p.bnorm or it == p.max_iterations then
            local actual, relative, backward = L.residual(p, x)
            last_relative, last_backward = relative, backward
            if L.accept(p, relative, backward) or it == p.max_iterations then
                local result, steps, stats = finish(it)
                if result then
                    return result, steps, stats
                end
                r, replace = steps, true
            else
                local drift = {}
                for i = 1, #r do
                    drift[i] = r[i] - actual[i]
                end
                if L.norm(drift) > 0.1 * L.norm(actual) then
                    r, replace = actual, true
                end
            end
        end
        z = pre(r)
        local next_rz = L.dot(r, z)
        if replace then
            p.stats.restarts = p.stats.restarts + 1
            for i = 1, #z do
                direction[i] = z[i]
            end
        else
            local beta = next_rz / rz
            for i = 1, #z do
                direction[i] = z[i] + beta * direction[i]
            end
        end
        rz = next_rz
    end
    error(
        string.format(
            "luafemm linear: CG exhausted iteration budget (relative %.6g, backward %.6g, target %.6g)",
            last_relative or -1,
            last_backward or -1,
            p.tolerance
        )
    )
end
return C
