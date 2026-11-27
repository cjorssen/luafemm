-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Private real-block GMRES with bounded complex incomplete factor recovery.
-- Reorthogonalisation and original-operator residual checks apply independently
-- of factor shifts. Complex symmetry or positive definiteness is not assumed.
-- @module luafemm-complex
local L = require("luafemm-linear")
local C = {}
local check, dot, norm = L.check, L.dot, L.norm
--- Solve (Kr+i Ki)x=br+i bi from borrowed one-based sparse rows.
-- Returns real/imaginary arrays, iterations, original relative residual, and
-- private telemetry. opts is the private algebra experiment contract in L.prepare.
function C.solve(Kr, Ki, br, bi, tolerance, opts)
    local p = L.prepare(Kr, Ki, br, bi, tolerance, opts, false)
    p.stats.solver = "gmres"
    local n, x = p.n, L.zeros(2 * p.n)
    if p.zero then
        return L.finish(p, x, 0, os.clock())
    end
    local pre = require("luafemm-factor").build(p)
    local function multiply(v)
        return L.multiply(p, v)
    end
    local total, started = 0, os.clock()
    while true do
        local residual, relative, backward = L.residual(p, x)
        if L.accept(p, relative, backward) or total == p.max_iterations then
            local ar, ai, it, res, stats = L.finish(p, x, total, started)
            if ar then
                return ar, ai, it, res, stats
            end
            x, residual = ai, it
        end
        check(
            total < p.max_iterations,
            "GMRES exhausted iteration budget before original-operator residual tolerance"
        )
        if total > 0 then
            p.stats.restarts = p.stats.restarts + 1
        end
        local r = pre(residual)
        local beta = norm(r)
        check(beta > 0 and beta < math.huge, "nonfinite preconditioned residual")
        local v, h, c, s, g = { r }, {}, {}, {}, { beta }
        for i = 1, 2 * n do
            r[i] = r[i] / beta
        end
        local steps = math.min(p.restart, 2 * n, p.max_iterations - total)
        for k = 1, steps do
            local w = pre(multiply(v[k]))
            local wnorm = norm(w)
            h[k] = {}
            -- Two modified Gram--Schmidt passes prevent loss of orthogonality.
            for pass = 1, 2 do
                for j = 1, k do
                    local d = dot(w, v[j])
                    h[k][j] = (h[k][j] or 0) + d
                    for i = 1, 2 * n do
                        w[i] = w[i] - d * v[j][i]
                    end
                end
            end
            h[k][k + 1] = norm(w)
            local breakdown = h[k][k + 1] <= 1e-14 * math.max(wnorm, 1e-300)
            if not breakdown then
                v[k + 1] = w
                for i = 1, 2 * n do
                    w[i] = w[i] / h[k][k + 1]
                end
            end
            for j = 1, k - 1 do
                h[k][j], h[k][j + 1] =
                    c[j] * h[k][j] + s[j] * h[k][j + 1], -s[j] * h[k][j] + c[j] * h[k][j + 1]
            end
            local d = norm({ h[k][k], h[k][k + 1] })
            check(d > 0 and d < math.huge, "GMRES breakdown")
            c[k], s[k] = h[k][k] / d, h[k][k + 1] / d
            h[k][k], h[k][k + 1] = d, 0
            g[k + 1], g[k] = -s[k] * g[k], c[k] * g[k]
            total = total + 1
            if breakdown or math.abs(g[k + 1]) <= beta * tolerance or k == steps then
                local y = {}
                for j = k, 1, -1 do
                    local value = g[j]
                    for l = j + 1, k do
                        value = value - h[l][j] * y[l]
                    end
                    y[j] = value / h[j][j]
                end
                for j = 1, k do
                    for i = 1, 2 * n do
                        x[i] = x[i] + v[j][i] * y[j]
                    end
                end
                break
            end
        end
    end
end
return C
