-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Private sparse algebra shared by magnetic and auxiliary solves.
-- Input rows/loads are borrowed read-only. Equilibration, ordering and factors
-- are owned by one solve; no numerical factor survives a changed operator.
-- Permutations concern algebraic unknowns only, never mesh or file identities.
-- @module luafemm-linear
local L = {}
local abs, max, sqrt = math.abs, math.max, math.sqrt
function L.check(ok, message)
    assert(ok, "luafemm linear: " .. message)
end
local check = L.check
function L.finite(x)
    return type(x) == "number" and x == x and abs(x) < math.huge
end
local finite = L.finite
function L.zeros(n)
    local a = {}
    for i = 1, n do
        a[i] = 0
    end
    return a
end
-- Scaled sum of squares avoids overflow/underflow of individual squares.
-- A mathematical norm beyond the finite range is an error, never convergence.
function L.norm(a)
    local scale, sum = 0, 1
    for i = 1, #a do
        check(finite(a[i]), "nonfinite vector")
        local v = abs(a[i])
        if v > 0 then
            if v > scale then
                sum, scale = 1 + sum * (scale / v) ^ 2, v
            else
                sum = sum + (v / scale) ^ 2
            end
        end
    end
    local value = scale * sqrt(sum)
    check(finite(value), "vector norm exceeds finite range")
    return value
end
function L.dot(a, b)
    local value = 0
    for i = 1, #a do
        value = value + a[i] * b[i]
    end
    check(finite(value), "nonfinite inner product")
    return value
end
function L.divide(ar, ai, br, bi)
    local scale = max(abs(br), abs(bi))
    check(scale > 0 and finite(scale), "singular or nonfinite factor pivot")
    br, bi = br / scale, bi / scale
    local d = br * br + bi * bi
    local r, s = (ar / scale * br + ai / scale * bi) / d, (ai / scale * br - ar / scale * bi) / d
    check(finite(r) and finite(s), "factor division exceeds finite range")
    return r, s
end
-- Form a*b*c/d without overflowing an intermediate product. Used only for
-- mapping solutions between units, outside iterative multiplication loops.
local function product(a, b, c, d)
    if a == 0 or b == 0 or c == 0 then
        return 0
    end
    local fa, ea = math.frexp(a)
    local fb, eb = math.frexp(b)
    local fc, ec = math.frexp(c)
    local fd, ed = math.frexp(d)
    local value = math.ldexp(fa * fb * fc / fd, ea + eb + ec - ed)
    check(finite(value) and value ~= 0, "solution scaling exceeds finite range")
    return value
end
local function array(a, n, label)
    check(type(a) == "table", label .. " must be an array")
    local count = 0
    for k in pairs(a) do
        check(
            type(k) == "number" and k % 1 == 0 and k >= 1 and k <= n,
            "invalid " .. label .. " index"
        )
        count = count + 1
    end
    check(count == n, "incomplete " .. label)
end
local function permutation(cols, n, choice)
    local order = {}
    if choice == "mesh" then
        for i = 1, n do
            order[i] = i
        end
        return order
    end
    -- Deterministic reverse Cuthill--McKee on the undirected sparsity graph.
    local adj, degree, roots, seen = {}, {}, {}, {}
    for i = 1, n do
        adj[i], roots[i] = {}, i
    end
    for i = 1, n do
        for _, j in ipairs(cols[i]) do
            if j ~= i then
                adj[i][j], adj[j][i] = true, true
            end
        end
    end
    for i = 1, n do
        degree[i] = 0
        for _ in pairs(adj[i]) do
            degree[i] = degree[i] + 1
        end
    end
    local function less(a, b)
        return degree[a] == degree[b] and a < b or degree[a] < degree[b]
    end
    table.sort(roots, less)
    for _, root in ipairs(roots) do
        if not seen[root] then
            local queue, at = { root }, 1
            seen[root] = true
            while at <= #queue do
                local i, nexts = queue[at], {}
                at = at + 1
                order[#order + 1] = i
                for j in pairs(adj[i]) do
                    if not seen[j] then
                        nexts[#nexts + 1] = j
                    end
                end
                table.sort(nexts, less)
                for _, j in ipairs(nexts) do
                    seen[j], queue[#queue + 1] = true, j
                end
            end
        end
    end
    local reverse = {}
    for i = 1, n do
        reverse[i] = order[n + 1 - i]
    end
    return reverse
end

--- Validate and finalise a borrowed operator with resolved kernel options.
-- Symmetric row/column scaling preserves the SPD class. Loads are normalised
-- independently, so uniformly tiny/large systems have ordinary-sized iterates.
-- opts: ordering mesh/rcm, storage packed/rows, factor auto/diagonal, fill 0..32,
-- max_iterations (4000), restart (40). Storage and fill remain private experiments;
-- public model policy is resolved below before calling this kernel.
function L.prepare(kr, ki, br, bi, tolerance, opts, spd)
    local start = os.clock()
    opts = opts or {}
    check(type(opts) == "table", "invalid solver options")
    local allowed = {
        ordering = true,
        storage = true,
        factor = true,
        fill = true,
        max_iterations = true,
        restart = true,
    }
    for k in pairs(opts) do
        check(allowed[k], "unknown solver option " .. tostring(k))
    end
    check(finite(tolerance) and tolerance > 0 and tolerance < 1, "tolerance must be in (0,1)")
    check(type(br) == "table", "invalid real load")
    local n = #br
    array(br, n, "real load")
    array(bi, n, "imaginary load")
    array(kr, n, "real rows")
    array(ki, n, "imaginary rows")
    local p = {
        n = n,
        tolerance = tolerance,
        spd = spd,
        cols = {},
        r = {},
        im = {},
        rows = {},
        irows = {},
        b = {},
        scale = {},
        stats = {
            factor_seconds = 0,
            iteration_seconds = 0,
            restarts = 0,
            retries = 0,
            shift = 0,
            factor_nonzeros = 0,
        },
    }
    p.ordering, p.storage, p.factor =
        opts.ordering or "mesh", opts.storage or "packed", opts.factor or "auto"
    p.fill, p.max_iterations, p.restart =
        opts.fill or 0, opts.max_iterations or 4000, opts.restart or 40
    check(p.ordering == "mesh" or p.ordering == "rcm", "invalid linear ordering")
    check(p.storage == "packed" or p.storage == "rows", "invalid sparse storage")
    check(p.factor == "auto" or p.factor == "diagonal", "invalid factor policy")
    for _, k in ipairs({ "fill", "max_iterations", "restart" }) do
        check(finite(p[k]) and p[k] % 1 == 0 and p[k] >= (k == "fill" and 0 or 1), "invalid " .. k)
    end
    check(
        p.fill <= 32 and p.restart <= 4000 and p.max_iterations <= 1000000,
        "solver budget exceeds limit"
    )
    local cols, scales, loadscale = {}, {}, 0
    for i = 1, n do
        check(finite(br[i]) and finite(bi[i]), "nonfinite right-hand side")
        check(not spd or bi[i] == 0, "CG requires a real operator and load")
        loadscale = max(loadscale, abs(br[i]), abs(bi[i]))
        local set, rowmax = { [i] = true }, 0
        for _, row in ipairs({ kr[i], ki[i] }) do
            check(type(row) == "table", "invalid sparse row")
            for j, value in pairs(row) do
                check(
                    type(j) == "number" and j % 1 == 0 and j >= 1 and j <= n,
                    "invalid sparse column"
                )
                check(finite(value), "nonfinite matrix coefficient")
                if value ~= 0 then
                    set[j] = true
                end
                rowmax = max(rowmax, abs(value))
            end
        end
        check(rowmax > 0, "singular operator: empty row " .. i)
        if spd then
            check((kr[i][i] or 0) > 0, "CG requires positive diagonal and SPD operator")
            for _, v in pairs(ki[i]) do
                check(v == 0, "CG requires a real operator")
            end
        end
        cols[i], scales[i] = {}, sqrt(rowmax)
        for j in pairs(set) do
            cols[i][#cols[i] + 1] = j
        end
        table.sort(cols[i])
    end
    p.order = permutation(cols, n, p.ordering)
    local inverse = {}
    for i, old in ipairs(p.order) do
        inverse[old] = i
    end
    p.original = { kr, ki, br, bi }
    p.loadscale = loadscale > 0 and loadscale or 1
    local rhsmax, scalemax = 0, 0
    for i, old in ipairs(p.order) do
        p.scale[i], p.cols[i], p.rows[i], p.irows[i] = scales[old], {}, {}, {}
        scalemax = max(scalemax, p.scale[i])
        p.b[i], p.b[n + i] =
            (br[old] / p.loadscale) / scales[old], (bi[old] / p.loadscale) / scales[old]
        check(br[old] == 0 or p.b[i] ~= 0, "load component lost during equilibration")
        check(bi[old] == 0 or p.b[n + i] ~= 0, "load component lost during equilibration")
        rhsmax = max(rhsmax, abs(p.b[i]), abs(p.b[n + i]))
        for _, col in ipairs(cols[old]) do
            local j = inverse[col]
            local r, s =
                (kr[old][col] or 0) / scales[old] / scales[col],
                (ki[old][col] or 0) / scales[old] / scales[col]
            check(finite(r) and finite(s), "equilibrated coefficient exceeds finite range")
            check((kr[old][col] or 0) == 0 or r ~= 0, "coefficient lost during equilibration")
            check((ki[old][col] or 0) == 0 or s ~= 0, "coefficient lost during equilibration")
            p.cols[i][#p.cols[i] + 1], p.rows[i][j], p.irows[i][j] = j, r, s
        end
        table.sort(p.cols[i])
        p.r[i], p.im[i] = {}, {}
        for k, j in ipairs(p.cols[i]) do
            p.r[i][k], p.im[i][k] = p.rows[i][j], p.irows[i][j]
        end
    end
    p.rhsmax, p.scalemax = rhsmax > 0 and rhsmax or 1, scalemax > 0 and scalemax or 1
    local weighted = {}
    for i = 1, n do
        local br, bi = p.b[i], p.b[n + i]
        p.b[i], p.b[n + i] = br / p.rhsmax, bi / p.rhsmax
        check(br == 0 or p.b[i] ~= 0, "load component lost during normalisation")
        check(bi == 0 or p.b[n + i] ~= 0, "load component lost during normalisation")
        weighted[i], weighted[n + i] =
            p.b[i] * (p.scale[i] / p.scalemax), p.b[n + i] * (p.scale[i] / p.scalemax)
        if spd then
            for _, j in ipairs(p.cols[i]) do
                local a, b = p.rows[i][j], p.rows[j][i] or 0
                check(
                    abs(a - b) <= 1e-12 * max(abs(a), abs(b), abs(p.rows[i][i]), abs(p.rows[j][j])),
                    "CG requires a symmetric operator"
                )
            end
        end
    end
    if spd then
        for i = 1, n do
            p.b[n + i] = nil
        end
    end
    p.bnorm, p.original_bnorm = L.norm(p.b), L.norm(weighted)
    p.zero = loadscale == 0
    check(p.zero or p.original_bnorm > 0, "load lost during equilibration")
    p.stats.finalize_seconds = os.clock() - start
    return p
end

function L.multiply(p, x)
    local n, y = p.n, {}
    local packed = p.storage == "packed"
    for i = 1, n do
        local r, s = 0, 0
        local cols, real, imag = p.cols[i], p.r[i], p.im[i]
        if p.spd then
            if packed then
                for k = 1, #cols do
                    r = r + real[k] * x[cols[k]]
                end
            else
                for k = 1, #cols do
                    local j = cols[k]
                    r = r + p.rows[i][j] * x[j]
                end
            end
        elseif packed then
            for k = 1, #cols do
                local j, a, b = cols[k], real[k], imag[k]
                r, s = r + a * x[j] - b * x[n + j], s + b * x[j] + a * x[n + j]
            end
        else
            for k = 1, #cols do
                local j = cols[k]
                local a, b = p.rows[i][j], p.irows[i][j]
                r, s = r + a * x[j] - b * x[n + j], s + b * x[j] + a * x[n + j]
            end
        end
        check(finite(r) and finite(s), "nonfinite matrix product")
        y[i] = r
        if not p.spd then
            y[n + i] = s
        end
    end
    return y
end

--- Recompute against the unshifted operator, undoing row scaling for its norm.
-- Backward error groups the real/imaginary parts of each complex row, using
-- |r_real|+|r_imag| over sum (|A_real|+|A_imag|)*max_j(|x_real|+|x_imag|)+|b|_1.
-- A row-scaled normwise denominator remains meaningful in near-zero regions;
-- componentwise relative requirements there can be unattainable after rounding.
-- Neither metric is a forward error or field discretisation error estimate.
function L.residual(p, x)
    local ax, weighted, backward = L.multiply(p, x), {}, 0
    local xmax = 0
    for i = 1, p.n do
        xmax = max(xmax, abs(x[i]) + abs(x[p.n + i] or 0))
    end
    for i = 1, p.n do
        local denom = abs(p.b[i]) + abs(p.b[p.n + i] or 0)
        for k, j in ipairs(p.cols[i]) do
            denom = denom + (abs(p.r[i][k]) + abs(p.im[i][k])) * xmax
        end
        ax[i] = p.b[i] - ax[i]
        local ri = p.spd and 0 or p.b[p.n + i] - ax[p.n + i]
        backward = max(backward, denom > 0 and (abs(ax[i]) + abs(ri)) / denom or 0)
        weighted[i] = ax[i] * (p.scale[i] / p.scalemax)
        if not p.spd then
            ax[p.n + i], weighted[p.n + i] = ri, ri * (p.scale[i] / p.scalemax)
        end
    end
    local relative = p.zero and L.norm(weighted) or L.norm(weighted) / p.original_bnorm
    check(finite(relative) and finite(backward), "nonfinite residual metric")
    return ax, relative, backward
end
function L.accept(p, relative, backward)
    return relative <= p.tolerance and backward <= max(p.tolerance, 1e-14)
end
--- Residual of the returned physical unknowns using the caller's original
-- coefficients, independently of the equilibrated multiplication/factorisation.
-- Divide products by the load scale using binary exponents, avoiding intermediate
-- overflow for uniformly extreme units. The zero-load case is handled exactly.
function L.original_residual(p, xr, xi)
    local kr, ki, br, bi = table.unpack(p.original)
    local residual, load, backward = {}, {}, 0
    local xmax = 0
    for i = 1, p.n do
        xmax = max(xmax, abs(xr[i]) + abs(xi[i]))
    end
    check(finite(xmax), "solution norm exceeds finite range")
    for i = 1, p.n do
        local r, s = br[i] / p.loadscale, bi[i] / p.loadscale
        local denom = abs(r) + abs(s)
        local cr, cs = 0, 0
        local function add(a, c, b)
            local t = a + b
            if abs(a) >= abs(b) then
                c = c + ((a - t) + b)
            else
                c = c + ((b - t) + a)
            end
            return t, c
        end
        load[i], load[p.n + i] = r, s
        local cols = {}
        for j in pairs(kr[i]) do
            cols[j] = true
        end
        for j in pairs(ki[i]) do
            cols[j] = true
        end
        local ids = {}
        for j in pairs(cols) do
            ids[#ids + 1] = j
        end
        table.sort(ids)
        for _, j in ipairs(ids) do
            local ar, ai = kr[i][j] or 0, ki[i][j] or 0
            local rr = product(ar, xr[j], 1, p.loadscale)
            local ii = product(ai, xi[j], 1, p.loadscale)
            local ir = product(ai, xr[j], 1, p.loadscale)
            local ri = product(ar, xi[j], 1, p.loadscale)
            r, cr = add(r, cr, -rr)
            r, cr = add(r, cr, ii)
            s, cs = add(s, cs, -ir)
            s, cs = add(s, cs, -ri)
            denom = denom
                + abs(product(ar, xmax, 1, p.loadscale))
                + abs(product(ai, xmax, 1, p.loadscale))
        end
        check(finite(denom), "backward-error denominator exceeds finite range")
        r, s = r + cr, s + cs
        residual[i], residual[p.n + i] = r, s
        backward = max(backward, denom > 0 and (abs(r) + abs(s)) / denom or 0)
    end
    return L.norm(residual) / (p.zero and 1 or L.norm(load)), backward, residual
end
function L.finish(p, x, iterations, started)
    local xr, xi, rounded = {}, {}, {}
    for i, old in ipairs(p.order) do
        xr[old] = product(x[i], p.loadscale, p.rhsmax, p.scale[i])
        xi[old] = product(x[p.n + i] or 0, p.loadscale, p.rhsmax, p.scale[i])
        rounded[i] = product(xr[old], p.scale[i], 1 / p.rhsmax, p.loadscale)
        if not p.spd then
            rounded[p.n + i] = product(xi[old], p.scale[i], 1 / p.rhsmax, p.loadscale)
        end
    end
    local relative, backward, actual = L.original_residual(p, xr, xi)
    if not L.accept(p, relative, backward) then
        local correction = {}
        for i, old in ipairs(p.order) do
            correction[i] = product(actual[old], 1, 1 / p.rhsmax, p.scale[i])
            if not p.spd then
                correction[p.n + i] = product(actual[p.n + old], 1, 1 / p.rhsmax, p.scale[i])
            end
        end
        return nil, rounded, correction, relative, backward
    end
    p.stats.iteration_seconds = os.clock() - started
    p.stats.relative_residual, p.stats.backward_error, p.stats.iterations =
        relative, backward, iterations
    return xr, xi, iterations, relative, p.stats
end

--- Accumulate private numerical telemetry without conflating Newton residuals.
function L.accumulate(total, item)
    total = total
        or {
            solves = 0,
            cg_solves = 0,
            gmres_solves = 0,
            ic_solves = 0,
            ilu_solves = 0,
            diagonal_solves = 0,
            finalize_seconds = 0,
            factor_seconds = 0,
            iteration_seconds = 0,
            restarts = 0,
            retries = 0,
            iterations = 0,
            factor_nonzeros = 0,
            relative_residual = 0,
            backward_error = 0,
            shift = 0,
        }
    total.solves = total.solves + 1
    if item.solver then
        total[item.solver .. "_solves"] = total[item.solver .. "_solves"] + 1
    end
    if item.preconditioner then
        local key = item.preconditioner .. "_solves"
        total[key] = total[key] + 1
    end
    for _, k in ipairs({
        "finalize_seconds",
        "factor_seconds",
        "iteration_seconds",
        "restarts",
        "retries",
        "iterations",
    }) do
        total[k] = total[k] + item[k]
    end
    for _, k in ipairs({ "factor_nonzeros", "relative_residual", "backward_error", "shift" }) do
        total[k] = max(total[k], item[k])
    end
    return total
end
--- Public numerical policy, owned by the model and included in cache identity.
L.option_defaults = {
    linear_solver = "auto",
    preconditioner = "auto",
    linear_tolerance = "auto",
    linear_max_iterations = 4000,
    gmres_restart = 40,
    linear_ordering = "mesh",
}
--- Apply validated numerical-only overrides to an owned derived model.
function L.override(m, values)
    check(values == nil or type(values) == "table", "invalid numerical overrides")
    for k, v in pairs(values or {}) do
        check(L.option_defaults[k] ~= nil, "unknown numerical override " .. tostring(k))
        m.options[k] = v
    end
    L.configure(m.options)
end
function L.configure(o)
    for k, v in pairs(L.option_defaults) do
        if o[k] == nil then
            o[k] = v
        end
    end
    check(
        o.linear_solver == "auto" or o.linear_solver == "cg" or o.linear_solver == "gmres",
        "unknown linear solver"
    )
    check(
        o.preconditioner == "auto"
            or o.preconditioner == "ic0"
            or o.preconditioner == "ilu0"
            or o.preconditioner == "diagonal",
        "unknown preconditioner"
    )
    check(o.linear_ordering == "mesh" or o.linear_ordering == "rcm", "unknown linear ordering")
    check(
        o.linear_tolerance == "auto"
            or finite(o.linear_tolerance) and o.linear_tolerance > 0 and o.linear_tolerance < 1,
        "linear tolerance must be auto or in (0,1)"
    )
    for k, limit in pairs({ linear_max_iterations = 1000000, gmres_restart = 4000 }) do
        check(
            finite(o[k]) and o[k] % 1 == 0 and o[k] >= 1 and o[k] <= limit,
            k .. " must be an integer in 1.." .. limit
        )
    end
end
--- Resolve by operator class, not frequency alone (fundamental laws are real).
-- Weight potentials always use SPD CG; only ordering, budget and diagonal
-- selection are inherited. Their measurement tolerance remains independent.
function L.policy(o, kind, tolerance)
    L.configure(o)
    local solver = kind == "general" and "gmres" or "cg"
    if kind ~= "weight" and o.linear_solver ~= "auto" then
        solver = o.linear_solver
    end
    check(kind ~= "general" or solver == "gmres", "CG cannot solve this general/complex operator")
    local factor = o.preconditioner
    if kind == "weight" and factor ~= "diagonal" then
        factor = "auto"
    end
    check(factor ~= "ic0" or solver == "cg", "ic0 requires CG")
    check(factor ~= "ilu0" or solver == "gmres", "ilu0 requires GMRES")
    if kind ~= "weight" and o.linear_tolerance ~= "auto" then
        tolerance = o.linear_tolerance
    end
    return solver,
        tolerance,
        {
            ordering = o.linear_ordering,
            max_iterations = o.linear_max_iterations,
            restart = o.gmres_restart,
            factor = factor == "diagonal" and "diagonal" or "auto",
        }
end
function L.real(m, K, b, tolerance, weight)
    local solver, tol, opts = L.policy(m.options, weight and "weight" or "spd", tolerance)
    if solver == "cg" then
        return require("luafemm").solve_spd(K, b, tol, opts)
    end
    local ki, bi = {}, {}
    for i = 1, #b do
        ki[i], bi[i] = {}, 0
    end
    local x, _, it, _, stats = require("luafemm-complex").solve(K, ki, b, bi, tol, opts)
    return x, it, stats
end
function L.complex(m, kr, ki, br, bi, tolerance)
    local _, tol, opts = L.policy(m.options, "general", tolerance or m.options.tolerance)
    return require("luafemm-complex").solve(kr, ki, br, bi, tol, opts)
end
--- Read numerical diagnostics without flattening or mutating cached records.
-- Imported answers have no local solve history. A zero-Newton native solve
-- has zero linear work, rather than unavailable diagnostics.
function L.statistic(m, key)
    local names = {
        solves = "solves",
        iterations = "iterations",
        restarts = "restarts",
        retries = "retries",
        nonzeros = "factor_nonzeros",
        shift = "shift",
        residual = "relative_residual",
        backward_error = "backward_error",
        finalize_seconds = "finalize_seconds",
        factor_seconds = "factor_seconds",
        iteration_seconds = "iteration_seconds",
        cg_solves = "cg_solves",
        gmres_solves = "gmres_solves",
        diagonal_solves = "diagonal_solves",
    }
    local field = names[key:match("^linear_(.+)$")]
    if not field then
        return nil
    end
    check(
        m.stats and not m.stats.imported,
        "linear statistics unavailable for an unsolved or imported answer"
    )
    return m.stats.linear and m.stats.linear[field] or 0
end

return L
