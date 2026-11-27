-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Private bounded incomplete factors; the physical operator is never changed.
-- IC(0) for SPD operators, complex ILU(0) for general operators. Experimental
-- fill adds at most p.fill two-hop candidates per row, ordered by connection
-- strength then column ID. This is bounded pattern fill, not a full ILUT port.
-- @module luafemm-factor
local S = require("luafemm-linear")
local F = {}
local abs, max, sqrt = math.abs, math.max, math.sqrt
local function pattern(p)
    local cols = {}
    for i = 1, p.n do
        local set, candidates = {}, {}
        for _, j in ipairs(p.cols[i]) do
            set[j] = true
        end
        if p.fill > 0 then
            for _, j in ipairs(p.cols[i]) do
                local w = abs(p.rows[i][j]) + abs(p.irows[i][j])
                for _, k in ipairs(p.cols[j]) do
                    if not set[k] and (not p.spd or k < i) then
                        candidates[k] = (candidates[k] or 0)
                            + w * (abs(p.rows[j][k]) + abs(p.irows[j][k]))
                    end
                end
            end
            local ids = {}
            for k in pairs(candidates) do
                ids[#ids + 1] = k
            end
            table.sort(ids, function(a, b)
                return candidates[a] == candidates[b] and a < b or candidates[a] > candidates[b]
            end)
            for k = 1, math.min(p.fill, #ids) do
                set[ids[k]] = true
            end
        end
        cols[i] = {}
        for j in pairs(set) do
            if not p.spd or j < i then
                cols[i][#cols[i] + 1] = j
            end
        end
        table.sort(cols[i])
    end
    return cols
end
local function incomplete(p, cols, shift)
    local n, lr, li, diag = p.n, {}, {}, {}
    local nonzeros = 0
    for i = 1, n do
        lr[i], li[i] = {}, {}
        if p.spd then
            for _, j in ipairs(cols[i]) do
                local value = p.rows[i][j] or 0
                for _, k in ipairs(cols[i]) do
                    if k >= j then
                        break
                    end
                    value = value - lr[i][k] * (lr[j][k] or 0)
                end
                lr[i][j] = value / diag[j]
                if not S.finite(lr[i][j]) then
                    return nil
                end
            end
            local value = p.rows[i][i] * (1 + shift)
            for _, j in ipairs(cols[i]) do
                value = value - lr[i][j] ^ 2
            end
            if not S.finite(value) or value <= 1e-14 * p.rows[i][i] then
                return nil
            end
            diag[i] = sqrt(value)
        else
            for _, j in ipairs(cols[i]) do
                lr[i][j], li[i][j] = p.rows[i][j] or 0, p.irows[i][j] or 0
            end
            local dr, di = lr[i][i], li[i][i]
            -- Preserve a nonzero diagonal's complex phase. Zero diagonals
            -- take a real row scale for the preconditioner alone.
            if max(abs(dr), abs(di)) > 0 then
                lr[i][i], li[i][i] = dr * (1 + shift), di * (1 + shift)
            else
                local row = 0
                for k = 1, #p.r[i] do
                    row = max(row, abs(p.r[i][k]), abs(p.im[i][k]))
                end
                lr[i][i] = shift * row
            end
            for _, j in ipairs(cols[i]) do
                if j >= i then
                    break
                end
                local ar, ai = S.divide(lr[i][j], li[i][j], lr[j][j], li[j][j])
                lr[i][j], li[i][j] = ar, ai
                for _, k in ipairs(cols[j]) do
                    if k > j and lr[i][k] ~= nil then
                        lr[i][k] = lr[i][k] - ar * lr[j][k] + ai * li[j][k]
                        li[i][k] = li[i][k] - ai * lr[j][k] - ar * li[j][k]
                        if not S.finite(lr[i][k]) or not S.finite(li[i][k]) then
                            return nil
                        end
                    end
                end
            end
            if max(abs(lr[i][i]), abs(li[i][i])) <= 1e-14 then
                return nil
            end
        end
        nonzeros = nonzeros + #cols[i] + (p.spd and 1 or 0)
    end
    local function apply(x)
        local y = {}
        for i = 1, #x do
            y[i] = x[i]
        end
        if p.spd then
            for i = 1, n do
                for _, j in ipairs(cols[i]) do
                    y[i] = y[i] - lr[i][j] * y[j]
                end
                y[i] = y[i] / diag[i]
            end
            for i = n, 1, -1 do
                y[i] = y[i] / diag[i]
                for _, j in ipairs(cols[i]) do
                    y[j] = y[j] - lr[i][j] * y[i]
                end
            end
        else
            for i = 1, n do
                for _, j in ipairs(cols[i]) do
                    if j >= i then
                        break
                    end
                    y[i], y[n + i] =
                        y[i] - lr[i][j] * y[j] + li[i][j] * y[n + j],
                        y[n + i] - li[i][j] * y[j] - lr[i][j] * y[n + j]
                end
            end
            for i = n, 1, -1 do
                for _, j in ipairs(cols[i]) do
                    if j > i then
                        y[i], y[n + i] =
                            y[i] - lr[i][j] * y[j] + li[i][j] * y[n + j],
                            y[n + i] - li[i][j] * y[j] - lr[i][j] * y[n + j]
                    end
                end
                y[i], y[n + i] = S.divide(y[i], y[n + i], lr[i][i], li[i][i])
            end
        end
        for i = 1, #x do
            S.check(S.finite(y[i]), "nonfinite preconditioned vector")
        end
        return y
    end
    return apply, nonzeros
end
--- Build once per operator. Failed attempts are discarded as a whole.
-- The finite ladder ends at a row-dominance shift, then a diagonal fallback.
-- A usable factor is not a promise that the Krylov budget will suffice.
function F.build(p)
    local start = os.clock()
    local cols = pattern(p)
    local dominance = 0
    for i = 1, p.n do
        local off = 0
        for k, j in ipairs(p.cols[i]) do
            if j ~= i then
                off = off + abs(p.r[i][k]) + abs(p.im[i][k])
            end
        end
        local d = max(abs(p.rows[i][i]), abs(p.irows[i][i]))
        dominance = max(dominance, d > 0 and 2 * off / d or 2)
    end
    local shifts = { 0, 0.001, 0.01, 0.1, 1 }
    if S.finite(dominance) and dominance > 1 then
        shifts[#shifts + 1] = dominance
    end
    if p.factor ~= "diagonal" then
        for attempt, shift in ipairs(shifts) do
            local ok, apply, count = pcall(incomplete, p, cols, shift)
            if ok and apply then
                p.stats.factor_seconds = os.clock() - start
                p.stats.retries, p.stats.shift, p.stats.factor_nonzeros = attempt - 1, shift, count
                p.stats.preconditioner = p.spd and "ic" or "ilu"
                return apply
            end
        end
    end
    p.stats.retries = p.factor == "diagonal" and 0 or #shifts
    p.stats.preconditioner, p.stats.factor_nonzeros = "diagonal", p.n
    p.stats.factor_seconds = os.clock() - start
    return function(x)
        local y = {}
        for i = 1, p.n do
            local dr, di = p.rows[i][i], p.irows[i][i]
            if max(abs(dr), abs(di)) <= 1e-14 then
                dr, di = 1, 0
            end
            if p.spd then
                y[i] = x[i] / dr
            else
                y[i], y[p.n + i] = S.divide(x[i], x[p.n + i], dr, di)
            end
        end
        return y
    end
end
return F
