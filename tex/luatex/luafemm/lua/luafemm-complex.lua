-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Sparse complex linear algebra for harmonic finite elements.
-- Real block GMRES with complex ILU(0), reorthogonalisation and true-residual
-- acceptance. No positive-definiteness or Hermitian assumption is made.
-- @module luafemm-complex
local C = {}
local function check(v, s)
    assert(v, "luafemm complex: " .. s)
end
local function dot(a, b)
    local s = 0
    for i = 1, #a do
        s = s + a[i] * b[i]
    end
    return s
end
local function norm(a)
    return math.sqrt(dot(a, a))
end
local function zeros(n)
    local a = {}
    for i = 1, n do
        a[i] = 0
    end
    return a
end
local function divide(ar, ai, br, bi)
    local scale = math.max(math.abs(br), math.abs(bi))
    check(scale > 0 and scale < math.huge, "singular or nonfinite ILU pivot")
    br, bi = br / scale, bi / scale
    local d = br * br + bi * bi
    return (ar / scale * br + ai / scale * bi) / d, (ai / scale * br - ar / scale * bi) / d
end
--- Solve (Kr+i Ki)(xr+i xi)=br+i bi using sparse one-based row tables.
-- Missing entries are zero. Returns real/imaginary arrays, iteration count
-- and relative true residual. Failure raises; unconverged fields are not saved.
-- @tparam table Kr Sparse real matrix rows.
-- @tparam table Ki Sparse imaginary matrix rows.
-- @tparam table br Real right-hand side.
-- @tparam table bi Imaginary right-hand side.
-- @tparam number tolerance Positive relative true-residual target.
-- @return Real solution, imaginary solution, iterations, relative residual.
function C.solve(Kr, Ki, br, bi, tolerance)
    local n = #br
    if n == 0 then
        return {}, {}, 0, 0
    end
    local cols, lr, li = {}, {}, {}
    for i = 1, n do
        local set = {}
        lr[i], li[i], cols[i] = {}, {}, {}
        for j in pairs(Kr[i]) do
            set[j] = true
        end
        for j in pairs(Ki[i]) do
            set[j] = true
        end
        for j in pairs(set) do
            cols[i][#cols[i] + 1] = j
            lr[i][j], li[i][j] = Kr[i][j] or 0, Ki[i][j] or 0
        end
        table.sort(cols[i])
    end
    for i = 1, n do
        for _, j in ipairs(cols[i]) do
            if j >= i then
                break
            end
            local ar, ai = divide(lr[i][j], li[i][j], lr[j][j] or 0, li[j][j] or 0)
            lr[i][j], li[i][j] = ar, ai
            for _, k in ipairs(cols[j]) do
                if k > j and lr[i][k] then
                    lr[i][k] = lr[i][k] - ar * lr[j][k] + ai * li[j][k]
                    li[i][k] = li[i][k] - ar * li[j][k] - ai * lr[j][k]
                end
            end
        end
        check(lr[i][i] and lr[i][i] ^ 2 + li[i][i] ^ 2 > 0, "zero ILU diagonal")
    end
    local function multiply(x)
        local y = zeros(2 * n)
        for i = 1, n do
            for _, j in ipairs(cols[i]) do
                local r, s = Kr[i][j] or 0, Ki[i][j] or 0
                y[i] = y[i] + r * x[j] - s * x[n + j]
                y[n + i] = y[n + i] + s * x[j] + r * x[n + j]
            end
        end
        return y
    end
    local function pre(x)
        local y = {}
        for i = 1, #x do
            y[i] = x[i]
        end
        for i = 1, n do
            for _, j in ipairs(cols[i]) do
                if j >= i then
                    break
                end
                y[i] = y[i] - lr[i][j] * y[j] + li[i][j] * y[n + j]
                y[n + i] = y[n + i] - li[i][j] * y[j] - lr[i][j] * y[n + j]
            end
        end
        for i = n, 1, -1 do
            for _, j in ipairs(cols[i]) do
                if j > i then
                    y[i] = y[i] - lr[i][j] * y[j] + li[i][j] * y[n + j]
                    y[n + i] = y[n + i] - li[i][j] * y[j] - lr[i][j] * y[n + j]
                end
            end
            y[i], y[n + i] = divide(y[i], y[n + i], lr[i][i], li[i][i])
        end
        return y
    end
    local b, x = {}, zeros(2 * n)
    for i = 1, n do
        b[i], b[n + i] = br[i], bi[i]
    end
    local scale = norm(b)
    if scale == 0 then
        return zeros(n), zeros(n), 0, 0
    end
    local total = 0
    while total < 4000 do
        local ax = multiply(x)
        local residual = {}
        for i = 1, 2 * n do
            residual[i] = b[i] - ax[i]
        end
        local relative = norm(residual) / scale
        if relative <= tolerance then
            local xr, xi = {}, {}
            for i = 1, n do
                xr[i], xi[i] = x[i], x[n + i]
            end
            return xr, xi, total, relative
        end
        local r = pre(residual)
        local beta = norm(r)
        check(beta > 0 and beta < math.huge, "nonfinite preconditioned residual")
        local v, h, c, s, g = { r }, {}, {}, {}, { beta }
        for i = 1, 2 * n do
            r[i] = r[i] / beta
        end
        local steps = math.min(40, 2 * n, 4000 - total)
        for k = 1, steps do
            local w = pre(multiply(v[k]))
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
            local breakdown = h[k][k + 1] < 1e-30
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
            local d = math.sqrt(h[k][k] ^ 2 + h[k][k + 1] ^ 2)
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
    error("luafemm complex: GMRES failed to reach the true residual tolerance")
end
return C
