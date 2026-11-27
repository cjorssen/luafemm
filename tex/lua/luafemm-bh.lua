-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Natural cubic H(B), including FEMM's smoothing and endpoint extrapolation.
-- Independent tridiagonal implementation; no upstream source is embedded.
-- @module luafemm-bh
local B = {}
--- Prepare an owned curve. Smoothing never mutates the input B-H table.
-- @tparam table input Strictly increasing {B,H} pairs starting at {0,0}, in SI.
-- @treturn table Owned smoothed points, nodal slopes and smoothing pass count.
function B.prepare(input)
    local t = {}
    for i, p in ipairs(input) do
        t[i] = { p[1], p[2] }
    end
    local n = #t
    for pass = 1, 1000 do
        local lower, diag, upper, rhs = {}, {}, {}, {}
        for i = 1, n do
            local l = i > 1 and t[i][1] - t[i - 1][1] or nil
            local r = i < n and t[i + 1][1] - t[i][1] or nil
            lower[i] = l and 2 / l or 0
            upper[i] = r and 2 / r or 0
            diag[i] = (l and 4 / l or 0) + (r and 4 / r or 0)
            rhs[i] = (l and 6 * (t[i][2] - t[i - 1][2]) / l ^ 2 or 0)
                + (r and 6 * (t[i + 1][2] - t[i][2]) / r ^ 2 or 0)
        end
        for i = 2, n do
            local q = lower[i] / diag[i - 1]
            diag[i] = diag[i] - q * upper[i - 1]
            rhs[i] = rhs[i] - q * rhs[i - 1]
        end
        local slopes = {}
        slopes[n] = rhs[n] / diag[n]
        for i = n - 1, 1, -1 do
            slopes[i] = (rhs[i] - upper[i] * slopes[i + 1]) / diag[i]
        end
        local good = true
        for i = 1, n - 1 do
            local l = t[i + 1][1] - t[i][1]
            local delta = (t[i + 1][2] - t[i][2]) / l
            local c0 = slopes[i]
            local c1 = 6 * delta - 4 * slopes[i] - 2 * slopes[i + 1]
            local c2 = 3 * (slopes[i] + slopes[i + 1] - 2 * delta)
            local function root(z)
                if z >= 0 and z <= 1 then
                    good = false
                end
            end
            if c2 == 0 then
                if c1 ~= 0 then
                    root(-c0 / c1)
                end
            else
                local d = c1 * c1 - 4 * c0 * c2
                if d > 0 then
                    root((-c1 - math.sqrt(d)) / (2 * c2))
                    root((-c1 + math.sqrt(d)) / (2 * c2))
                end
            end
        end
        if good then
            return { points = t, slopes = slopes, smoothing_passes = pass - 1 }
        end
        local nextt = { { t[1][1], t[1][2] } }
        for i = 2, n - 1 do
            nextt[i] = {
                (t[i - 1][1] + t[i][1] + t[i + 1][1]) / 3,
                (t[i - 1][2] + t[i][2] + t[i + 1][2]) / 3,
            }
        end
        nextt[n] = { t[n][1], t[n][2] }
        t = nextt
    end
    error("luafemm B-H: smoothing did not converge")
end
--- Return H and dH/dB in SI for nonnegative B.
-- @tparam table curve Prepared cubic law.
-- @tparam number b Nonnegative flux-density magnitude in tesla.
-- @return H in A/m and its derivative in A/(m T).
function B.evaluate(curve, b)
    local t, s = curve.points, curve.slopes
    local n = #t
    if b > t[n][1] then
        return t[n][2] + s[n] * (b - t[n][1]), s[n]
    end
    for i = 1, n - 1 do
        if b <= t[i + 1][1] then
            local l = t[i + 1][1] - t[i][1]
            local z = (b - t[i][1]) / l
            local h = (1 - 3 * z * z + 2 * z ^ 3) * t[i][2]
                + z * (1 - z) ^ 2 * l * s[i]
                + z * z * (3 - 2 * z) * t[i + 1][2]
                + z * z * (z - 1) * l * s[i + 1]
            local d = 6 * z * (z - 1) * t[i][2] / l
                + (1 - 4 * z + 3 * z * z) * s[i]
                + 6 * z * (1 - z) * t[i + 1][2] / l
                + z * (3 * z - 2) * s[i + 1]
            return h, d
        end
    end
end
return B
