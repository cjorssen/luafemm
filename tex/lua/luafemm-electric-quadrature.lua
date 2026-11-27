-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Bounded adaptive vector integration on a scalar interval.
-- Positive four-point Gauss rules compare a parent with two children. Each
-- component has its own absolute-integrand scale: units and cancellation do
-- not mix in a vector norm. The shared budget counts rule evaluations.
local Q = {}
local nodes = {
    { 0.06943184420297371, 0.17392742256872693 },
    { 0.33000947820757187, 0.32607257743127307 },
    { 0.66999052179242813, 0.32607257743127307 },
    { 0.93056815579702629, 0.17392742256872693 },
}
function Q.new(tolerance, depth, budget, label)
    local count = 0
    label = label or "electric quadrature"
    local function rule(fn, a, b)
        count = count + 1
        assert(count <= budget, "luafemm " .. label .. ": quadrature cell budget exceeded")
        local values, scale = {}, {}
        for _, v in ipairs(nodes) do
            local y = fn(a + (b - a) * v[1])
            for i, x in ipairs(y) do
                assert(
                    x == x and math.abs(x) < math.huge,
                    "luafemm " .. label .. ": nonfinite integrand"
                )
                local z = (b - a) * v[2] * x
                values[i] = (values[i] or 0) + z
                scale[i] = (scale[i] or 0) + math.abs(z)
            end
        end
        return values, scale
    end
    local function refine(fn, a, b, d, coarse)
        coarse = coarse or rule(fn, a, b)
        local mid = (a + b) / 2
        local left, ls = rule(fn, a, mid)
        local right, rs = rule(fn, mid, b)
        local result, ok = {}, true
        for i, x in ipairs(coarse) do
            result[i] = left[i] + right[i]
            ok = ok and math.abs(result[i] - x) <= tolerance * math.max(ls[i] + rs[i], 1e-300)
        end
        if ok then
            return result
        end
        assert(
            d < depth,
            "luafemm " .. label .. ": quadrature tolerance not reached within maximum depth"
        )
        left = refine(fn, a, mid, d + 1, left)
        right = refine(fn, mid, b, d + 1, right)
        for i in ipairs(result) do
            result[i] = left[i] + right[i]
        end
        return result
    end
    return function(fn, a, b)
        return refine(fn, a or 0, b or 1, 0)
    end, function()
        return count
    end
end
return Q
