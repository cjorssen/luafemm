-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Signed Galerkin restriction for equipotentials and translated periodic pairs.
-- No penalty parameter: V_i = sign_i * U_root + offset_i exactly. The optional
-- passive offset supports singularity subtraction without evaluating its kernel
-- at unconstrained singular mesh vertices. Constraint cycles may
-- force a zero root; incompatible voltage data are errors, including corners.
-- @module luafemm-electric-constraints
local C = {}
local function check(v, s)
    assert(v, "luafemm electric constraints: " .. s)
end
function C.solve(f, K, b, epsilon0, offset)
    local rhs = b
    if offset then
        rhs = {}
        for i, row in ipairs(K) do
            rhs[i] = b[i]
            for j, k in pairs(row) do
                rhs[i] = rhs[i] - k * (offset[j] or 0)
            end
        end
    end
    local parent, sign, zero = {}, {}, {}
    for i = 1, #f.nodes do
        parent[i], sign[i] = i, 1
    end
    local function find(i)
        if parent[i] ~= i then
            local r, s = find(parent[i])
            sign[i], parent[i] = sign[i] * s, r
        end
        return parent[i], sign[i]
    end
    local function join(a, b, s)
        local ra, sa = find(a)
        local rb, sb = find(b)
        if ra == rb then
            if sa ~= s * sb then
                zero[ra] = true
            end
        else
            parent[ra], sign[ra] = rb, s * sb * sa
            zero[rb] = zero[rb] or zero[ra]
        end
    end
    local representatives = {}
    for i, name in pairs(f.floating) do
        if representatives[name] then
            join(i, representatives[name], 1)
        else
            representatives[name] = i
        end
    end
    for _, p in ipairs(f.periodic_pairs or {}) do
        join(p[1], p[2], p[3])
    end
    local fixed, roots, count = {}, {}, 0
    for i = 1, #f.nodes do
        local r, s = find(i)
        if zero[r] then
            fixed[r] = 0
        end
        if f.fixed[i] ~= nil then
            local value = (f.fixed[i] - (offset and offset[i] or 0)) * s
            check(
                fixed[r] == nil
                    or math.abs(value - fixed[r]) <= 1e-10 * math.max(1, math.abs(value)),
                "incompatible potentials on a periodic pair or corner"
            )
            fixed[r] = value
        end
    end
    local anchored = next(fixed)
    for _, p in ipairs(f.periodic_pairs or {}) do
        anchored = anchored or p[3] == -1
    end
    for _, e in ipairs(f.edges) do
        anchored = anchored or (e.condition == "mixed" and e.coefficient > 0)
    end
    check(anchored, "at least one prescribed potential or positive mixed reference is required")
    for i = 1, #f.nodes do
        local r = find(i)
        if fixed[r] == nil and not roots[r] then
            count = count + 1
            roots[r] = count
        end
    end
    local kr, br = {}, {}
    for j = 1, count do
        kr[j], br[j] = {}, 0
    end
    for u = 1, #f.nodes do
        local ru, su = find(u)
        local j = roots[ru]
        if j then
            br[j] = br[j] + su * rhs[u]
            for v, a in pairs(K[u]) do
                local rv, sv = find(v)
                local k = roots[rv]
                if k then
                    kr[j][k] = (kr[j][k] or 0) + su * sv * a
                else
                    br[j] = br[j] - su * sv * a * fixed[rv]
                end
            end
        end
    end
    local charged_roots = {}
    for _, name in ipairs(f.electrodes) do
        local p = f.electrodes[name]
        if p.condition == "charge" then
            check(representatives[name], "charge-controlled electrode has no mesh nodes")
            local r, s = find(representatives[name])
            local j = roots[r]
            check(
                j and not charged_roots[r],
                "charge-controlled electrodes cannot be joined by a periodic constraint or reference"
            )
            charged_roots[r] = true
            br[j] = br[j] + s * p.charge / (epsilon0 * f.depth)
        end
    end
    local x, iterations, stats = {}, 0, { relative_residual = 0, backward_error = 0 }
    if count > 0 then
        x, iterations, stats =
            require("luafemm-cg").solve(kr, br, f.tolerance, { max_iterations = f.max_iterations })
    end
    local solution = {}
    for i = 1, #f.nodes do
        local r, s = find(i)
        solution[i] = s * (fixed[r] or x[roots[r]]) + (offset and offset[i] or 0)
    end
    return solution, iterations, stats
end
return C
