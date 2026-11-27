-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Reversible DC homogenisation of isotropic nonlinear iron and insulating sheets.
-- Across sheets B_n is continuous; along sheets H_t is continuous. The local
-- scalar solve minimises the constituent energy at fixed macroscopic B. All
-- returned laws/tangents are scaled by mu0, as required by the FE assembly.
local L = {}
local function check(v, s)
    assert(v, "luafemm lamination: " .. s)
end
local function finite(x)
    return type(x) == "number" and x == x and math.abs(x) < math.huge
end
local function frame(p)
    local a = math.rad(p.material_angle + (p.lamination == "parallel y" and 90 or 0))
    return math.cos(a), math.sin(a)
end
--- Return iron tangential induction and local law/tangent in sheet coordinates.
-- A bracketed Newton step is accepted only inside the monotone bracket.
-- Residual tolerance is relative to the nonzero imposed tangential induction;
-- neither a small field nor the saturation limit is replaced by zero.
local function state(p, bt, bn)
    local f, g = p.lamination_fill, 1 - p.lamination_fill
    local target = math.abs(bt)
    local function evaluate(q)
        local b = math.sqrt(q * q + bn * bn)
        local nu, d = require("luafemm").constitutive(p, b)
        local t = b > 0 and q / b or 0
        local n = b > 0 and bn / b or 0
        return f * q + g * nu * q - target,
            nu,
            nu + (d - nu) * t * t,
            (d - nu) * t * n,
            nu + (d - nu) * n * n
    end
    local lo, hi = 0, target / f
    check(finite(hi), "field/fill exceeds numerical range")
    local q = target
    if target == 0 then
        q = 0
    else
        local converged = false
        for _ = 1, 120 do
            local r, nu, tt = evaluate(q)
            if math.abs(r) <= 2e-14 * target then
                converged = true
                break
            end
            if r > 0 then
                hi = q
            else
                lo = q
            end
            local nextq = q - r / (f + g * tt)
            if not finite(nextq) or nextq <= lo or nextq >= hi then
                nextq = lo + (hi - lo) / 2
            end
            if nextq == q then
                break
            end
            q = nextq
        end
        local r = evaluate(q)
        check(
            converged or math.abs(r) <= 5e-13 * target,
            "local constitutive solve did not converge"
        )
    end
    q = bt < 0 and -q or q
    local b = math.sqrt(q * q + bn * bn)
    local nu, d = require("luafemm").constitutive(p, b)
    local t = b > 0 and q / b or 0
    local n = b > 0 and bn / b or 0
    return q, nu, nu + (d - nu) * t * t, (d - nu) * t * n, nu + (d - nu) * n * n
end
local function rotate(c, s, tt, tn, nn)
    return c * c * tt - 2 * c * s * tn + s * s * nn,
        c * s * (tt - nn) + (c * c - s * s) * tn,
        s * s * tt + 2 * c * s * tn + c * c * nn
end
--- Scaled H and consistent symmetric tangent in model coordinates.
function L.law(p, bx, by)
    check(finite(bx) and finite(by), "non-finite induction")
    local c, s = frame(p)
    local bt, bn = c * bx + s * by, -s * bx + c * by
    local q, nu, tt, tn, nn = state(p, bt, bn)
    local f, g = p.lamination_fill, 1 - p.lamination_fill
    local a = f + g * tt
    local ht, hn = nu * q, (f * nu + g) * bn
    local xx, xy, yy = rotate(c, s, tt / a, f * tn / a, f * nn + g - f * g * tn * tn / a)
    check(
        finite(ht) and finite(hn) and finite(xx) and finite(xy) and finite(yy),
        "constitutive range exceeded"
    )
    return c * ht - s * hn, s * ht + c * hn, xx, xy, yy
end
--- Energy of both constituents, including the insulating parallel path.
function L.energy(p, bx, by)
    local c, s = frame(p)
    local bt, bn = c * bx + s * by, -s * bx + c * by
    local q, nu = state(p, bt, bn)
    local f = p.lamination_fill
    local mu0 = require("luafemm").mu0
    return f * require("luafemm-post").primitive(p, math.sqrt(q * q + bn * bn))
        + (1 - f) * ((nu * q) ^ 2 + bn * bn) / (2 * mu0)
end
--- Reversible secant tensor reproducing H at the captured DC bias.
function L.secant(p, bx, by)
    local c, s = frame(p)
    local _, nu = state(p, c * bx + s * by, -s * bx + c * by)
    local f, g = p.lamination_fill, 1 - p.lamination_fill
    return rotate(c, s, nu / (f + g * nu), 0, f * nu + g)
end
return L
