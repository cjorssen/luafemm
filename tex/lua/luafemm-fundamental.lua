-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Lossless first-harmonic projection of an isotropic piecewise-linear H(B).
-- Assume B(t)=Bhat*cos(omega*t), retain the fundamental of H, and solve that
-- real nonlinear spatial problem. This is not transient hysteresis, nor FEMM's
-- alternative B(H) projection. Exact segment integrals avoid angular sampling
-- noise in Newton tangents. No physical material tables are mutated.
-- @module luafemm-fundamental
local F = {}
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm fundamental: " .. s)
end
--- Validate the explicit material-law selector after the static B-H setup.
-- @tparam table p Owned material; the original B-H table remains unchanged.
function F.configure(p)
    p.harmonic_law = p.harmonic_law or "standard"
    check(
        p.harmonic_law == "standard"
            or p.harmonic_law == "fundamental"
            or p.harmonic_law == "equivalent",
        "unknown harmonic law"
    )
    if p.harmonic_law == "equivalent" then
        check(
            p.bh
                and p.interpolation == "linear"
                and p.hc == 0
                and p.lamination == "none"
                and p.loss_angle_x == p.loss_angle_y,
            "equivalent law requires isotropic piecewise-linear B-H data and a common loss angle"
        )
    end
    if p.harmonic_law == "fundamental" then
        check(
            p.bh and p.interpolation == "linear",
            "fundamental law requires a linear-interpolated B-H table"
        )
        check(
            p.hc == 0 and not require("luafemm-losses").active(p),
            "fundamental law requires a lossless material without permanent magnetisation"
        )
    end
end
--- Whether a model requests any lossless nonlinear fundamental material.
-- @tparam table m Model containing configured materials.
-- @treturn boolean Whether to use the real nonlinear harmonic solve.
function F.active(m)
    for _, p in pairs(m.materials) do
        if p.harmonic_law == "fundamental" then
            return true
        end
    end
    return false
end
--- Reject combinations whose spatial phasors need a complex nonlinear solver.
-- @tparam table m Unmeshed model, checked before cache reuse.
function F.validate(m)
    if not F.active(m) then
        return
    end
    check(m.options.frequency > 0, "fundamental law requires frequency>0")
    check(
        not m.operating and m.options.problem_type ~= "axial",
        "fundamental law requires an unbiased in-plane model"
    )
    check(
        not require("luafemm-ac-circuits").active(m),
        "fundamental law requires prescribed series currents without wire homogenisation"
    )
    for _, p in pairs(m.materials) do
        check(not p.nonlinear_sheet, "fundamental law excludes nonlinear directional sheets")
        check(
            not require("luafemm-losses").active(p),
            "fundamental law excludes complex magnetic losses"
        )
    end
    for _, c in pairs(m.circuits or {}) do
        check(c.phase == 0, "fundamental law requires phase=0 excitations")
    end
    check(m.materials.air.conductivity == 0, "fundamental law excludes conductive background air")
    for _, r in ipairs(m.regions) do
        check((r.current_phase or 0) == 0, "fundamental law requires phase=0 excitations")
        check(
            m.materials[r.material].conductivity == 0 or r.circuit or r.conductor == "stranded",
            "fundamental law excludes resolved eddy currents"
        )
    end
end
--- Fundamental H amplitude, its tangent, and average stored energy (SI).
-- Integrate each linear H=a*B+c segment over the corresponding cosine interval.
-- The continuation beyond the table uses vacuum slope, like the static law.
-- @tparam table p Isotropic material with a validated piecewise-linear B-H table.
-- @tparam number b Nonnegative peak flux-density magnitude in tesla.
-- @return H1 in A/m, derivative in A/(m T), and mean energy density in J/m^3.
function F.scalar(p, b)
    local points = p.bh
    if b <= 1e-15 then
        return 0, points[2][2] / points[2][1], 0
    end
    local h, tangent, energy, previous_energy = 0, 0, 0, 0
    for i = 1, #points do
        local lo, hlo = points[i][1], points[i][2]
        local hi = points[i + 1] and points[i + 1][1] or math.huge
        local slope = points[i + 1] and (points[i + 1][2] - hlo) / (hi - lo) or 1 / mu0
        local offset = hlo - slope * lo
        local constant = previous_energy - slope * lo * lo / 2 - offset * lo
        if b > lo then
            local first, last = math.acos(math.min(hi / b, 1)), math.acos(lo / b)
            local c1 = math.sin(last) - math.sin(first)
            local c2 = (last - first) / 2 + (math.sin(2 * last) - math.sin(2 * first)) / 4
            h = h + slope * b * c2 + offset * c1
            tangent = tangent + slope * c2
            energy = energy + slope * b * b * c2 / 2 + offset * b * c1 + constant * (last - first)
        end
        if hi >= b then
            break
        end
        previous_energy = previous_energy + (hlo + points[i + 1][2]) * (hi - lo) / 2
    end
    return 4 * h / math.pi, 4 * tangent / math.pi, 2 * energy / math.pi
end
--- Scaled vector law and consistent symmetric tangent for real Newton assembly.
-- @tparam table p Configured nonlinear material.
-- @tparam number bx Real peak Bx, in tesla.
-- @tparam number by Real peak By, in tesla.
-- @return mu0*Hx, mu0*Hy, tangent xx, xy, yy.
function F.law(p, bx, by)
    local b = math.sqrt(bx * bx + by * by)
    local h, d = F.scalar(p, b)
    local nu = b > 1e-15 and mu0 * h / b or mu0 * d
    local k = b > 1e-15 and (mu0 * d - nu) / (b * b) or 0
    return nu * bx, nu * by, nu + k * bx * bx, k * bx * by, nu + k * by * by
end
return F
