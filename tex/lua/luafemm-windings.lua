-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Homogenised round-wire windings with peak exp(+i omega t) phasors.
-- The hexagonal-array fits follow Meeker's effective wire model as used in
-- FEMM/xfemm GetFillFactor. This independent implementation stores resistivity
-- (ohm metre), avoiding ambiguous effective conductivity units. The local
-- energy correction is essential: macroscopic A already includes part of it.
-- @module luafemm-windings
local W = {}
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm windings: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
--- Validate geometry-independent material data; diameter is always in mm.
-- @tparam table p Owned material with conductivity and permeability prepared.
function W.configure(p)
    p.wire_model = p.wire_model or "none"
    p.strands = p.strands or 1
    p.wire_diameter = p.wire_diameter or 0
    check(
        p.wire_model == "none"
            or p.wire_model == "round"
            or p.wire_model == "stranded"
            or p.wire_model == "litz"
            or p.wire_model == "square"
            or p.wire_model == "cca10"
            or p.wire_model == "cca15",
        "unknown wire model"
    )
    check(
        finite(p.strands) and p.strands >= 1 and p.strands % 1 == 0,
        "strands must be a positive integer"
    )
    check(finite(p.wire_diameter) and p.wire_diameter >= 0, "invalid wire diameter")
    if p.wire_model == "none" then
        check(p.wire_diameter == 0 and p.strands == 1, "wire data require a wire model")
        return
    end
    check(p.wire_diameter > 0, "wire diameter must be positive")
    check(p.wire_model ~= "round" or p.strands == 1, "round wire requires strands=1")
    check(p.wire_model ~= "square" or p.strands == 1, "square wire requires strands=1")
    check(
        p.conductivity > 0 and p.mur_x == 1 and p.mur_y == 1 and not p.bh and p.hc == 0,
        "wire models require a linear nonmagnetic conductor"
    )
    check(
        not require("luafemm-losses").active(p),
        "wire and magnetic loss models cannot be combined"
    )
end
--- Check circuit ownership before cache lookup; fill is derived, never guessed.
-- @tparam table m Physical model before cache lookup or meshing.
function W.validate(m)
    local used = {}
    for _, r in ipairs(m.regions) do
        local p = m.materials[r.material]
        if p.wire_model ~= "none" then
            used[p] = true
            check(m.options.frequency > 0, "wire homogenisation requires frequency>0")
            check(
                not m.operating and m.options.problem_type ~= "axial",
                "wire homogenisation requires an unbiased planar or axisymmetric model"
            )
            check(
                r.circuit and r.conductor ~= "solid" and r.current_constraint ~= "zero",
                "wire regions require a named stranded circuit"
            )
            check(
                r.fill == 1 and r.length == 0,
                "wire regions derive winding fill and length from geometry"
            )
        end
    end
    for _, p in pairs(m.materials) do
        check(p.wire_model == "none" or used[p], "wire material requires a circuit region")
    end
end
--- Return effective relative reluctivity, resistivity and geometric copper fill.
-- radius is the strand radius except for an untransposed bundle, whose equal-
-- area radius is used. This model does not resolve individual strand currents.
-- @tparam table p Configured material.
-- @tparam number frequency Positive frequency in Hz.
-- @tparam number fill Copper area / block area, in (0, pi/sqrt(12)].
function W.effective(p, frequency, fill)
    check(finite(frequency) and frequency > 0, "positive frequency required")
    check(
        finite(fill)
            and fill > 0
            and fill <= (p.wire_model == "square" and 1 or math.pi / math.sqrt(12)),
        "wire fill exceeds hexagonal packing or is nonpositive"
    )
    if p.wire_model == "square" then
        -- A square array is represented by an equivalent insulated foil stack.
        -- Use one consistent diffusion frequency in PDE and terminal properties.
        local fraction = math.sqrt(fill)
        local d = p.wire_diameter * 0.001
        local pitch, omega = d / fraction, 2 * math.pi * frequency
        local w = omega * p.conductivity * mu0 * d * d * fill / 4
        local fr, fi = require("luafemm-losses").sheet(math.sqrt(w / 2), math.sqrt(w / 2))
        local mr, mi = 1 - fraction + fraction * fr, fraction * fi
        local den = mr * mr + mi * mi
        local dc = 1 / (p.conductivity * fill)
        local rr = dc * fr / (fr * fr + fi * fi) + omega * mu0 * pitch * pitch * mi / 12
        local ri = -dc * fi / (fr * fr + fi * fi)
            + omega * mu0 * pitch * pitch * ((1 - fraction) / 4 - mr / 12)
        check(rr > 0 and finite(rr) and finite(ri), "invalid equivalent-foil coefficients")
        return {
            real = { mr / den, 0, mr / den },
            imaginary = { -mi / den, 0, -mi / den },
            hysteresis = { 0, 0, 0 },
            lamination = { 0, 0, 0 },
            resistivity_real = rr,
            resistivity_imaginary = ri,
            dc_resistivity = dc,
            fill = fill,
        }
    end
    local radius = p.wire_diameter * 0.0005
    if p.wire_model == "stranded" then
        radius = radius * math.sqrt(p.strands)
    end
    local omega = 2 * math.pi * frequency
    local w = omega * p.conductivity * mu0 * radius ^ 2 / 2
    local c1 = 0.7756067409818643
        + fill
            * (0.6873854335408803 + fill * (0.06841584481674128 - 0.07143732702512284 * fill))
    local c2 = 1.5 * fill / c1
    local c3 = 0.8824642871525136
        + fill
            * (-0.008605512994838827 + fill * (0.7223208744682307 - 0.2157183942377177 * fill))
    local c4 = math.log(1.5299240194394943 / math.sqrt(fill)) - c3 / 3
    -- FEMM 5's published CCA fits use the bulk DC conductivity and a copper
    -- cross-section fraction of 10% or 15%. Coefficients describe microscopic
    -- response, not a replacement of aluminium by homogeneous copper.
    if p.wire_model == "cca10" or p.wire_model == "cca15" then
        check(fill >= 0.05, "CCA fits require fill >= 0.05")
        local a = p.wire_model == "cca10"
                and {
                    {
                        0.7270741505617485,
                        0.8902950067721367,
                        0.11894736885885195,
                        -0.12247276254503957,
                    },
                    {
                        0.006784920229549677,
                        1.8942880489198526,
                        -1.3631438759519217,
                        0.504431701685587,
                    },
                    {
                        0.7879151117538901,
                        -0.7216375417404758,
                        2.266562622770727,
                        -1.1349603497265566,
                    },
                    {
                        1.6190358786364027,
                        -3.967461431323324,
                        4.068458362148543,
                        -1.7975534375694646,
                    },
                }
            or {
                {
                    0.7486913529860821,
                    0.9042845510838825,
                    0.1361040321433224,
                    -0.10652380745682069,
                },
                {
                    0.006790468527313965,
                    1.8945509985370095,
                    -1.3643501010185972,
                    0.5036765577982594,
                },
                {
                    0.7519323306379909,
                    -0.7184130853903536,
                    2.2551599408130594,
                    -1.1439852083597326,
                },
                { 1.6204876714534728, -3.966250827502031, 4.064413095233383, -1.7940297437703776 },
            }
        local function poly(t)
            return t[1] + fill * (t[2] + fill * (t[3] + fill * t[4]))
        end
        c1, c2, c3, c4 = poly(a[1]), poly(a[2]), poly(a[3]), poly(a[4])
    end
    local L = require("luafemm-losses")
    local z = math.sqrt(c1 * w / 2)
    local fr, fi = L.sheet(z, z)
    local mr, mi = 1 - c2 + c2 * fr, c2 * fi
    local den = mr * mr + mi * mi
    local nr, ni = mr / den, -mi / den
    z = math.sqrt(c3 * w / 2)
    fr, fi = L.sheet(z, z)
    den = fr * fr + fi * fi
    local dc = 1 / (p.conductivity * fill)
    local pitch = 1.6494541661869013 * radius / math.sqrt(fill)
    local correction = omega * mu0 * pitch ^ 2 / 12
    local rr = dc * fr / den + correction * mi
    local ri = dc * (c4 * w - fi / den) - correction * mr
    check(
        finite(rr) and finite(ri) and finite(nr) and finite(ni) and rr > 0 and ni >= 0,
        "non-passive or non-finite wire coefficients"
    )
    return {
        real = { nr, 0, nr },
        imaginary = { ni, 0, ni },
        hysteresis = { 0, 0, 0 },
        lamination = { 0, 0, 0 },
        resistivity_real = rr,
        resistivity_imaginary = ri,
        dc_resistivity = dc,
        fill = fill,
    }
end
--- Rebuild per-region coefficients from the meshed circuit areas, also on restore.
-- @tparam table m Meshed model with static circuit areas and moments.
function W.prepare(m)
    m.wire_regions = {}
    for id, r in ipairs(m.regions) do
        local p = m.materials[r.material]
        if p.wire_model ~= "none" then
            local area = m.static_circuits[r.circuit].areas[id]
            local fill = math.abs(r.turns)
                * p.strands
                * (p.wire_model == "square" and 4 or math.pi)
                * (p.wire_diameter * 0.0005) ^ 2
                / area
            m.wire_regions[id] = W.effective(p, m.options.frequency, fill)
        end
    end
end
--- Microscopic transport/proximity dissipation and internal reactive energy.
-- J is the macroscopic winding current density; B is the homogenised flux.
-- The returned signed energy correction must be added exactly once.
-- @tparam table m Solved model supplying positive frequency in Hz.
-- @tparam table e Element, optionally carrying effective wire coefficients.
-- @tparam table f SI peak J, Bx and By phasor pairs.
-- @return Transport W/m^3, proximity W/m^3 and local energy J/m^3.
function W.density(m, e, f)
    if not e.wire then
        return 0, 0, 0
    end
    local t = e.wire
    local j2 = f.J[1] ^ 2 + f.J[2] ^ 2
    local b2 = f.Bx[1] ^ 2 + f.Bx[2] ^ 2 + f.By[1] ^ 2 + f.By[2] ^ 2
    local omega = 2 * math.pi * m.options.frequency
    return t.resistivity_real * j2 / 2,
        omega * t.imaginary[1] * b2 / (2 * mu0),
        t.resistivity_imaginary * j2 / (4 * omega)
end
return W
