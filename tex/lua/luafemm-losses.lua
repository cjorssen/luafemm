-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Passive linear complex permeability and insulated-sheet homogenisation.
-- Peak phasors use exp(+i omega t): passive permeability has negative imaginary
-- part and reluctivity has positive imaginary part. Thickness is in millimetres,
-- independently of model units. No material is mutated during frequency sweeps.
-- @module luafemm-losses
local L = {}
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm losses: " .. s)
end
local function number(v, default, name)
    if v == nil then
        v = default
    end
    check(type(v) == "number" and v == v and math.abs(v) < math.huge, "invalid " .. name)
    return v
end
--- Validate material inputs before preparing static or harmonic tensors.
-- @tparam table p Owned material, with lamination and intrinsic permeabilities.
function L.configure(p)
    p.loss_angle = number(p.loss_angle, 0, "loss angle")
    p.loss_angle_x = number(p.loss_angle_x, p.loss_angle, "loss angle x")
    p.loss_angle_y = number(p.loss_angle_y, p.loss_angle, "loss angle y")
    for _, a in ipairs({ p.loss_angle, p.loss_angle_x, p.loss_angle_y }) do
        check(a >= 0 and a < 90, "loss angles must be in [0,90) degrees")
    end
    p.lamination_thickness = number(p.lamination_thickness, 0, "lamination thickness")
    check(p.lamination_thickness >= 0, "lamination thickness must be nonnegative")
    check(
        p.lamination_thickness == 0 or p.lamination == "stacked",
        "thickness requires lamination=stacked"
    )
end
--- Whether frequency-dependent material coefficients are requested.
-- @tparam table p Configured material.
-- @treturn boolean Whether to prepare an effective harmonic tensor.
function L.active(p)
    return p.loss_angle_x > 0 or p.loss_angle_y > 0 or p.lamination == "stacked"
end
--- Reject unsupported constitutive combinations before cache reuse.
-- @tparam table m Physical model before meshing.
-- @tparam table p Configured material; not mutated.
function L.validate(m, p)
    check(
        p.lamination == "none" or p.lamination == "stacked",
        "harmonic laminations require lamination=stacked"
    )
    if L.active(p) then
        check(
            not p.bh or p.harmonic_law == "equivalent",
            "loss models require a linear material or equivalent law"
        )
        check(not m.operating, "loss models on captured operating points are not implemented")
        check(
            m.options.problem_type ~= "axial",
            "axial diffusion does not support complex material laws"
        )
    end
    if p.library then
        local f = p.library.fields
        check(
            ((f.LamType or 0) == 0 and (f.NStrands or 0) == 0 and (f.WireD or 0) == 0)
                or ((f.LamType >= 3 and f.LamType <= 8) and p.wire_model ~= "none"),
            "catalogue lamination/strand models are unsupported in harmonic mode"
        )
        check((f.J_im or 0) == 0, "unsupported catalogue harmonic property J_im")
        check(
            p.library_harmonic
                or ((f.Phi_h or 0) == 0 and (f.Phi_hx or 0) == 0 and (f.Phi_hy or 0) == 0),
            "catalogue loss metadata require a supported linear conversion"
        )
        -- Native properties must carry any nonzero catalogue loss data. A raw
        -- library table is provenance, not an alternative hidden material law.
        check(
            (f.d_lam or 0) == 0 or p.lamination == "stacked" and p.lamination_thickness > 0,
            "catalogue lamination data require an explicit supported sheet model"
        )
    end
end
local function div(ar, ai, br, bi)
    local scale = math.max(math.abs(br), math.abs(bi))
    check(scale > 0 and scale < math.huge, "non-finite effective permeability")
    br, bi = br / scale, bi / scale
    local d = br * br + bi * bi
    return (ar / scale * br + ai / scale * bi) / d, (ai / scale * br - ar / scale * bi) / d
end
local function sinc(x)
    return math.abs(x) < 1e-4 and 1 - x * x / 6 + x ^ 4 / 120 or math.sin(x) / x
end
--- Return tanh(z)/z and the mean squared microscopic H/H_surface. Exponentials
-- only have nonpositive real arguments; the small-z series avoids cancellation.
-- @tparam number x Nonnegative real part of z.
-- @tparam number y Imaginary part of z.
-- @return Real and imaginary transfer, and mean-square microscopic H ratio.
function L.sheet(x, y)
    local t = x * x + y * y
    if t < 1e-8 then
        local ar, ai = x * x - y * y, 2 * x * y
        return 1 - ar / 3 + 2 * (ar * ar - ai * ai) / 15,
            -ai / 3 + 4 * ar * ai / 15,
            1 - 2 * ar / 3 + (6 * ar * ar - 2 * ai * ai) / 15
    end
    local e = math.exp(-2 * x)
    local qr, qi = e * math.cos(2 * y), -e * math.sin(2 * y)
    local tr, ti = div(1 - qr, -qi, 1 + qr, qi)
    local fr, fi = div(tr, ti, x, y)
    local q
    if x > 20 then
        q = ((1 - e * e) / (2 * x) + 2 * e * sinc(2 * y)) / (1 + e * e + 2 * e * math.cos(2 * y))
    else
        local sh = math.abs(x) < 1e-4 and 1 + (2 * x) ^ 2 / 6 + (2 * x) ^ 4 / 120
            or math.sinh(2 * x) / (2 * x)
        q = (sh + sinc(2 * y)) / (math.cosh(2 * x) + math.cos(2 * y))
    end
    return fr, fi, q
end
local function rotate(x, y, angle)
    local c, s = math.cos(math.rad(angle)), math.sin(math.rad(angle))
    return { c * c * x + s * s * y, c * s * (x - y), s * s * x + c * c * y }
end
--- Effective relative reluctivity and loss tensors at one positive frequency.
-- @tparam table p Prepared linear material (read-only).
-- @tparam number frequency Hertz.
-- @treturn table Rotated real/imaginary reluctivity and intrinsic/sheet loss tensors.
function L.material(p, frequency)
    local nr, ni, hy, ed = {}, {}, {}, {}
    local mus = { p.intrinsic_mur_x, p.intrinsic_mur_y }
    local angles = { p.loss_angle_x, p.loss_angle_y }
    local fill = p.lamination == "stacked" and p.lamination_fill or 1
    for i = 1, 2 do
        local a = math.rad(angles[i])
        local mr, mi = mus[i] * math.cos(a), -mus[i] * math.sin(a)
        local fr, fi, q = 1, 0, 1
        if p.lamination == "stacked" and p.lamination_thickness > 0 and p.conductivity > 0 then
            local k = p.lamination_thickness
                * 0.001
                / 2
                * math.sqrt(2 * math.pi * frequency * mu0 * mus[i] * p.conductivity)
            fr, fi, q =
                L.sheet(k * math.cos(math.pi / 4 - a / 2), k * math.sin(math.pi / 4 - a / 2))
        end
        local er, ei = 1 - fill + fill * (mr * fr - mi * fi), fill * (mr * fi + mi * fr)
        nr[i], ni[i] = div(1, 0, er, ei)
        local magnitude = nr[i] * nr[i] + ni[i] * ni[i]
        hy[i] = fill * -mi * q * magnitude
        if p.lamination ~= "stacked" or p.lamination_thickness == 0 or p.conductivity == 0 then
            hy[i] = ni[i]
        end
        ed[i] = ni[i] - hy[i]
        check(ed[i] >= -1e-10 * math.max(ni[i], 1e-30), "non-passive sheet loss")
        ed[i] = math.max(0, ed[i])
        check(nr[i] > 0 and ni[i] >= 0, "non-passive effective permeability")
    end
    return {
        real = rotate(nr[1], nr[2], p.material_angle),
        imaginary = rotate(ni[1], ni[2], p.material_angle),
        hysteresis = rotate(hy[1], hy[2], p.material_angle),
        lamination = rotate(ed[1], ed[2], p.material_angle),
    }
end
--- Apply a symmetric tensor {xx,xy,yy} to a real vector.
-- @tparam table t Three independent tensor coefficients.
-- @tparam number x First vector component.
-- @tparam number y Second vector component.
-- @return Two transformed components.
function L.apply(t, x, y)
    return t[1] * x + t[2] * y, t[2] * x + t[3] * y
end
--- Cycle-average intrinsic and microscopic sheet losses per volume (W/m^3).
-- @tparam table m Model supplying the frequency in hertz.
-- @tparam table e Element with a prepared ac_material tensor, or a lossless element.
-- @tparam table f Point phasors Bx/By in tesla.
-- @return Intrinsic magnetic loss density and microscopic sheet Joule loss density.
function L.density(m, e, f)
    if not e.ac_material then
        return 0, 0
    end
    local h, l = 0, 0
    for part = 1, 2 do
        local x, y = f.Bx[part], f.By[part]
        local hx, hy = L.apply(e.ac_material.hysteresis, x, y)
        local lx, ly = L.apply(e.ac_material.lamination, x, y)
        h = h + x * hx + y * hy
        l = l + x * lx + y * ly
    end
    local w = math.pi * m.options.frequency / mu0
    return w * math.max(0, h), w * math.max(0, l)
end
return L
