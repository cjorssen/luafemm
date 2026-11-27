-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Representability checks and peak exp(+i omega t) magnetic FEM conversions.
-- FEMM identifies wound labels by |turns|>1 or LamType>2. Never let that
-- implicit rule silently replace a resolved conductor with a homogenised one.
-- @module luafemm-harmonic-exchange
local X = {}
local function check(v, s)
    assert(v, "luafemm harmonic exchange: " .. s)
end
--- Convert a Cartesian peak phasor into the native amplitude/degree contract.
function X.source(re, im)
    return { current = math.sqrt(re * re + im * im), phase = math.deg(math.atan(im, re)) }
end
--- Translate supported passive FEM material fields into native declarations.
function X.import_material(p, q, response)
    local kind = p.lamtype or 0
    check(kind == 0 or kind >= 3 and kind <= 8, "harmonic directional laminations are unsupported")
    check(
        response or not p.bh or #p.bh == 0,
        "FEMM nonlinear harmonic laws are not the native fundamental law"
    )
    check(response or (p.h_c or 0) == 0, "harmonic permanent magnets require a separate DC bias")
    local thickness, fill = p.d_lam or 0, p.lamfill or 1
    q.lamination = kind == 0 and (thickness > 0 or fill ~= 1) and "stacked" or "none"
    q.lamination_thickness = thickness
    q.loss_angle_x, q.loss_angle_y = p.phi_hx or p.phi_h or 0, p.phi_hy or p.phi_h or 0
    check(
        (p.phi_h or 0) == 0 or (q.loss_angle_x == p.phi_h and q.loss_angle_y == p.phi_h),
        "conflicting isotropic/directional loss angles"
    )
    if kind > 2 then
        q.wire_model = ({
            [3] = "round",
            [4] = "stranded",
            [5] = "litz",
            [6] = "square",
            [7] = "cca10",
            [8] = "cca15",
        })[kind]
        q.wire_diameter, q.strands = p.wired or 0, p.nstrands or 1
        check(fill == 1 and thickness == 0, "wire records cannot also carry sheet data")
    else
        check(
            (p.nstrands or 0) == 0 and (p.wired or 0) == 0,
            "wire metadata require a supported wire LamType"
        )
        check(
            q.lamination ~= "stacked" or thickness > 0 or (p.sigma or 0) == 0,
            "conductive sheet records require positive thickness"
        )
    end
end
--- FEMM's documented wound/solid label distinction, used only in AC mode.
function X.conductor(p, label)
    return ((p.lamtype or 0) > 2 or math.abs(label[8]) > 1 or (p.sigma or 0) == 0) and "stranded"
        or "solid"
end
--- Reject nonrepresentable laws and controllers before constructing output.
function X.validate(m)
    if m.operating then
        require("luafemm-response-exchange").validate(m)
    end
    require("luafemm-harmonic").validate(m)
    for _, p in pairs(m.materials) do
        check(
            p.harmonic_law == "standard",
            "fundamental H(B) laws cannot be exported as FEMM B(H) laws"
        )
        check(
            p.lamination ~= "stacked" or p.lamination_thickness > 0 or p.conductivity == 0,
            "conductive sheet export requires positive thickness"
        )
        if require("luafemm-losses").active(p) then
            local t = require("luafemm-losses").material(p, m.options.frequency)
            for _, v in ipairs({ t.real, t.imaginary }) do
                check(
                    math.abs(v[2]) <= 1e-12 * math.max(math.abs(v[1]), math.abs(v[3]), 1e-30),
                    "rotated complex anisotropy cannot be represented by FEMM"
                )
            end
            check(
                p.material_angle % 90 == 0
                    or (p.intrinsic_mur_x == p.intrinsic_mur_y and p.loss_angle_x == p.loss_angle_y),
                "rotated intrinsic sheet laws cannot be represented by FEMM"
            )
        end
    end
    for _, c in pairs(m.circuits or {}) do
        if c.type ~= "series" or c.drive ~= "current" then
            check(
                m.stats and m.solution_signature == require("luafemm-cache").signature(m),
                "voltage/parallel export requires an unchanged solved model; output freezes branch currents"
            )
        end
    end
end
local function phasor(current, phase)
    phase = math.rad(phase or 0)
    return current * math.cos(phase), current * math.sin(phase)
end
--- Complete the common geometry document with complex physical records.
-- Voltage/parallel controllers become independent resolved series branches.
-- Material values stay intrinsic: FEMM performs its own wire/sheet reduction.
function X.document(m, d, t, labels, circuit_ids)
    d.header.comment =
        "luafemm peak exp(+i omega t); polygonal computational geometry; voltage/parallel controllers frozen to solved branch currents"
    for name, ids in pairs(circuit_ids) do
        local c = m.circuits[name]
        for winding, id in pairs(ids) do
            local re, im = phasor(c.current, c.phase)
            if c.drive ~= "current" or c.type ~= "series" then
                local value = require("luafemm-ac-circuits").result(m, name, winding)
                re, im = value.current_real, value.current_imaginary
            end
            d.circuits[id].totalamps_re, d.circuits[id].totalamps_im = re, im
        end
    end
    for _, f in ipairs(t.faces) do
        if not f.hole then
            local id = labels[f.id] + 1
            local label, p = d.labels[id], m.materials[f.material]
            local r, q = f.region and m.regions[f.region], d.materials[id]
            local wire = p.wire_model ~= "none"
            local sheet = p.lamination == "stacked"
            if wire or sheet or require("luafemm-losses").active(p) then
                local swap = p.material_angle % 180 == 90
                q.mu_x, q.mu_y = p.intrinsic_mur_x, p.intrinsic_mur_y
                q.phi_hx, q.phi_hy = p.loss_angle_x, p.loss_angle_y
                if swap then
                    q.mu_x, q.mu_y, q.phi_hx, q.phi_hy = q.mu_y, q.mu_x, q.phi_hy, q.phi_hx
                end
            end
            if wire then
                q.lamtype = (
                    { round = 3, stranded = 4, litz = 5, square = 6, cca10 = 7, cca15 = 8 }
                )[p.wire_model]
                q.wired, q.nstrands, q.sigma = p.wire_diameter, p.strands, p.conductivity / 1e6
            elseif sheet then
                q.d_lam, q.lamfill = p.lamination_thickness, p.lamination_fill
                q.sigma = p.conductivity / 1e6
            end
            if r and r.current_constraint == "zero" then
                check(label[5] == 0, "conflicting constrained circuit")
                d.circuits[#d.circuits + 1] = {
                    circuitname = "zero-current region " .. id,
                    totalamps_re = 0,
                    totalamps_im = 0,
                    circuittype = 1,
                }
                label[5], label[8], q.j_re = #d.circuits, 1, 0
            end
            if label[5] ~= 0 then
                local stranded = r.conductor ~= "solid" and r.current_constraint ~= "zero"
                check(
                    not stranded or wire or math.abs(label[8]) > 1 or p.conductivity == 0,
                    "FEMM would resolve this single-turn stranded conductor; use a physical wire model"
                )
                check(
                    stranded or math.abs(label[8]) == 1,
                    "split solid regions cannot be exported as separate fractional-turn FEMM circuits"
                )
            else
                check(
                    not r or r.conductor ~= "stranded" or p.conductivity == 0,
                    "stranded conducting regions require a named circuit for FEMM exchange"
                )
                q.j_re, q.j_im = phasor(f.current / 1e6, r and r.current_phase or 0)
            end
        end
    end
end
--- Per-label answer state. FEMM case 0 stores minus the voltage gradient;
-- case 1 stores additional uniform J in MA/m2. Axisymmetric gradients use V/r.
function X.answer_state(m, face, label)
    if label[5] == 0 then
        return 1, 0, 0
    end
    local r = m.regions[face.region]
    local entry = m.ac_network and m.ac_network.by_region[face.region]
    if entry and entry.solid then
        local k = entry.solid.index
        local length = m.options.problem_type == "axisymmetric" and 2 * math.pi
            or m.options.depth * m.options.unit
        return 0, -m.ac_network.real[k] / length, -m.ac_network.imaginary[k] / length
    end
    local re, im
    if entry then
        re, im = entry.branch.ir, entry.branch.ii
    else
        re, im = phasor(m.circuits[r.circuit].current, m.circuits[r.circuit].phase)
    end
    local density = r.turns / m.static_circuits[r.circuit].areas[face.region] / 1e6
    return 1, density * re, density * im
end
return X
