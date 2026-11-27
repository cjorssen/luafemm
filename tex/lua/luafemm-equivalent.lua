-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Passive nonlinear equivalent-amplitude harmonic material law.
-- The amplitude is sqrt(|Bx|^2+|By|^2), i.e. sqrt(2) times vector RMS.
-- Project the reversible H(B) curve for a linearly polarised sinusoid, then
-- multiply reluctivity by exp(i delta). Elliptical waveforms use this explicitly
-- equivalent amplitude; this is not harmonic balance or hysteresis history.
local N = {}
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm equivalent harmonic: " .. s)
end
function N.active(m)
    for _, p in pairs(m.materials) do
        if p.harmonic_law == "equivalent" then
            return true
        end
    end
    return false
end
function N.validate(m)
    if not N.active(m) then
        return
    end
    check(
        m.options.frequency > 0
            and not m.operating
            and not m.transient
            and not m.sliding
            and m.options.problem_type ~= "axial",
        "requires an unbiased fixed in-plane AC model"
    )
    check(
        not require("luafemm-fundamental").active(m),
        "cannot mix real fundamental and equivalent laws"
    )
    check(
        not require("luafemm-ac-circuits").active(m)
            and not m.network
            and not require("luafemm-flux-drive").active(m),
        "requires prescribed stranded current sources"
    )
    for _, p in pairs(m.materials) do
        check(
            not p.bh or p.harmonic_law == "equivalent",
            "nonlinear materials must select equivalent law"
        )
        check(
            p.conductivity == 0 and p.wire_model == "none" and p.hc == 0,
            "equivalent law excludes diffusion, wire homogenisation and permanent magnets"
        )
    end
end
function N.material(p, b)
    local h, d = require("luafemm-fundamental").scalar(p, b)
    local nu = mu0 * (b > 1e-15 and h / b or d)
    local angle = math.rad(p.loss_angle_x)
    local nr, ni = nu * math.cos(angle), nu * math.sin(angle)
    return {
        real = { nr, 0, nr },
        imaginary = { ni, 0, ni },
        hysteresis = { ni, 0, ni },
        lamination = { 0, 0, 0 },
        iteration_ratio = b > 1e-15 and d * b / math.max(h, 1e-300) or 1,
    }
end
--- Element-centroid equivalent law; the ordered nodal arrays own the state.
function N.element(m, e)
    local bx, by, ix, iy = 0, 0, 0, 0
    if m.A_im then
        if m.options.problem_type == "axisymmetric" then
            local x, y = 0, 0
            for _, id in ipairs(e.ids) do
                x = x + m.nodes[id][1] / 3
                y = y + m.nodes[id][2] / 3
            end
            bx, by = require("luafemm-axisymmetric").field(m, e, x, y)
            local old = m.A
            m.A = m.A_im
            ix, iy = require("luafemm-axisymmetric").field(m, e, x, y)
            m.A = old
        else
            for j, id in ipairs(e.ids) do
                bx, by = bx + e.gy[j] * m.A[id], by - e.gx[j] * m.A[id]
                ix, iy = ix + e.gy[j] * m.A_im[id], iy - e.gx[j] * m.A_im[id]
            end
        end
    end
    return N.material(e.p, math.sqrt(bx * bx + by * by + ix * ix + iy * iy))
end
--- Damped coefficient iteration with an independently measured fixed-point
-- defect. Inner linear residual and nonlinear defect must both meet tolerance.
function N.solve(m, assemble)
    N.validate(m)
    local H = require("luafemm-harmonic")
    local cache = m.options.cache
    local tensors = {}
    for i, e in ipairs(m.triangles) do
        if e.p.harmonic_law == "equivalent" then
            tensors[i] = N.element(m, e)
        end
    end
    m.equivalent_tensors = tensors
    m.options.cache = { mode = "off" }
    local ok, result = pcall(function()
        for iteration = 1, m.options.max_newton do
            -- The linear harmonic assembler expects zero free-node initial data.
            for i, p in ipairs(m.nodes) do
                if p.dof then
                    m.A[i] = 0
                end
            end
            m.A_im = nil
            H.solve(m, assemble)
            local defect, ratio = 0, 1
            local next_tensors = {}
            for i, t in pairs(tensors) do
                local fresh = N.element(m, m.triangles[i])
                defect = math.max(defect, math.abs(fresh.real[1] - t.real[1]) / fresh.real[1])
                ratio = math.max(ratio, fresh.iteration_ratio)
                next_tensors[i] = fresh
            end
            m.stats.nonlinear_iterations, m.stats.nonlinear_defect = iteration, defect
            if m.options.log then
                m.options.log(
                    string.format("Equivalent iteration %d: defect %.3g", iteration, defect)
                )
            end
            if defect <= m.options.tolerance then
                return m
            end
            for i, t in pairs(tensors) do
                local fresh = next_tensors[i]
                for _, key in ipairs({ "real", "imaginary", "hysteresis" }) do
                    local relaxation = math.min(0.5, 1 / ratio)
                    for j = 1, 3 do
                        t[key][j] = (1 - relaxation) * t[key][j] + relaxation * fresh[key][j]
                    end
                end
            end
        end
        error("luafemm equivalent harmonic: nonlinear iteration did not converge")
    end)
    m.options.cache, m.equivalent_tensors = cache, nil
    if not ok then
        m.stats = nil
        error(result, 0)
    end
    require("luafemm").update_fields(m)
    require("luafemm-cache").store(m)
    return result
end
return N
