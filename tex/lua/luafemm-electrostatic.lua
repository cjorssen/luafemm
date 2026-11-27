-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Electric model construction, ownership and provider dispatch.
-- Coordinates are in model units; charges are C, epsilon is F/m, samples are
-- SI. The Coulomb kernel is genuinely three dimensional, not logarithmic.
-- Analytical exclusion disks belong to views and never alter the kernel.
-- Planar finite elements delegate geometry, assembly and raw gradients to
-- luafemm-electric-fe; they never reinterpret a Coulomb point charge.
-- @module luafemm-electrostatic
local E = {}
local C = require("luafemm-physical")
-- Numerical vacuum permittivity in SI, sufficient for teaching diagrams.
E.epsilon0 = 8.8541878128e-12
local function check(v, s)
    assert(v, "luafemm electrostatic: " .. s)
end
function E.number(v, name)
    v = tonumber(v)
    check(v and v == v and math.abs(v) < math.huge, "finite number required for " .. name)
    return v
end
function E.positive(v, name)
    v = E.number(v, name)
    check(v > 0, name .. " must be positive")
    return v
end
function E.integer(v, name, lo, hi)
    v = E.number(v, name)
    check(v == math.floor(v) and v >= lo and v <= hi, name .. " outside integer range")
    return v
end
local scales = { m = 1, cm = 0.01, mm = 0.001, um = 1e-6 }
--- Create a mutable source builder. Freeze before rendering or profiling.
function E.new(o)
    o = o or {}
    local unit = o.length_unit or "cm"
    check(scales[unit], "unknown length unit")
    local D = require("luafemm-dielectrics")
    local finite_elements = o.model == "finite elements" or o.model == "source correction"
    local dielectric = finite_elements
            and o.material
            and o.material ~= "custom"
            and o.material ~= ""
            and D.get(o.material)
        or nil
    local tensor = D.tensor(o)
    if not finite_elements then
        dielectric = D.scalar(o.material)
        check(
            tensor.xx == tensor.yy and math.abs(tensor.xy) < 1e-14,
            "analytical models require isotropic permittivity"
        )
    end
    local er =
        E.positive(o.epsilon_r or (dielectric and dielectric.ex) or 1, "relative permittivity")
    if not finite_elements then
        er = tensor.xx
    end
    local m = {
        sources = {},
        material = dielectric and dielectric.id or nil,
        unit = scales[unit],
        epsilon = E.epsilon0 * er,
        tensor = tensor,
        bounds = {
            E.number(o.xmin or (o.problem_type == "axisymmetric" and 0 or -3), "xmin"),
            E.number(o.ymin or -2, "ymin"),
            E.number(o.xmax or 3, "xmax"),
            E.number(o.ymax or 2, "ymax"),
        },
    }
    check(m.bounds[1] < m.bounds[3] and m.bounds[2] < m.bounds[4], "empty view window")
    m.descriptor = {
        schema = 1,
        physics = "electrostatics",
        backend = "analytic",
        geometry = "coplanar point charges",
        physical_dimension = 3,
        evaluation_dimension = 2,
        length_scale = m.unit,
        length_unit = unit,
        potential = "electric scalar potential",
        reference = "zero at infinity",
        medium = "homogeneous linear isotropic",
        epsilon = m.epsilon,
        validity = "all finite points except nonzero sources",
        units = { length = "m", charge = "C", V = "V", E = "V/m", D = "C/m^2" },
    }
    local kind = o.model or "point charges"
    check(
        o.problem_type == nil or o.problem_type == "planar" or finite_elements,
        "problem type requires a finite-element model"
    )
    if finite_elements then
        require("luafemm-electric-fe").configure(m, o)
        if kind == "source correction" then
            require("luafemm-electric-hybrid").configure(m, o)
        end
    elseif kind ~= "point charges" then
        require("luafemm-electric-classics").configure(m, o)
    else
        for _, key in ipairs({ "volume_density", "surface_density", "voltage", "reference_radius" }) do
            check(o[key] == nil or o[key] == "", key .. " requires a symmetry model")
        end
    end
    return m
end
--- Add one charge. Coincident sources are summed before evaluation.
-- A zero net source is removable, including its singularity and view mask.
function E.add(m, x, y, q, r)
    check(not m.classic and not m.fe, "point charges cannot be added to this model")
    check(not m.frozen, "declare all charges before rendering or capturing profiles")
    x, y, q = E.number(x, "source x"), E.number(y, "source y"), E.number(q, "charge")
    r = E.positive(r or 0.22, "exclusion radius")
    m.sources[#m.sources + 1] = { x = x, y = y, q = q, radius = r }
end
function E.freeze(m)
    if m.fe then
        return require("luafemm-electric-fe").solve(m)
    end
    if m.frozen then
        return m
    end
    local merged = {}
    for _, s in ipairs(m.sources) do
        local target
        for _, t in ipairs(merged) do
            if s.x == t.x and s.y == t.y then
                target = t
                break
            end
        end
        if target then
            target.q = target.q + s.q
            target.radius = math.max(target.radius, s.radius)
        else
            merged[#merged + 1] = C.copy(s)
        end
    end
    m.sources = {}
    for _, s in ipairs(merged) do
        if s.q ~= 0 then
            m.sources[#m.sources + 1] = s
        end
    end
    m.frozen = true
    return m
end
--- Return an independent frozen snapshot, with no reference to a mutable builder.
function E.snapshot(m)
    return E.freeze(C.copy(m))
end
function E.sample(m, x, y, side)
    x, y = E.number(x, "x"), E.number(y, "y")
    E.freeze(m)
    if m.fe then
        return require("luafemm-electric-fe").sample(m, x, y, side)
    end
    if m.classic then
        return require("luafemm-electric-classics").sample(m, x, y, side)
    end
    local v, ex, ey = 0, 0, 0
    local k = 1 / (4 * math.pi * m.epsilon)
    for _, s in ipairs(m.sources) do
        local dx, dy = (x - s.x) * m.unit, (y - s.y) * m.unit
        local rr = dx * dx + dy * dy
        if rr == 0 then
            return nil, "point-charge singularity"
        end
        local r = math.sqrt(rr)
        local a = k * s.q / r
        v = v + a
        ex = ex + a * dx / rr
        ey = ey + a * dy / rr
    end
    local mag = math.sqrt(ex * ex + ey * ey)
    if not (math.abs(v) < math.huge and mag < math.huge) then
        return nil, "Coulomb evaluation overflow"
    end
    return {
        V = v,
        Ex = ex,
        Ey = ey,
        E = mag,
        Dx = m.epsilon * ex,
        Dy = m.epsilon * ey,
        D = m.epsilon * mag,
    }
end
--- Minimum clearance from all source-exclusion disks, in model units.
function E.clearance(m, x, y)
    local d = math.huge
    for _, s in ipairs(m.sources) do
        d = math.min(d, math.sqrt((x - s.x) ^ 2 + (y - s.y) ^ 2) - s.radius)
    end
    return d
end
--- Return SI integral properties with explicit per-length/per-area names.
function E.properties(m, electrode)
    E.freeze(m)
    if m.fe then
        return require("luafemm-electric-fe").properties(m, electrode)
    end
    check(not electrode or electrode == "", "electrode quantities require finite elements")
    if m.classic then
        return require("luafemm-electric-classics").properties(m)
    end
    local charge = 0
    for _, s in ipairs(m.sources) do
        charge = charge + s.q
    end
    return { charge = charge }
end
return E
