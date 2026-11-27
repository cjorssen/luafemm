-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Analytic source kernels with a regular scalar P1 correction.
-- Geometry and SI conventions match the electric FE provider. Only correction
-- coefficients are discretized; total gradients are evaluated at each point.
-- The reference medium must occupy a positive ball around every source.
-- Affine total-potential constraints and analytic reference fluxes preserve
-- conductor charges and periodic traces. Divergent total self-energy and native
-- answer conversion remain unavailable. Snapshots contain passive data only.
-- @module luafemm-electric-hybrid
local H = {}
local E = require("luafemm-electrostatic")
local F = require("luafemm-electric-fe")
local P = require("luafemm-physical")
local metric = require("luafemm-electric-metric")
local pi, abs, sqrt = math.pi, math.abs, math.sqrt
local function check(ok, message)
    assert(ok, "luafemm electric source correction: " .. message)
end

local function kernel(model, sources, er, reference, x, y)
    local v, ex, ey = 0, 0, 0
    for _, s in ipairs(sources) do
        if s.q ~= 0 then
            local dx, dy = (x - s.x) * model.unit, (y - s.y) * model.unit
            local r2 = dx * dx + dy * dy
            if r2 == 0 then
                return nil, "analytical source singularity"
            end
            local r = sqrt(r2)
            if model.fe.axisymmetric then
                local a = s.q / (4 * pi * E.epsilon0 * er)
                v, ex, ey = v + a / r, ex + a * dx / (r2 * r), ey + a * dy / (r2 * r)
            else
                local a = s.q / (2 * pi * E.epsilon0 * er)
                v = v - a * math.log(r / (reference * model.unit))
                ex, ey = ex + a * dx / r2, ey + a * dy / r2
            end
        end
    end
    check(
        v == v and math.abs(v) < math.huge and math.abs(ex) < math.huge and math.abs(ey) < math.huge,
        "source field overflow"
    )
    return { V = v, Ex = ex, Ey = ey }
end

local function segment_distance(s, a, b)
    local dx, dy = b[1] - a[1], b[2] - a[2]
    local t =
        math.max(0, math.min(1, ((s.x - a[1]) * dx + (s.y - a[2]) * dy) / (dx * dx + dy * dy)))
    return sqrt((s.x - a[1] - t * dx) ^ 2 + (s.y - a[2] - t * dy) ^ 2)
end

local function value(v, fallback)
    if v == nil or v == "" then
        return fallback
    end
    return v
end
function H.configure(m, o)
    check(
        m.tensor.xx == m.tensor.yy and m.tensor.xy == 0,
        "source correction requires isotropic dielectrics"
    )
    check(
        not m.fe.axisymmetric or value(o.reference_radius, nil) == nil,
        "reference radius is only available for planar source correction"
    )
    m.hybrid = {
        quadrature_tolerance = E.positive(
            value(o.source_quadrature_tolerance, 1e-6),
            "source quadrature tolerance"
        ),
        quadrature_depth = E.integer(
            value(o.source_quadrature_depth, 8),
            "maximum source quadrature depth",
            1,
            12
        ),
        quadrature_cells = E.integer(
            value(o.source_quadrature_cells, 250000),
            "maximum source quadrature cells",
            1,
            2000000
        ),
        epsilon_r = m.tensor.xx,
        reference_radius = E.positive(value(o.reference_radius, 1), "reference radius"),
    }
    check(m.hybrid.quadrature_tolerance < 1, "source quadrature tolerance must be less than one")
    local d = m.descriptor
    d.backend = "source correction"
    d.epsilon = E.epsilon0 * m.hybrid.epsilon_r
    d.validity =
        "finite planar or axisymmetric domain; singular analytic sources and raw P1 correction"
    d.units.line_charge_density = "C/m"
end
--- Add a physical source in model coordinates, before any mesh or evaluation.
-- kind is explicit: a line strength is never silently read as a 3D charge.
function H.add(m, x, y, o)
    check(m.hybrid, "electric source nodes require model=source correction")
    check(
        not m.fe.nodes and not m.frozen,
        "declare all electric sources before meshing or evaluation"
    )
    o = o or {}
    local axial = m.fe.axisymmetric
    check(
        o.kind == (axial and "axis point charge" or "line charge"),
        "source kind does not match problem type"
    )
    local strength, forbidden
    if axial then
        strength, forbidden = o.charge, o.charge_per_length
    else
        strength, forbidden = o.charge_per_length, o.charge
    end
    check(value(forbidden, nil) == nil, "incompatible source strength unit")
    local q = E.number(value(strength, nil), axial and "charge" or "charge per length")
    x, y = E.number(x, "source x"), E.number(y, "source y")
    check(not axial or x == 0, "axis point charge must lie on r=0")
    m.sources[#m.sources + 1] = {
        x = x,
        y = y,
        q = q,
        kind = o.kind,
        guard = E.positive(value(o.guard_radius, 0.1), "source guard radius"),
        radius = E.positive(value(o.radius, 0.22), "exclusion radius"),
    }
end
--- Validate before meshing. Coincident source strengths combine exactly; a
-- zero total removes both the singularity and its numerical/display guards.
function H.validate(m)
    local f = m.fe
    check(not f.native, "source correction does not accept an imported native mesh")
    check(not f.axisymmetric or f.bounds[1] == 0, "axial source domain must include r=0")
    check(f.density == 0, "additional volume charge is unavailable with source correction")
    for _, r in ipairs(f.regions) do
        check(r.density == 0, "additional volume charge is unavailable with source correction")
    end
    for _, p in ipairs(f.paths) do
        check(
            p.kind == "electrode" or p.kind == "boundary",
            "additional surface charge is unavailable with source correction"
        )
    end
    local merged = {}
    for _, s in ipairs(m.sources) do
        local found
        for _, t in ipairs(merged) do
            if s.x == t.x and s.y == t.y then
                t.q, t.guard, t.radius =
                    t.q + s.q, math.max(t.guard, s.guard), math.max(t.radius, s.radius)
                found = true
                break
            end
        end
        if not found then
            merged[#merged + 1] = P.copy(s)
        end
    end
    m.sources = {}
    for _, s in ipairs(merged) do
        E.number(s.q, "combined source strength")
        if s.q ~= 0 then
            local b = f.bounds
            check(
                s.y - s.guard > b[2]
                    and s.y + s.guard < b[4]
                    and s.x + s.guard < b[3]
                    and (f.axisymmetric or s.x - s.guard > b[1]),
                "source guard must be inside the physical domain"
            )
            m.sources[#m.sources + 1] = s
        end
    end
end
--- Certify homogeneous source balls against the actual constrained mesh and
-- every thin electrode. A rendering exclusion radius never changes this test.
function H.check_mesh(m)
    local f, er = m.fe, m.hybrid.epsilon_r
    for _, s in ipairs(m.sources) do
        for _, p in ipairs(f.paths) do
            if p.kind == "electrode" then
                for j = 1, p.closed and #p.points or #p.points - 1 do
                    local a, b = p.points[j], p.points[j % #p.points + 1]
                    if not (f.axisymmetric and a[1] == 0 and b[1] == 0) then
                        check(
                            segment_distance(s, a, b) > s.guard,
                            "electrode intersects a source guard"
                        )
                    end
                end
            end
        end
        for _, t in ipairs(f.triangles) do
            if t.electrode or t.tensor.xx ~= er or t.tensor.yy ~= er or t.tensor.xy ~= 0 then
                local w = F.barycentric(m, t, s.x, s.y)
                check(
                    math.min(w[1], w[2], w[3]) < -1e-12,
                    "source lies in a conductor or a different dielectric"
                )
                for j = 1, 3 do
                    check(
                        segment_distance(s, f.nodes[t.ids[j]], f.nodes[t.ids[j % 3 + 1]]) > s.guard,
                        "material interface intersects a source guard"
                    )
                end
            end
        end
    end
end
function H.base(m, x, y)
    return kernel(m, m.sources, m.hybrid.epsilon_r, m.hybrid.reference_radius, x, y)
end
--- Assemble only the regular correction. The sparse stiffness is shared with
-- the ordinary scalar formulation, but total and correction traces stay apart.
function H.solve(m)
    local h = m.hybrid
    if h.solution then
        return m
    end
    local f = F.mesh(m)
    local er = h.epsilon_r
    local function base(x, y)
        return H.base(m, x, y)
    end
    local integrate, cell_count = require("luafemm-electric-source-quadrature").integrator(m, base)
    local K, b = {}, {}
    for i = 1, #f.nodes do
        K[i], b[i] = {}, 0
    end
    local contrast_elements = 0
    for _, t in ipairs(f.triangles) do
        if not t.electrode then
            local contrast = t.tensor.xx ~= er or t.tensor.yy ~= er or t.tensor.xy ~= 0
            if contrast then
                contrast_elements = contrast_elements + 1
            end
            local weight = metric.triangle(m, t)
            for i, u in ipairs(t.ids) do
                for j, v in ipairs(t.ids) do
                    K[u][v] = (K[u][v] or 0)
                        + t.area
                            * weight
                            * (t.gx[i] * (t.tensor.xx * t.gx[j] + t.tensor.xy * t.gy[j]) + t.gy[i] * (t.tensor.xy * t.gx[j] + t.tensor.yy * t.gy[j]))
                end
            end
            if contrast then
                local ix, iy = integrate(t)
                for j, id in ipairs(t.ids) do
                    b[id] = b[id] + ix * t.gx[j] + iy * t.gy[j]
                end
            end
        end
    end

    local boundary = require("luafemm-electric-hybrid-boundary")
    local boundary_cells = boundary.assemble(m, K, b, base, h.quadrature_cells - cell_count())
    local fixed, offset = {}, {}
    for id in pairs(f.owners) do
        local p = f.nodes[id]
        offset[id] = -assert(base(p[1], p[2])).V
    end
    for id in pairs(f.periodic_nodes) do
        local p = f.nodes[id]
        offset[id] = -assert(base(p[1], p[2])).V
    end
    for id, potential in pairs(f.fixed) do
        fixed[id] = potential + (offset[id] or 0)
    end
    local reduced = {}
    for k, v in pairs(f) do
        reduced[k] = v
    end
    reduced.fixed = fixed
    local u, iterations, stats
    if h.imported_solution then
        -- Extended answers own correction coefficients, never total nodal V.
        -- Reconstruct gradients/reactions without invoking a linear solver.
        u, iterations, stats =
            require("luafemm-physical").copy(h.imported_solution),
            0,
            { relative_residual = h.imported_residual or 0 }
        for id, v in pairs(fixed) do
            check(
                abs(u[id] - v) <= 1e-8 * math.max(1, abs(v)),
                "restored correction contradicts a prescribed potential"
            )
        end
    else
        u, iterations, stats =
            require("luafemm-electric-constraints").solve(reduced, K, b, E.epsilon0, offset)
    end
    local flux_cells =
        boundary.reactions(m, K, b, u, base, h.quadrature_cells - cell_count() - boundary_cells)
    for _, t in ipairs(f.triangles) do
        local ex, ey = 0, 0
        if not t.electrode then
            for j, id in ipairs(t.ids) do
                ex = ex - t.gx[j] * u[id] / m.unit
                ey = ey - t.gy[j] * u[id] / m.unit
            end
        end
        check(abs(ex) < math.huge and abs(ey) < math.huge, "correction field overflow")
        t.correction_ex, t.correction_ey = ex, ey
    end
    h.solution, h.iterations, h.residual = u, iterations, stats.relative_residual
    h.contrast_elements = contrast_elements
    h.quadrature_evaluations = cell_count() + boundary_cells + flux_cells
    m.frozen = true
    return m
end
--- Evaluate the exact kernel plus the raw correction within one chosen cell.
-- An explicit cell is used by one-sided profiles and local streamline stages.
function H.sample(m, x, y, element)
    H.solve(m)
    local i, w
    if type(element) == "number" then
        i = element
        w = F.barycentric(m, m.fe.triangles[i], x, y)
    else
        i, w = F.locate(m, x, y)
    end
    if not i then
        return nil, "outside source-correction domain"
    end
    local t = m.fe.triangles[i]
    local v, ex, ey
    if t.electrode then
        v, ex, ey = m.fe.voltages[t.electrode], 0, 0
    else
        local s, why = H.base(m, x, y)
        if not s then
            return nil, why
        end
        v, ex, ey = s.V, s.Ex + t.correction_ex, s.Ey + t.correction_ey
        for j, id in ipairs(t.ids) do
            v = v + w[j] * m.hybrid.solution[id]
        end
    end
    local e = sqrt(ex * ex + ey * ey)
    local dx = E.epsilon0 * (t.tensor.xx * ex + t.tensor.xy * ey)
    local dy = E.epsilon0 * (t.tensor.xy * ex + t.tensor.yy * ey)
    local d = sqrt(dx * dx + dy * dy)
    check(abs(v) < math.huge and e < math.huge and d < math.huge, "total field overflow")
    return {
        V = v,
        Ex = ex,
        Ey = ey,
        E = e,
        Dx = dx,
        Dy = dy,
        D = d,
        Er = m.fe.axisymmetric and ex or nil,
        Ez = m.fe.axisymmetric and ey or nil,
        Dr = m.fe.axisymmetric and dx or nil,
        Dz = m.fe.axisymmetric and dy or nil,
        element = i,
        region = t.region,
        electrode = t.electrode,
    }
end
function H.properties(m, name)
    H.solve(m)
    local f = m.fe
    if name and name ~= "" then
        check(f.charges[name] ~= nil, "unknown electrode: " .. name)
        return {
            charge = f.charges[name] * f.depth,
            ["charge per length"] = not f.axisymmetric and f.charges[name] or nil,
            voltage = f.voltages[name],
        }
    end
    return {
        nodes = #f.nodes,
        triangles = #f.triangles,
        iterations = m.hybrid.iterations,
        ["relative residual"] = m.hybrid.residual,
        charge = f.source_charge * f.depth,
        ["charge per length"] = not f.axisymmetric and f.source_charge or nil,
        ["charge balance"] = f.balance * f.depth,
        ["charge balance per length"] = not f.axisymmetric and f.balance or nil,
        ["periodic flux"] = f.periodic_flux * f.depth,
    }
end
return H
