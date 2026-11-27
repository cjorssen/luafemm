-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Axial H diffusion in cross-sections of infinitely long isolated conductors.
-- The surrounding connected insulator has uniform imposed H_z. The stored
-- unknown is u=H_z/H_applied-1, dimensionless; it is NOT a vector potential.
-- J=(dH_z/dy,-dH_z/dx,0). Displacement current and end effects are neglected.
-- @module luafemm-axial
local X = {}
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm axial: " .. s)
end
--- Validate the reduced long-solenoid model, never silently reinterpret sources.
-- @tparam table m Model before meshing or cache restoration.
function X.validate(m)
    local o = m.options
    check(o.frequency > 0, "axial diffusion requires frequency>0")
    check(o.mesher == "delaunay", "axial diffusion requires the constrained delaunay mesher")
    check(not m.import_document, "axial diffusion does not import FEM problems")
    check(not next(m.circuits or {}), "use applied field, not a meridional winding circuit")
    check(
        not next(m.frontiers or {}),
        "axial diffusion uses imposed H on conductor surfaces, not named frontiers"
    )
    for _, b in pairs(m.boundaries) do
        check(
            b.type == "dirichlet" and b.value == 0 and b.dx == 0 and b.dy == 0,
            "use applied field instead of vector-potential boundary conditions"
        )
    end
    for _, p in pairs(m.materials) do
        check(p.mur_x == p.mur_y, "axial diffusion requires isotropic permeability")
        check(p.conductivity > 0 or p.mur_x == 1, "insulators must be nonmagnetic")
    end
    for _, r in ipairs(m.regions) do
        check(
            r.current == 0 and (r.current_phase or 0) == 0 and r.conductor ~= "stranded",
            "axial diffusion accepts passive solid conductors only"
        )
    end
    for _, f in ipairs(m.topology and m.topology.faces or {}) do
        check(
            not f.hole,
            "excluded holes are unsupported; draw connected insulating regions instead"
        )
    end
end
--- Compute fixed nodes using triangle conductivity and exterior connectivity.
-- Enclosed insulating cavities need a coupled flux constraint and are rejected.
-- @tparam table nodes SI mesh nodes.
-- @tparam table triangles Classified triangles with ids and material p.
-- @tparam table exterior Exterior node set.
-- @treturn table Zero-valued fixed node map.
function X.fixed(nodes, triangles, exterior)
    local fixed, edges, adj, queue, seen = {}, {}, {}, {}, {}
    for id in pairs(exterior) do
        fixed[id] = 0
    end
    for i, e in ipairs(triangles) do
        adj[i] = {}
        for j, a in ipairs(e.ids) do
            local b = e.ids[j % 3 + 1]
            local key = math.min(a, b) .. ":" .. math.max(a, b)
            local edge = edges[key] or { a = a, b = b }
            edges[key] = edge
            if edge.first then
                adj[i][#adj[i] + 1] = edge.first
                adj[edge.first][#adj[edge.first] + 1] = i
                edge.second = i
            else
                edge.first = i
            end
        end
        if e.p.conductivity == 0 then
            for _, id in ipairs(e.ids) do
                fixed[id] = 0
            end
        end
    end
    for _, edge in pairs(edges) do
        if not edge.second then
            fixed[edge.a], fixed[edge.b] = 0, 0
            local i = edge.first
            if triangles[i].p.conductivity == 0 and not seen[i] then
                seen[i] = true
                queue[#queue + 1] = i
            end
        end
    end
    local k = 1
    while k <= #queue do
        for _, j in ipairs(adj[queue[k]]) do
            if triangles[j].p.conductivity == 0 and not seen[j] then
                seen[j] = true
                queue[#queue + 1] = j
            end
        end
        k = k + 1
    end
    for i, e in ipairs(triangles) do
        check(
            e.p.conductivity > 0 or seen[i],
            "enclosed insulating cavity: its axial field cannot be prescribed independently"
        )
    end
    return fixed
end
--- Replace vector-potential degrees by the conductor-interior H perturbation.
-- @tparam table m Prepared model; scalar boundary values remain zero in caches.
function X.prepare(m)
    X.validate(m)
    local fixed = X.fixed(m.nodes, m.triangles, m.exterior_nodes)
    m.free, m.fixed_values = {}, fixed
    for id, p in ipairs(m.nodes) do
        p.dof = nil
        if fixed[id] == nil then
            m.free[#m.free + 1] = id
            p.dof = #m.free
        end
    end
    for _, e in ipairs(m.triangles) do
        for i, id in ipairs(e.ids) do
            e.dofs[i] = m.nodes[id].dof or 0
        end
    end
end
--- Assemble curl-resistivity diffusion in H, with a consistent complex mass.
-- @tparam table m Prepared axial model.
-- @treturn table Solved model storing dimensionless perturbation phasors.
function X.solve(m)
    local start = os.clock()
    local K, L, br, bi = {}, {}, {}, {}
    for i = 1, #m.free do
        K[i], L[i], br[i], bi[i] = {}, {}, 0, 0
    end
    local omega = 2 * math.pi * m.options.frequency
    for _, e in ipairs(m.triangles) do
        if e.p.conductivity > 0 then
            local eta = 1 / e.p.conductivity
            local c = omega * mu0 * e.p.mur_x
            for i, row in ipairs(e.dofs) do
                if row > 0 then
                    bi[row] = bi[row] - c * e.area / 3
                    for j, col in ipairs(e.dofs) do
                        if col > 0 then
                            K[row][col] = (K[row][col] or 0)
                                + eta * e.area * (e.gx[i] * e.gx[j] + e.gy[i] * e.gy[j])
                            L[row][col] = (L[row][col] or 0) + c * e.area * (i == j and 2 or 1) / 12
                        end
                    end
                end
            end
        end
    end
    local assembly_seconds = os.clock() - start
    local real, imag, it, res, telemetry = require("luafemm-linear").complex(m, K, L, br, bi)
    m.A_im = {}
    for id, p in ipairs(m.nodes) do
        m.A[id] = p.dof and real[p.dof] or 0
        m.A_im[id] = p.dof and imag[p.dof] or 0
    end
    m.stats = {
        newton = 0,
        cg = 0,
        gmres = it,
        linear = require("luafemm-linear").accumulate(nil, telemetry),
        assembly_seconds = assembly_seconds,
        residual = res,
        seconds = os.clock() - start,
        nodes = #m.nodes,
        triangles = #m.triangles,
    }
    require("luafemm").update_fields(m)
    require("luafemm-cache").store(m)
    return m
end
--- Return axial magnetic and transverse current phasors at an SI point.
-- @tparam table m Solved axial model.
-- @tparam table e Containing triangle.
-- @tparam number x SI horizontal coordinate.
-- @tparam number y SI vertical coordinate.
-- @treturn table Bz, Hz, Jx, Jy pairs; in-plane B/H and out-of-plane J vanish.
function X.field(m, e, x, y)
    local p = m.nodes[e.ids[1]]
    local b = e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2])
    local c = e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
    local basis = { 1 - b - c, b, c }
    local hr, hi, jxr, jxi, jyr, jyi = 1, 0, 0, 0, 0, 0
    for i, id in ipairs(e.ids) do
        hr, hi = hr + basis[i] * m.A[id], hi + basis[i] * m.A_im[id]
        if e.p.conductivity > 0 then
            jxr, jxi = jxr + e.gy[i] * m.A[id], jxi + e.gy[i] * m.A_im[id]
            jyr, jyi = jyr - e.gx[i] * m.A[id], jyi - e.gx[i] * m.A_im[id]
        end
    end
    local a = math.rad(m.options.applied_phase)
    local h = m.options.applied_field / mu0
    local function rotate(r, i)
        return { h * (r * math.cos(a) - i * math.sin(a)), h * (r * math.sin(a) + i * math.cos(a)) }
    end
    local hz = rotate(hr, hi)
    return {
        Bz = { mu0 * e.p.mur_x * hz[1], mu0 * e.p.mur_x * hz[2] },
        Hz = hz,
        Jx = rotate(jxr, jxi),
        Jy = rotate(jyr, jyi),
        Bx = { 0, 0 },
        By = { 0, 0 },
        Hx = { 0, 0 },
        Hy = { 0, 0 },
        J = { 0, 0 },
        Jeddy = { 0, 0 },
    }
end
--- Refresh axial diagnostics without interpreting the scalar unknown as A_z.
-- @tparam table m Solved model.
function X.update(m)
    local peak = 0
    for _, e in ipairs(m.triangles) do
        local x, y = 0, 0
        for _, id in ipairs(e.ids) do
            x = x + m.nodes[id][1] / 3
            y = y + m.nodes[id][2] / 3
        end
        local f = X.field(m, e, x, y)
        e.bx, e.by, e.hx, e.hy = 0, 0, 0, 0
        e.b = math.sqrt(f.Bz[1] ^ 2 + f.Bz[2] ^ 2)
        peak = math.max(peak, e.b)
    end
    m.stats.peak = peak
    m.solution_signature = require("luafemm-cache").signature(m)
end
--- Integrate physical volume, cycle-average energy and eddy Joule losses.
-- @tparam table m Solved model with positive depth.
-- @tparam table chosen Triangle selection map.
-- @treturn table SI observables; no fictitious longitudinal current is reported.
function X.integrals(m, chosen)
    local r = { area = 0, volume = 0, energy = 0, joule_loss = 0, eddy_loss = 0 }
    for i, e in ipairs(m.triangles) do
        if chosen[i] then
            r.area = r.area + e.area
            for _, q in ipairs(require("luafemm-harmonic").quadrature(m, e)) do
                local f = X.field(m, e, q.r, q.z)
                local v = q.weight * m.options.depth * m.options.unit
                r.volume = r.volume + v
                r.energy = r.energy + v * (f.Bz[1] * f.Hz[1] + f.Bz[2] * f.Hz[2]) / 4
                if e.p.conductivity > 0 then
                    r.joule_loss = r.joule_loss
                        + v
                            * (f.Jx[1] ^ 2 + f.Jx[2] ^ 2 + f.Jy[1] ^ 2 + f.Jy[2] ^ 2)
                            / (2 * e.p.conductivity)
                end
            end
        end
    end
    r.hysteresis_loss, r.lamination_loss, r.magnetic_loss, r.total_loss = 0, 0, 0, r.joule_loss
    r.eddy_loss, r.coenergy = r.joule_loss, r.energy
    return r
end
return X
