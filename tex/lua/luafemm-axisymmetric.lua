-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Axisymmetric A_phi finite elements in the physical (r,z) half-plane.
-- A_phi is continuous P1; Br=-dA/dz and Bz=dA/dr+A/r. Volume integrals
-- use 2*pi*r dr dz. Assembly omits the common 2*pi factor, not the radius.
-- @module luafemm-axisymmetric
local X = {}
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm axisymmetric: " .. s)
end
--- Check combinations before any meshing or solving.
-- @tparam table m Model with validated scalar options.
function X.validate(m)
    check(m.options.xmin >= 0, "radius must be nonnegative (set xmin=0)")
    check(m.options.field_model == "open", "ideal confined mode is planar only")
    for _, p in ipairs(m.pairs or {}) do
        check(p.kelvin, "axisymmetric pairing requires a Kelvin exterior interface")
    end
    for _, p in pairs(m.frontiers or {}) do
        check(
            not p.condition or not p.condition.open_radius,
            "circular open boundaries are planar only"
        )
    end
end
-- Positive symmetric degree-five triangle rule; weights sum to one.
local rule = { { 1 / 3, 1 / 3, 1 / 3, 0.225 } }
for _, p in ipairs({
    { 0.059715871789770, 0.470142064105115, 0.132394152788506 },
    { 0.797426985353087, 0.101286507323456, 0.125939180544827 },
}) do
    for k = 1, 3 do
        local q = { p[2], p[2], p[2], p[3] }
        q[k] = p[1]
        rule[#rule + 1] = q
    end
end
--- Prepare SI quadrature and the curl of each local azimuthal basis function.
-- @tparam table m Model; nodes are supplied separately during reconstruction.
-- @tparam table e Prepared triangle with ids, gradients and material.
-- @tparam table nodes SI nodal coordinates.
function X.prepare(m, e, nodes)
    e.quadrature = {}
    e.meridional_current = 0
    for _, p in ipairs(rule) do
        local r, z = 0, 0
        for i, id in ipairs(e.ids) do
            r = r + p[i] * nodes[id][1]
            z = z + p[i] * nodes[id][2]
        end
        check(r > 0, "element has no positive radius")
        local q =
            { r = r, z = z, n = { p[1], p[2], p[3] }, br = {}, bz = {}, weight = e.area * p[4] * r }
        for i = 1, 3 do
            q.br[i] = -e.gy[i]
            q.bz[i] = e.gx[i] + p[i] / r
        end
        q.current = m.options.source and m.options.source(r, z) or e.current
        check(
            type(q.current) == "number"
                and q.current == q.current
                and math.abs(q.current) < math.huge,
            "invalid source"
        )
        e.meridional_current = e.meridional_current + e.area * p[4] * q.current
        e.quadrature[#e.quadrature + 1] = q
    end
end
--- Evaluate physical components at an SI point in a known triangle.
-- The axis uses the regular one-sided limit A_phi=O(r).
-- @tparam table m Prepared model owning the field.
-- @tparam table e Containing triangle.
-- @tparam number r Nonnegative SI radius.
-- @tparam number z SI axial coordinate.
-- @tparam[opt] table A Nodal potential; defaults to the solved m.A.
-- @return Br, Bz in T, A_phi in T m, Hr, Hz in A/m.
function X.field(m, e, r, z, A)
    A = A or m.A
    local p = m.nodes[e.ids[1]]
    local b = e.gx[2] * (r - p[1]) + e.gy[2] * (z - p[2])
    local c = e.gx[3] * (r - p[1]) + e.gy[3] * (z - p[2])
    local n = { 1 - b - c, b, c }
    local a, ar, az = 0, 0, 0
    for i, id in ipairs(e.ids) do
        a = a + n[i] * A[id]
        ar = ar + e.gx[i] * A[id]
        az = az + e.gy[i] * A[id]
    end
    local br, bz = -az, ar + (r > 0 and a / r or ar)
    if r == 0 then
        br = 0
        a = 0
    end
    local hr, hz = require("luafemm-static").element_law(m, e, br, bz, r, z)
    return br,
        bz,
        a,
        hr / mu0 - require("luafemm-static").coercivity(m, e) * e.mx,
        hz / mu0 - require("luafemm-static").coercivity(m, e) * e.my
end
--- Add nonlinear volume residuals and symmetric material tangents.
-- @tparam table m Prepared model.
-- @tparam table A Nodal A_phi values.
-- @tparam table R Residual, updated in place.
-- @tparam[opt] table K Sparse tangent, updated in place.
function X.assemble(m, A, R, K)
    local law = require("luafemm-static").element_law
    for _, e in ipairs(m.triangles) do
        for _, q in ipairs(e.quadrature) do
            local br, bz = 0, 0
            for i, id in ipairs(e.ids) do
                br = br + q.br[i] * A[id]
                bz = bz + q.bz[i] * A[id]
            end
            local hr, hz, rr, rz, zz = law(m, e, br, bz, q.r, q.z)
            hr, hz =
                hr - mu0 * require("luafemm-static").coercivity(m, e) * e.mx,
                hz - mu0 * require("luafemm-static").coercivity(m, e) * e.my
            for i, row in ipairs(e.dofs) do
                if row > 0 then
                    R[row] = R[row]
                        + q.weight
                            * (m.nodes[e.ids[i]].dof_sign or 1)
                            * (hr * q.br[i] + hz * q.bz[i] - mu0 * q.current * q.n[i])
                    if K then
                        for j, col in ipairs(e.dofs) do
                            if col > 0 then
                                K[row][col] = (K[row][col] or 0)
                                    + q.weight
                                        * (m.nodes[e.ids[i]].dof_sign or 1)
                                        * (m.nodes[e.ids[j]].dof_sign or 1)
                                        * (rr * q.br[i] * q.br[j] + rz * (q.br[i] * q.bz[j] + q.bz[i] * q.br[j]) + zz * q.bz[i] * q.bz[j])
                            end
                        end
                    end
                end
            end
        end
    end
end
--- Add exact radially weighted affine edge integrals.
-- Boundary convention: -H.t = c*A_phi+g on the meridional CCW tangent.
-- @tparam table m Model with prepared weak edges and nodal degrees.
-- @tparam table A Current nodal azimuthal potential.
-- @tparam table R Residual updated in place.
-- @tparam[opt] table K Tangent updated in place.
function X.boundary(m, A, R, K)
    for _, e in ipairs(m.boundary_edges) do
        check(not e.open_radius, "circular open boundaries are planar only")
        local ra, rb = m.nodes[e.ids[1]][1], m.nodes[e.ids[2]][1]
        local mass =
            { { (3 * ra + rb) / 12, (ra + rb) / 12 }, { (ra + rb) / 12, (ra + 3 * rb) / 12 } }
        for i, id in ipairs(e.ids) do
            local row = m.nodes[id].dof
            if row then
                for j, jd in ipairs(e.ids) do
                    local w = mu0 * e.length * mass[i][j] * (m.nodes[id].dof_sign or 1)
                    R[row] = R[row] + w * (e.coefficient * A[jd] + e.g[j])
                    local col = m.nodes[jd].dof
                    if K and col then
                        K[row][col] = (K[row][col] or 0)
                            + w * e.coefficient * (m.nodes[jd].dof_sign or 1)
                    end
                end
            end
        end
    end
end
--- Radial first moment of a face in model units cubed, including holes.
-- @tparam table face Arrangement face; first contour is outer, others are holes.
-- @treturn number Integral of radius over meridional area, in model units cubed.
function X.moment(face)
    local total = 0
    for k, p in ipairs(face.contours) do
        local area, moment = 0, 0
        for i, a in ipairs(p) do
            local b = p[i % #p + 1]
            local c = a[1] * b[2] - b[1] * a[2]
            area = area + c
            moment = moment + (a[1] + b[1]) * c / 6
        end
        total = total + (k == 1 and 1 or -1) * moment * (area < 0 and -1 or 1)
    end
    return total
end
--- Exact integral of r*A_phi over an element (two affine factors).
-- @tparam table m Model with nodal A_phi in T m.
-- @tparam table e Prepared triangle.
-- @treturn number Integral of full flux over meridional area, in Wb m^2.
function X.linkage(m, e)
    local total = 0
    for i, id in ipairs(e.ids) do
        for j, jd in ipairs(e.ids) do
            total = total + m.nodes[id][1] * m.A[jd] * (i == j and 2 or 1) * e.area / 12
        end
    end
    return 2 * math.pi * total
end
--- Extract flux contours on a uniformly subdivided drawing mesh.
-- Subdivision affects rendering only, never the solve. Shared edge midpoint
-- identifiers preserve stitching; the quadratic flux is sampled exactly.
-- @tparam table m Solved model; left unchanged by drawing subdivision.
-- @tparam integer count Requested automatic line count.
-- @tparam[opt] table levels Explicit full-flux levels, in Wb.
-- @treturn table Oriented paths in model units.
function X.contours(m, count, levels)
    local nodes, values, tri = {}, {}, {}
    for i, p in ipairs(m.nodes) do
        nodes[i] = { p[1], p[2] }
        values[i] = m.A[i]
    end
    for _, e in ipairs(m.triangles) do
        tri[#tri + 1] = { ids = e.ids, source = e }
    end
    for _ = 1, m.options.contour_subdivisions do
        local mids, nexttri = {}, {}
        local function midpoint(a, b)
            local key = math.min(a, b) .. ":" .. math.max(a, b)
            if not mids[key] then
                local id = #nodes + 1
                mids[key] = id
                nodes[id] = { (nodes[a][1] + nodes[b][1]) / 2, (nodes[a][2] + nodes[b][2]) / 2 }
                values[id] = (values[a] + values[b]) / 2
            end
            return mids[key]
        end
        for _, t in ipairs(tri) do
            local a, b, c = table.unpack(t.ids)
            local ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            for _, ids in ipairs({ { a, ab, ca }, { ab, b, bc }, { ca, bc, c }, { ab, bc, ca } }) do
                nexttri[#nexttri + 1] = { ids = ids, source = t.source }
            end
        end
        tri = nexttri
    end
    for i, a in ipairs(values) do
        values[i] = 2 * math.pi * nodes[i][1] * a
    end
    for _, t in ipairs(tri) do
        local r, z = 0, 0
        for _, id in ipairs(t.ids) do
            r = r + nodes[id][1] / 3
            z = z + nodes[id][2] / 3
        end
        t.bx, t.by = X.field(m, t.source, r, z)
        t.source = nil
    end
    return require("luafemm").contours({
        stats = {},
        nodes = nodes,
        A = values,
        triangles = tri,
        options = { unit = m.options.unit, field_model = "open" },
    }, count, levels)
end
--- Accumulate physical full-revolution regional integrals using assembly quadrature.
-- Internal helper for the named postprocessor; radial resultant and torque vanish.
-- @tparam table m Solved model.
-- @tparam table e Selected triangle with prepared quadrature.
-- @tparam table r Integral accumulator, updated in place.
-- @tparam string quantity integrals or force.
-- @tparam string method Selected force method (Lorentz is local).
function X.regional(m, e, r, quantity, method)
    local post = require("luafemm-post")
    r.area = r.area + e.area
    r.current = r.current + (m.transient and e.time_current or e.meridional_current)
    for _, q in ipairs(e.quadrature) do
        local br, bz, _, hr, hz = X.field(m, e, q.r, q.z)
        local b = math.sqrt(br * br + bz * bz)
        local volume = 2 * math.pi * q.weight
        local fields = {
            Bx = br,
            By = bz,
            B = b,
            Hx = hr,
            Hy = hz,
            H = math.sqrt(hr * hr + hz * hz),
            x = q.r,
            y = q.z,
        }
        for key, value in pairs(fields) do
            r[key] = r[key] + volume * value
        end
        r.volume = r.volume + volume
        local energy
        if e.p.nonlinear_sheet and not m.operating then
            energy = require("luafemm-lamination").energy(e.p, br, bz)
        elseif e.p.bh and not m.operating then
            energy = post.primitive(e.p, b)
        else
            local h1, h2 = require("luafemm-static").element_law(m, e, br, bz, q.r, q.z)
            energy = (br * h1 + bz * h2) / (2 * mu0)
        end
        energy = energy - require("luafemm-static").coercivity(m, e) * (e.mx * br + e.my * bz)
        r.energy = r.energy + volume * energy
        r.coenergy = r.coenergy + volume * (br * hr + bz * hz - energy)
        if quantity == "force" and method == "lorentz" then
            check(
                not e.p.bh and e.p.mur_x == 1 and e.p.mur_y == 1 and e.p.hc == 0,
                "Lorentz body must be nonmagnetic"
            )
            local current = m.transient and require("luafemm-time").field(m, e, q.r, q.z).J
                or q.current
            r.Fy = r.Fy - volume * current * br
        end
    end
end
return X
