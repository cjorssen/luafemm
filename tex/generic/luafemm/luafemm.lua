-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Planar magnetostatics with P1 triangles and a nonlinear isotropic law.
-- This module is usable by texlua without a running TeX interpreter.
-- Public lengths use options.unit metres per model unit; assembly uses SI.
-- @module luafemm
local M = { version = "0.8.0-dev", mu0 = 4 * math.pi * 1e-7 }
local boundary_conditions = require("luafemm-boundary")
local geometry = require("luafemm-geometry")
local abs, sqrt, max = math.abs, math.sqrt, math.max
local function check(x, message)
    assert(x, "luafemm: " .. message)
end
local function finite(x)
    return type(x) == "number" and x == x and abs(x) < math.huge
end
local function dot(a, b)
    local s = 0
    for i = 1, #a do
        s = s + a[i] * b[i]
    end
    return s
end
local function norm(a)
    return sqrt(dot(a, a))
end
local function zeros(n)
    local a = {}
    for i = 1, n do
        a[i] = 0
    end
    return a
end

-- Models own their input tables. A caller may reuse a configuration, material
-- or polygon without subsequent preparation changing that caller's data.
local function copy(value)
    if type(value) ~= "table" then
        return value
    end
    local result = {}
    for key, item in pairs(value) do
        result[key] = copy(item)
    end
    return result
end

--- Construct an unsolved model, copying and validating the options.
-- @tparam[opt] table o Domain, mesh, solver and optional source/boundary callbacks.
-- @treturn table Model owning materials, regions and options; no mesh yet.
function M.new(o)
    o = copy(o or {})
    local m = { materials = {}, regions = {}, options = o }
    o.field_model = o.field_model or "open"
    check(o.field_model == "open" or o.field_model == "ideal", "unknown field model")
    o.unit = o.unit or 0.001
    o.h = o.h or 2
    check(o.depth == nil or (finite(o.depth) and o.depth > 0), "depth must be positive")
    o.xmin = o.xmin or -90
    o.xmax = o.xmax or 90
    o.ymin = o.ymin or -75
    o.ymax = o.ymax or 95
    o.tolerance = o.tolerance or 1e-7
    o.max_newton = o.max_newton or 60
    o.mesher = o.mesher or "grid"
    o.min_angle, o.max_area, o.max_steiner =
        o.min_angle or 0, o.max_area or 0, o.max_steiner or 5000
    for _, k in ipairs({ "unit", "h", "xmin", "xmax", "ymin", "ymax", "tolerance", "max_newton" }) do
        check(finite(o[k]), "non-finite problem option " .. k)
    end
    check(
        not o.mesh_tolerance or (finite(o.mesh_tolerance) and o.mesh_tolerance >= 0),
        "invalid mesh tolerance"
    )
    check(o.mesher == "grid" or o.mesher == "delaunay", "unknown mesher")
    check(
        finite(o.min_angle) and o.min_angle >= 0 and o.min_angle <= 30,
        "minimum angle must be between 0 and 30 degrees"
    )
    check(finite(o.max_area) and o.max_area >= 0, "invalid maximum area")
    check(
        finite(o.max_steiner) and o.max_steiner >= 0 and o.max_steiner % 1 == 0,
        "max Steiner points must be a nonnegative integer"
    )
    check(
        o.mesher == "delaunay" or (o.min_angle == 0 and o.max_area == 0),
        "quality refinement requires the delaunay mesher"
    )
    check(
        o.unit > 0 and o.h > 0 and o.xmax > o.xmin and o.ymax > o.ymin,
        "invalid domain or mesh size"
    )
    check(
        o.tolerance > 0 and o.max_newton >= 1 and o.max_newton % 1 == 0,
        "invalid solver tolerance or iteration limit"
    )
    boundary_conditions.configure(m)
    require("luafemm-cache").configure(m)
    M.material(m, "air", { mur = 1 })
    return m
end

--- Declare Dirichlet, Neumann or Robin data on one side of the rectangle.
-- All values use SI units; see luafemm-boundary.set for the option table.
-- @tparam table m Unmeshed model.
-- @tparam string side left, right, bottom, top or all.
-- @tparam[opt] table condition Boundary type and affine coefficients.
function M.boundary(m, side, condition)
    boundary_conditions.set(m, side, condition)
end

--- Define a material before meshing; explicit properties override a library base.
-- @tparam table m Model returned by new.
-- @tparam string name Model-local material name.
-- @tparam[opt] table p mur, bh (T/A per metre), coercivity, remanence, current, library.
function M.material(m, name, p)
    check(not m.nodes, "declare materials before meshing")
    p = copy(p or {})
    if type(p.library) == "string" then
        local base = require("luafemm-materials").properties(p.library)
        for k, v in pairs(p) do
            if k ~= "library" then
                base[k] = v
            end
        end
        p = base
    end
    p.mur = p.mur or 1
    p.current = p.current or 0
    check(finite(p.current), "material current density must be finite")
    check(finite(p.mur) and p.mur > 0, "relative permeability must be positive")
    p.remanence = p.remanence or 0
    p.coercivity = p.coercivity or 0
    check(
        finite(p.remanence) and p.remanence >= 0 and finite(p.coercivity) and p.coercivity >= 0,
        "remanence and coercivity must be finite and nonnegative"
    )
    check(
        p.remanence == 0 or (p.coercivity == 0 and not p.bh),
        "remanence requires a linear recoil law; use coercivity with a shifted B/H table"
    )
    p.hc = p.remanence > 0 and p.remanence / (M.mu0 * p.mur) or p.coercivity
    if p.bh then
        check(#p.bh >= 2 and p.bh[1][1] == 0 and p.bh[1][2] == 0, "B/H table must start at 0/0")
        for i = 2, #p.bh do
            check(
                finite(p.bh[i][1])
                    and finite(p.bh[i][2])
                    and p.bh[i][1] > p.bh[i - 1][1]
                    and p.bh[i][2] > p.bh[i - 1][2],
                "B and H must be finite and strictly increasing"
            )
        end
    end
    p.interpolation = p.interpolation or m.options.interpolation or "linear"
    check(p.interpolation == "linear" or p.interpolation == "femm", "unknown B-H interpolation")
    if p.bh and p.interpolation == "femm" then
        p.curve = require("luafemm-bh").prepare(p.bh)
    end
    m.materials[name] = p
end

--- Evaluate the scaled secant and differential material law.
-- Piecewise-linear H(B) uses the vacuum slope beyond the last point, keeping
-- magnetisation saturated. Coercivity is a separate source in the weak form.
-- @tparam table p Prepared material returned in model.materials.
-- @tparam number b Nonnegative flux-density magnitude, in tesla.
-- @treturn number mu0*H_table/B, using the initial slope at B=0.
-- @treturn number mu0*dH_table/dB.
function M.constitutive(p, b)
    if not p.bh then
        return 1 / p.mur, 1 / p.mur
    end
    if p.curve then
        local h, d = require("luafemm-bh").evaluate(p.curve, b)
        return b > 1e-15 and M.mu0 * h / b or M.mu0 * d, M.mu0 * d
    end
    local t = p.bh
    local lo, hi = t[1], t[2]
    if b > t[#t][1] then
        local h = t[#t][2] + (b - t[#t][1]) / M.mu0
        return M.mu0 * h / b, 1
    end
    for i = 2, #t do
        if b <= t[i][1] then
            lo, hi = t[i - 1], t[i]
            break
        end
    end
    local slope = (hi[2] - lo[2]) / (hi[1] - lo[1])
    local h = lo[2] + slope * (b - lo[1])
    return b > 1e-15 and M.mu0 * h / b or M.mu0 * slope, M.mu0 * slope
end

--- Add one simple closed polygon; later regions take material precedence.
-- @tparam table m Unmeshed model.
-- @tparam string material Previously defined material name.
-- @tparam table points Array of {x,y} vertices, in model units, without closure repeat.
-- @tparam[opt] number current Signed Jz in A/m^2; nil inherits the material value.
-- @tparam[opt=0] number angle Magnetisation direction in model-frame degrees.
-- @tparam[opt=0] number mesh_size Local target spacing; zero inherits the global value.
function M.region(m, material, points, current, angle, mesh_size)
    return M.region_contours(m, material, { points }, current, angle, mesh_size, "nonzero")
end

--- Add one compound region; holes preserve the previously declared material.
-- Contours are individually simple, closed implicitly and copied on input.
-- @tparam table m Unmeshed model.
-- @tparam string material Previously defined material name.
-- @tparam table contours Array of vertex arrays in model units.
-- @tparam[opt] number current Signed Jz in A/m^2; nil inherits the material source.
-- @tparam[opt=0] number angle Magnetisation direction in model-frame degrees.
-- @tparam[opt=0] number mesh_size Local target spacing; zero inherits the global value.
-- @tparam[opt="nonzero"] string fill_rule "nonzero" or "even odd", as in PGF.
function M.region_contours(m, material, contours, current, angle, mesh_size, fill_rule)
    check(not m.nodes, "declare regions before meshing")
    check(not m.topology, "imported topology is immutable; edit the FEM source before importing")
    check(m.materials[material], "unknown material " .. material)
    fill_rule = fill_rule or "nonzero"
    check(fill_rule == "nonzero" or fill_rule == "even odd", "unknown fill rule")
    check(type(contours) == "table" and #contours > 0, "a region needs at least one contour")
    for _, points in ipairs(contours) do
        for _, p in ipairs(points) do
            check(finite(p[1]) and finite(p[2]), "invalid vertex")
            check(
                p[1] >= m.options.xmin
                    and p[1] <= m.options.xmax
                    and p[2] >= m.options.ymin
                    and p[2] <= m.options.ymax,
                "region outside domain"
            )
        end
        require("luafemm-mesh").validate(points)
    end
    if current == nil then
        current = m.materials[material].current
    end
    check(finite(current), "invalid current density")
    angle = angle or 0
    check(finite(angle), "invalid magnetization angle")
    mesh_size = mesh_size or 0
    check(finite(mesh_size) and mesh_size >= 0, "invalid local mesh size")
    local owned = copy(contours)
    m.regions[#m.regions + 1] = {
        material = material,
        contours = owned,
        fill_rule = fill_rule,
        current = current or 0,
        angle = angle,
        mx = math.cos(math.rad(angle)),
        my = math.sin(math.rad(angle)),
        mesh_size = mesh_size,
    }
    return m.regions[#m.regions]
end

local function axis(a, b, h, anchors)
    anchors[#anchors + 1] = a
    anchors[#anchors + 1] = b
    table.sort(anchors)
    local unique = {}
    for _, v in ipairs(anchors) do
        if #unique == 0 or v - unique[#unique] > 1e-9 then
            unique[#unique + 1] = v
        end
    end
    local out = { unique[1] }
    for i = 2, #unique do
        local n = math.ceil((unique[i] - unique[i - 1]) / h)
        for j = 1, n do
            out[#out + 1] = unique[i - 1] + (unique[i] - unique[i - 1]) * j / n
        end
    end
    return out
end

--- Mesh and classify a model. Further material/region declarations are invalid.
-- @tparam table m Model.
-- @treturn table The same model with nodes, triangles, free degrees and potential.
function M.mesh(m)
    if m.nodes then
        return m
    end
    require("luafemm-ideal").compile(m)
    local cache = require("luafemm-cache")
    if cache.restore(m, false) then
        return m
    end
    if m.options.mesher == "delaunay" then
        local points, elements, stats = require("luafemm-mesh").generate(m.options, m.regions)
        if m.topology then
            points, elements, stats =
                require("luafemm-topology").carve(m.topology, points, elements, stats)
        end
        m.mesh_stats = stats
        M.import_mesh(m, points, elements)
        cache.store(m)
        return m
    end
    local o = m.options
    local xs, ys = {}, {}
    for _, r in ipairs(m.regions) do
        for _, points in ipairs(r.contours) do
            for i, p in ipairs(points) do
                local q = points[i % #points + 1]
                check(
                    abs(p[1] - q[1]) < 1e-8 or abs(p[2] - q[2]) < 1e-8,
                    "grid mesher requires axis-aligned edges; choose mesher=delaunay"
                )
                xs[#xs + 1] = p[1]
                ys[#ys + 1] = p[2]
            end
        end
    end
    xs = axis(o.xmin, o.xmax, o.h, xs)
    ys = axis(o.ymin, o.ymax, o.h, ys)
    local nodes, tri, boundary = {}, {}, {}
    local nx, ny = #xs, #ys
    check(nx * ny <= 150000, "mesh exceeds the vertex limit; increase mesh size")
    for j, y in ipairs(ys) do
        for i, x in ipairs(xs) do
            local id = #nodes + 1
            nodes[id] = { x * o.unit, y * o.unit }
            if i == 1 or i == nx or j == 1 or j == ny then
                boundary[id] = true
            end
        end
    end
    local function triangle(ids)
        tri[#tri + 1] = ids
    end
    for j = 1, ny - 1 do
        for i = 1, nx - 1 do
            local a = (j - 1) * nx + i
            local b = a + 1
            local c = a + nx
            local d = c + 1
            -- Alternating diagonals reduce directional bias.
            if (i + j) % 2 == 0 then
                triangle({ a, b, d })
                triangle({ a, d, c })
            else
                triangle({ a, b, c })
                triangle({ b, d, c })
            end
        end
    end
    m.xs, m.ys = xs, ys
    M.prepare_mesh(m, nodes, tri, boundary)
    cache.store(m)
    return m
end

--- Internal import of an audited mesh, including its tagged exterior boundary.
-- Nodes are reordered for the IC(0) preconditioner; all mesh tags are remapped.
-- @tparam table m Model with mesh_stats from luafemm-mesh.generate.
-- @tparam table points Coordinates in model units.
-- @tparam table elements Counterclockwise triples of one-based vertex indices.
function M.import_mesh(m, points, elements)
    local order = {}
    for i = 1, #points do
        order[i] = i
    end
    table.sort(order, function(a, b)
        if points[a][2] == points[b][2] then
            return points[a][1] < points[b][1]
        end
        return points[a][2] < points[b][2]
    end)
    local map, nodes, tri = {}, {}, {}
    for id, old in ipairs(order) do
        map[old] = id
        nodes[id] = { points[old][1] * m.options.unit, points[old][2] * m.options.unit }
    end
    for _, el in ipairs(elements) do
        tri[#tri + 1] = { map[el[1]], map[el[2]], map[el[3]] }
    end
    for _, segment in ipairs(m.mesh_stats.segments) do
        segment[1], segment[2] = map[segment[1]], map[segment[2]]
    end
    local boundary = {}
    for old in pairs(m.mesh_stats.boundary_nodes) do
        boundary[map[old]] = true
    end
    m.mesh_stats.boundary_nodes = boundary
    return M.prepare_mesh(m, nodes, tri, boundary)
end

--- Internal reconstruction shared by fresh meshes and passive cache records.
-- Coordinates are in SI; vertex/element order is preserved exactly so that
-- saved potentials, interface sampling and contour connectivity remain stable.
-- Material references and excitations always come from the current model.
-- @tparam table m Model with current physical declarations.
-- @tparam table nodes SI coordinate pairs, without prepared degrees of freedom.
-- @tparam table elements Counterclockwise triples of one-based node indices.
-- @tparam table boundary Set of exterior node indices.
function M.prepare_mesh(m, nodes, elements, boundary)
    local o = m.options
    local fixed, edges, gauge, gauges, components, representatives, loads
    if o.field_model == "ideal" then
        fixed, representatives, loads, m.ideal_loops =
            require("luafemm-ideal").prepare(m, nodes, elements)
        edges = {}
    else
        fixed, edges, gauge, gauges, components =
            boundary_conditions.prepare(m, nodes, elements, boundary)
    end
    local tri, free, A = {}, {}, {}
    for id, p in ipairs(nodes) do
        p.dof = nil
        if
            fixed[id] == nil
            and (not representatives or not representatives[id] or representatives[id] == id)
        then
            free[#free + 1] = id
            p.dof = #free
        end
        A[id] = fixed[id] or 0
    end
    if representatives then
        for id, first in pairs(representatives) do
            nodes[id].dof = nodes[first].dof
        end
        m.ideal_loads = {}
        for id, value in pairs(loads) do
            m.ideal_loads[nodes[id].dof] = value
        end
    end
    local function triangle(ids)
        local p, q, r = nodes[ids[1]], nodes[ids[2]], nodes[ids[3]]
        local twice = (q[1] - p[1]) * (r[2] - p[2]) - (r[1] - p[1]) * (q[2] - p[2])
        check(twice > 0, "nonpositive triangle area")
        local gx = { (q[2] - r[2]) / twice, (r[2] - p[2]) / twice, (p[2] - q[2]) / twice }
        local gy = { (r[1] - q[1]) / twice, (p[1] - r[1]) / twice, (q[1] - p[1]) / twice }
        local x, y = (p[1] + q[1] + r[1]) / 3, (p[2] + q[2] + r[2]) / 3
        local material, current, mx, my = "air", 0, 1, 0
        -- Last declared region wins, like a paint operation.
        for _, reg in ipairs(m.regions) do
            if geometry.contains(reg, x / o.unit, y / o.unit) then
                material, current, mx, my = reg.material, reg.current, reg.mx, reg.my
            end
        end
        local face
        if m.topology then
            face = require("luafemm-topology").locate(m.topology, x / o.unit, y / o.unit)
            check(face and not face.hole, "triangle outside computational domain")
            material, current, mx, my = face.material, face.current, face.mx, face.my
        end
        if o.source then
            current = o.source(x, y)
        end
        local dofs = {}
        for k, id in ipairs(ids) do
            dofs[k] = nodes[id].dof or 0
        end
        tri[#tri + 1] = {
            ids = ids,
            dofs = dofs,
            gx = gx,
            gy = gy,
            area = twice / 2,
            face = face and face.id,
            material = material,
            p = m.materials[material],
            current = current,
            mx = mx,
            my = my,
        }
    end
    for _, ids in ipairs(elements) do
        triangle(ids)
    end
    boundary_conditions.compatibility({
        triangles = tri,
        boundary_edges = edges,
        gauge_node = gauge,
        gauge_nodes = gauges,
        components = components,
    })
    m.nodes, m.triangles, m.free, m.A = nodes, tri, free, A
    m.exterior_nodes, m.fixed_values = boundary, fixed
    m.boundary_edges, m.gauge_node = edges, gauge
    m.gauge_nodes, m.components = gauges, components
    return m
end

local function assemble(m, A, tangent)
    local n = #m.free
    local R = zeros(n)
    local K = tangent and {} or nil
    if K then
        for i = 1, n do
            K[i] = {}
        end
    end
    for _, e in ipairs(m.triangles) do
        local ax, ay = 0, 0
        for i = 1, 3 do
            ax = ax + e.gx[i] * A[e.ids[i]]
            ay = ay + e.gy[i] * A[e.ids[i]]
        end
        local b = sqrt(ax * ax + ay * ay)
        local nu, d = M.constitutive(e.p, b)
        local v = {}
        for i = 1, 3 do
            v[i] = ax * e.gx[i] + ay * e.gy[i]
        end
        local correction = b > 1e-15 and (d - nu) / (b * b) or 0
        for i = 1, 3 do
            local row = e.dofs[i]
            if row > 0 then
                -- H = nu(B) B - Hc m; test curl(Ni ez)=(Ni,y,-Ni,x).
                local permanent = M.mu0 * e.p.hc * (e.mx * e.gy[i] - e.my * e.gx[i])
                R[row] = R[row] + e.area * (nu * v[i] - permanent - M.mu0 * e.current / 3)
                if K then
                    for j = 1, 3 do
                        local col = e.dofs[j]
                        if col > 0 then
                            K[row][col] = (K[row][col] or 0)
                                + e.area
                                    * (nu * (e.gx[i] * e.gx[j] + e.gy[i] * e.gy[j]) + correction * v[i] * v[j])
                        end
                    end
                end
            end
        end
    end
    if m.ideal_loads then
        for row, value in pairs(m.ideal_loads) do
            R[row] = R[row] - M.mu0 * value
        end
    else
        boundary_conditions.assemble(m, A, R, K, M.mu0)
    end
    return R, K
end

-- IC(0)-preconditioned conjugate gradients. On a non-M-matrix Newton
-- tangent, retry the entire factorization with a diagonal shift instead
-- of clipping individual pivots (which can create an unstable factor).
local function pcg(K, rhs, tol)
    local n = #rhs
    local L, diag, indices = {}, {}, {}
    for i = 1, n do
        local cols = {}
        for j in pairs(K[i]) do
            if j < i then
                cols[#cols + 1] = j
            end
        end
        table.sort(cols)
        indices[i] = cols
    end
    local function factor(shift)
        L, diag = {}, {}
        for i = 1, n do
            local cols = indices[i]
            local li = {}
            L[i] = li
            for _, j in ipairs(cols) do
                local s = K[i][j]
                for _, k in ipairs(cols) do
                    if k >= j then
                        break
                    end
                    s = s - li[k] * (L[j][k] or 0)
                end
                li[j] = s / diag[j]
            end
            local d = K[i][i] * (1 + shift)
            for _, j in ipairs(cols) do
                d = d - li[j] ^ 2
            end
            if not finite(d) or d < K[i][i] * 0.01 then
                return false
            end
            diag[i] = sqrt(d)
        end
        return true
    end
    local factored = false
    for _, shift in ipairs({ 0, 0.01, 0.1, 1, 10, 100 }) do
        if factor(shift) then
            factored = true
            break
        end
    end
    check(factored, "incomplete Cholesky failed even with diagonal shift")
    local function pre(r)
        local z = {}
        for i = 1, n do
            local s = r[i]
            for _, j in ipairs(indices[i]) do
                s = s - L[i][j] * z[j]
            end
            z[i] = s / diag[i]
        end
        for i = n, 1, -1 do
            z[i] = z[i] / diag[i]
            for _, j in ipairs(indices[i]) do
                z[j] = z[j] - L[i][j] * z[i]
            end
        end
        return z
    end
    local x = zeros(n)
    local r = {}
    for i = 1, n do
        r[i] = rhs[i]
    end
    local target = max(norm(rhs) * tol, 1e-30)
    if norm(r) <= target then
        return x, 0
    end
    local z = pre(r)
    local p = {}
    for i = 1, n do
        p[i] = z[i]
    end
    local rz = dot(r, z)
    for iteration = 1, 4000 do
        local q = {}
        for i = 1, n do
            local s = 0
            for j, v in pairs(K[i]) do
                s = s + v * p[j]
            end
            q[i] = s
        end
        local pq = dot(p, q)
        check(pq > 0, "non-positive CG curvature")
        local alpha = rz / pq
        for i = 1, n do
            x[i] = x[i] + alpha * p[i]
            r[i] = r[i] - alpha * q[i]
        end
        if norm(r) <= target then
            return x, iteration
        end
        z = pre(r)
        local next_rz = dot(r, z)
        local beta = next_rz / rz
        for i = 1, n do
            p[i] = z[i] + beta * p[i]
        end
        rz = next_rz
    end
    error("luafemm: CG failed to converge")
end

--- Solve with damped Newton and IC(0)-preconditioned conjugate gradients.
-- The mesh is generated on demand. Failure to converge raises an error.
-- @tparam table m Model.
-- @treturn table The same model with A, element fields and convergence statistics.
function M.solve(m)
    if m.stats then
        return m
    end
    require("luafemm-ideal").compile(m)
    local cache = require("luafemm-cache")
    if cache.restore(m, true) then
        return m
    end
    if not m.nodes then
        M.mesh(m)
    end
    local start = os.clock()
    local A = m.A
    local initial
    local total_cg = 0
    for iteration = 0, m.options.max_newton do
        local R, K = assemble(m, A, true)
        local residual = norm(R)
        initial = initial or max(residual, 1e-30)
        if m.options.log then
            m.options.log(
                string.format("Newton %d: relative residual %.3g", iteration, residual / initial)
            )
        end
        if residual <= m.options.tolerance * initial or residual < 1e-25 then
            m.A = A
            m.stats = {
                newton = iteration,
                cg = total_cg,
                residual = residual / initial,
                seconds = os.clock() - start,
                nodes = #m.nodes,
                triangles = #m.triangles,
            }
            M.update_fields(m)
            cache.store(m)
            return m
        end
        check(iteration < m.options.max_newton, "Newton failed to converge")
        for i = 1, #R do
            R[i] = -R[i]
        end
        local delta, it = pcg(K, R, 1e-9)
        total_cg = total_cg + it
        local alpha = 1
        local accepted = false
        for backtrack = 0, 24 do
            local trial = {}
            for i = 1, #A do
                trial[i] = A[i]
            end
            for id, node in ipairs(m.nodes) do
                if node.dof then
                    trial[id] = A[id] + alpha * delta[node.dof]
                end
            end
            local r = assemble(m, trial, false)
            if norm(r) < residual * (1 - 1e-4 * alpha) then
                A = trial
                accepted = true
                break
            end
            alpha = alpha / 2
        end
        check(accepted, "Newton line search failed")
    end
end

--- Refresh piecewise-constant element fields and aggregate statistics from A.
-- @tparam table m Solved model with statistics allocated by solve.
function M.update_fields(m)
    local peak = 0
    local current = 0
    for _, e in ipairs(m.triangles) do
        local bx, by = 0, 0
        for i = 1, 3 do
            bx = bx + e.gy[i] * m.A[e.ids[i]]
            by = by - e.gx[i] * m.A[e.ids[i]]
        end
        e.bx, e.by, e.b = bx, by, sqrt(bx * bx + by * by)
        local nu = M.constitutive(e.p, e.b)
        e.hx, e.hy = nu * bx / M.mu0 - e.p.hc * e.mx, nu * by / M.mu0 - e.p.hc * e.my
        peak = max(peak, e.b)
        current = current + e.current * e.area
    end
    m.stats.peak = peak
    m.stats.net_current = current
    m.solution_signature = require("luafemm-cache").signature(m)
end

--- Evaluate a solved model, using the first containing triangle at an interface.
-- @tparam table m Solved model.
-- @tparam number x Model coordinate (not metres unless options.unit=1).
-- @tparam number y Model coordinate.
-- @return Bx, By (T), Az (T m), Hx, Hy (A/m); error outside the domain.
function M.sample(m, x, y)
    check(m and m.stats, "solve before sampling the field")
    check(finite(x) and finite(y), "sample coordinates must be finite")
    x = x * m.options.unit
    y = y * m.options.unit
    for _, e in ipairs(m.triangles) do
        local p = m.nodes[e.ids[1]]
        local dx, dy = x - p[1], y - p[2]
        local b = e.gx[2] * dx + e.gy[2] * dy
        local c = e.gx[3] * dx + e.gy[3] * dy
        if b >= -1e-10 and c >= -1e-10 and b + c <= 1 + 1e-10 then
            return e.bx,
                e.by,
                (1 - b - c) * m.A[e.ids[1]] + b * m.A[e.ids[2]] + c * m.A[e.ids[3]],
                e.hx,
                e.hy
        end
    end
    error("luafemm: sample outside domain")
end

-- Marching triangles, stitched by mesh-edge identifiers, not rounded
-- coordinates. Orient every connected contour with B=(Ay,-Ax).
--- Extract oriented potential contours without smoothing the finite-element field.
-- @tparam table m Solved model.
-- @tparam number count Number of evenly spaced levels if levels is absent.
-- @tparam[opt] table levels Explicit potential levels in T m.
-- @treturn table Paths with x/y points in model units and a level property.
function M.contours(m, count, levels)
    check(m and m.stats, "solve before tracing field lines")
    if not levels then
        check(
            finite(count) and count >= 1 and count == math.floor(count),
            "line count must be a positive integer"
        )
        local lo, hi = math.huge, -math.huge
        for id, a in ipairs(m.A) do
            -- Confined circuits use the transported-flux range. Permanent
            -- magnets can also create local recirculation extrema inside a part.
            if m.options.field_model ~= "ideal" or m.exterior_nodes[id] then
                lo = math.min(lo, a)
                hi = max(hi, a)
            end
        end
        levels = {}
        if hi - lo < 1e-20 then
            return {}
        end
        for i = 1, count do
            levels[i] = lo + (hi - lo) * (i - 0.5) / count
        end
    end
    local paths = {}
    for _, level in ipairs(levels) do
        local segments, adj = {}, {}
        for _, e in ipairs(m.triangles) do
            local hits = {}
            for i = 1, 3 do
                local a, b = e.ids[i], e.ids[i % 3 + 1]
                local va, vb = m.A[a], m.A[b]
                if (va > level) ~= (vb > level) then
                    local t = (level - va) / (vb - va)
                    local p, q = m.nodes[a], m.nodes[b]
                    hits[#hits + 1] = {
                        x = (p[1] + t * (q[1] - p[1])) / m.options.unit,
                        y = (p[2] + t * (q[2] - p[2])) / m.options.unit,
                        key = math.min(a, b) .. ":" .. max(a, b),
                    }
                end
            end
            if #hits == 2 then
                if (hits[2].x - hits[1].x) * e.bx + (hits[2].y - hits[1].y) * e.by < 0 then
                    hits[1], hits[2] = hits[2], hits[1]
                end
                local id = #segments + 1
                segments[id] = hits
                for _, p in ipairs(hits) do
                    adj[p.key] = adj[p.key] or {}
                    table.insert(adj[p.key], id)
                end
            end
        end
        local used = {}
        local function trace(id, start)
            local path = { level = level }
            local key = start
            while id and not used[id] do
                used[id] = true
                local seg = segments[id]
                local a, b = seg[1], seg[2]
                if a.key ~= key then
                    a, b = b, a
                end
                if #path == 0 then
                    path[1] = a
                    path.forward = (a == seg[1])
                end
                path[#path + 1] = b
                key = b.key
                id = nil
                for _, other in ipairs(adj[key]) do
                    if not used[other] then
                        id = other
                        break
                    end
                end
            end
            if not path.forward then
                local rev = { level = level }
                for i = #path, 1, -1 do
                    rev[#rev + 1] = path[i]
                end
                path = rev
            end
            if #path > 1 then
                paths[#paths + 1] = path
            end
        end
        for id, s in ipairs(segments) do
            if not used[id] then
                if #adj[s[1].key] == 1 then
                    trace(id, s[1].key)
                elseif #adj[s[2].key] == 1 then
                    trace(id, s[2].key)
                end
            end
        end
        for id, s in ipairs(segments) do
            if not used[id] then
                trace(id, s[1].key)
            end
        end
    end
    return paths
end

--- Parse the comma-separated numeric a/b pairs used by the TeX interface.
-- @tparam string s Pair list; an empty string returns an empty table.
-- @treturn table Array of two-number arrays; malformed input is an error.
function M.parse_pairs(s)
    local out = {}
    for pair in s:gmatch("[^,]+") do
        local a, b = pair:match("^%s*([^%s/]+)%s*/%s*([^%s/]+)%s*$")
        a, b = tonumber(a), tonumber(b)
        check(finite(a) and finite(b), "invalid numeric pair")
        out[#out + 1] = { a, b }
    end
    return out
end

--- Add an ideal-window excitation in model coordinates, with signed NI in amperes.
function M.excitation(m, x, y, turns)
    return require("luafemm-ideal").excitation(m, x, y, turns)
end

--- Flux through an oriented section, along its right normal, in webers.
-- A positive constant problem depth is required; lengths use model units.
-- @tparam table m Source model, solved on demand.
-- @tparam table points Section polyline, left unchanged.
-- @tparam[opt=0] number snap_tolerance Endpoint projection distance in model units.
-- @treturn number Signed flux in Wb.
function M.flux(m, points, snap_tolerance)
    return require("luafemm-ideal").flux(m, points, snap_tolerance)
end

--- Import a magnetic FEMM problem without solving it.
function M.import_fem(file, options)
    local codec = require("luafemm-fem")
    return codec.model(codec.read(file), options)
end

--- Export the computational geometry without triggering a solve.
function M.export_fem(m, file)
    return require("luafemm-fem").write(m, file)
end

--- Export the existing converged mesh and potential to a static answer file.
function M.export_ans(m, file)
    return require("luafemm-ans").write(m, file)
end

return M
