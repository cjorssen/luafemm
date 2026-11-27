-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Planar and axisymmetric scalar Poisson elements with prescribed dielectric and electrode data.
-- The equation is -div(epsilon grad V)=rho. Stiffness and loads are divided
-- by epsilon0 for conditioning; lengths remain model coordinates during mesh
-- assembly, with explicit SI factors on volume and edge loads. All returned
-- field/integral values are SI. No magnetic constitutive law is reused.
-- @module luafemm-electric-fe
local F = {}
local E = require("luafemm-electrostatic")
local P = require("luafemm-physical")
local G = require("luafemm-geometry")
local mesh = require("luafemm-mesh")
local metric = require("luafemm-electric-metric")
local abs, sqrt, min, max = math.abs, math.sqrt, math.min, math.max
local function check(v, s)
    assert(v, "luafemm electric FE: " .. s)
end
-- Empty TeX strings have the same meaning as omitted Lua options.
local function option(v, d)
    if v == nil or v == "" then
        return d
    end
    return v
end
local function number(v, name, d)
    return E.number(option(v, d), name)
end
local function positive(v, name, d)
    return E.positive(option(v, d), name)
end
local function key(a, b)
    if a > b then
        a, b = b, a
    end
    return a .. ":" .. b
end
local function mutable(m)
    check(m.fe, "finite-element model required")
    check(not m.fe.native_locked, "imported model owns its physical declarations")
    check(
        not m.fe.nodes,
        "declare all electric regions and boundaries before meshing or evaluation"
    )
    return m.fe
end
function F.configure(m, o)
    for _, k in ipairs({ "surface_density", "voltage", "reference_radius" }) do
        check(
            o[k] == nil or o[k] == "" or k == "reference_radius" and o.model == "source correction",
            k .. " is not a planar FE model input; use a physical path"
        )
    end
    local b = m.bounds
    local geometry = option(o.problem_type, "planar")
    check(geometry == "planar" or geometry == "axisymmetric", "unknown electric problem type")
    local axial = geometry == "axisymmetric"
    check(not axial or b[1] >= 0, "axisymmetric domain requires r >= 0")
    check(not axial or o.depth == nil or o.depth == "", "depth is not an axisymmetric input")
    local boundary = option(o.outer_boundary, "voltage")
    check(boundary == "voltage" or boundary == "insulated", "unknown outer boundary")
    local tolerance = positive(o.solve_tolerance, "solver tolerance", 1e-10)
    check(tolerance < 1, "solver tolerance must be less than one")
    m.fe = {
        axisymmetric = axial,
        regions = {},
        paths = {},
        electrodes = {},
        bounds = P.copy(b),
        density = number(o.volume_density, "volume charge density", 0),
        boundary_tolerance = number(o.boundary_tolerance, "boundary snap tolerance", 0.0001),
        outer_boundary = boundary,
        outer_potential = number(o.outer_potential, "boundary potential", 0),
        depth = axial and 1 or positive(o.depth, "depth", 1) * m.unit,
        options = {
            xmin = b[1],
            ymin = b[2],
            xmax = b[3],
            ymax = b[4],
            h = positive(o.mesh_size, "mesh size", 0.25),
            mesh_tolerance = number(o.mesh_tolerance, "mesh tolerance", 0),
            min_angle = number(o.min_angle, "minimum angle", 0),
            max_area = number(o.max_area, "maximum area", 0),
            max_steiner = E.integer(
                option(o.max_steiner, 10000),
                "maximum Steiner points",
                0,
                1000000
            ),
        },
        tolerance = tolerance,
        max_iterations = E.integer(
            option(o.max_iterations, 10000),
            "maximum iterations",
            1,
            1000000
        ),
    }
    check(
        m.fe.boundary_tolerance >= 0
            and m.fe.boundary_tolerance < 0.01 * min(b[3] - b[1], b[4] - b[2]),
        "boundary snap tolerance must be nonnegative and smaller than 1% of the domain span"
    )
    local q = m.fe.options
    check(
        q.mesh_tolerance >= 0 and q.min_angle >= 0 and q.min_angle < 34 and q.max_area >= 0,
        "invalid mesh quality target"
    )
    local d = m.descriptor
    d.backend = "finite elements"
    d.geometry = geometry
    d.symmetry = axial and "revolution about r=0; x=r, y=z"
        or "translation perpendicular to the section"
    d.reference = "prescribed potential, positive mixed term or antiperiodic constraint"
    d.medium = "piecewise linear symmetric positive dielectric tensor"
    d.relative_permittivity_tensor = P.copy(m.tensor)
    if m.tensor.xx ~= m.tensor.yy or m.tensor.xy ~= 0 then
        d.epsilon = nil
    end
    d.validity = "finite rectangular domain; raw element gradients, deterministic one-sided samples"
    d.depth = not axial and m.fe.depth or nil
    d.units.volume_charge_density = "C/m^3"
    d.units.surface_charge_density = "C/m^2"
end
local function vertices(m, contours, closed)
    check(type(contours) == "table" and #contours > 0, "empty physical path")
    local owned = P.copy(contours)
    local b = m.fe.bounds
    -- Coordinate snapping here is restricted to roundoff, not a geometric repair.
    local tol = 1e-10 * max(b[3] - b[1], b[4] - b[2])
    for _, points in ipairs(owned) do
        check(#points >= (closed and 3 or 2), "physical path has too few vertices")
        for _, p in ipairs(points) do
            p[1], p[2] = E.number(p[1], "vertex x"), E.number(p[2], "vertex y")
            check(
                p[1] >= b[1] - tol
                    and p[1] <= b[3] + tol
                    and p[2] >= b[2] - tol
                    and p[2] <= b[4] + tol,
                "physical path outside FE domain"
            )
            p[1] = max(b[1], min(b[3], p[1]))
            p[2] = max(b[2], min(b[4], p[2]))
        end
        if closed then
            mesh.validate(points)
        end
    end
    return owned
end
--- Last declared dielectric wins inside its PGF fill rule; holes retain earlier data.
function F.region(m, contours, o)
    local f = mutable(m)
    o = o or {}
    local D = require("luafemm-dielectrics")
    local dielectric = o.material
            and o.material ~= ""
            and o.material ~= "custom"
            and D.get(o.material)
        or nil
    local rule = option(o.fill_rule, "nonzero")
    check(rule == "nonzero" or rule == "even odd", "unknown fill rule")
    local h = number(o.mesh_size, "local mesh size", 0)
    check(h >= 0, "local mesh size must be nonnegative")
    f.regions[#f.regions + 1] = {
        material = dielectric and dielectric.id or nil,
        tensor = D.tensor(o, m.tensor),
        contours = vertices(m, contours, true),
        fill_rule = rule,
        mesh_size = h,
        epsilon = positive(
            o.epsilon_r,
            "relative permittivity",
            dielectric and dielectric.ex or m.epsilon / E.epsilon0
        ) * E.epsilon0,
        density = number(o.density, "volume charge density", dielectric and dielectric.qv or 0),
    }
end
--- One open/closed physical polyline: equipotential electrode, prescribed sheet
-- charge or exterior condition. A solid electrode is one simple closed loop.
function F.path(m, points, closed, kind, o)
    local f = mutable(m)
    o = o or {}
    check(
        kind == "electrode" or kind == "charge" or kind == "boundary",
        "unknown physical path kind"
    )
    local owned = P.copy(points)
    if
        closed
        and #owned > 1
        and owned[1][1] == owned[#owned][1]
        and owned[1][2] == owned[#owned][2]
    then
        table.remove(owned)
    end
    owned = vertices(m, { owned }, closed)[1]
    if f.axisymmetric and kind ~= "boundary" then
        local positive_radius = false
        for _, p in ipairs(owned) do
            positive_radius = positive_radius or p[1] > 0
        end
        check(positive_radius, "an axis-only path has zero revolved area")
    end
    local h = number(o.mesh_size, "local mesh size", 0)
    check(h >= 0, "local mesh size must be nonnegative")
    local path = { points = owned, closed = closed, kind = kind, mesh_size = h }
    if kind == "electrode" then
        path.condition = option(o.condition, "voltage")
        check(
            path.condition == "voltage" or path.condition == "charge",
            "unknown electrode condition"
        )
        if path.condition == "voltage" then
            path.potential = number(o.potential, "electrode potential")
            check(option(o.charge, nil) == nil, "charge requires a charge-controlled electrode")
        else
            check(
                option(o.potential, nil) == nil,
                "a charge-controlled electrode cannot prescribe potential"
            )
            path.charge = number(o.charge, "electrode charge", 0)
        end
        path.name = option(o.name, nil)
        if not path.name then
            local index = #f.electrodes + 1
            repeat
                path.name = "electrode " .. index
                index = index + 1
            until not f.electrodes[path.name]
        end
        check(
            type(path.name) == "string" and path.name:match("%S"),
            "electrode needs a nonempty name"
        )
        check(path.name ~= "outer boundary", "reserved electrode name: outer boundary")
        local previous = f.electrodes[path.name]
        check(
            not previous
                or (
                    previous.condition == path.condition
                    and previous.potential == path.potential
                    and previous.charge == path.charge
                ),
            "inconsistent data for connected electrode parts"
        )
        path.solid = closed and option(o.interior, "conductor") == "conductor"
        check(
            option(o.interior, "conductor") == "conductor" or o.interior == "dielectric",
            "unknown electrode interior"
        )
        if not previous then
            f.electrodes[path.name] = path
            f.electrodes[#f.electrodes + 1] = path.name
        end
    elseif kind == "charge" then
        path.density = number(o.density, "surface charge density")
    else
        path.condition = option(o.condition, "insulated")
        check(
            path.condition == "voltage"
                or path.condition == "insulated"
                or path.condition == "displacement"
                or path.condition == "mixed"
                or path.condition == "periodic"
                or path.condition == "antiperiodic",
            "unknown electric boundary condition"
        )
        if path.condition == "voltage" then
            path.potential = number(o.potential, "boundary potential")
        end
        path.displacement = (path.condition == "displacement" or path.condition == "mixed")
                and number(
                    o.displacement,
                    "normal displacement",
                    path.condition == "mixed" and 0 or nil
                )
            or 0
        path.pair = o.pair
        path.coefficient = path.condition == "mixed"
                and number(o.coefficient, "mixed coefficient", 0)
            or 0
        check(path.coefficient >= 0, "mixed coefficient must be nonnegative")
    end
    f.paths[#f.paths + 1] = path
end
local function onsegment(p, a, b, tol)
    local dx, dy = b[1] - a[1], b[2] - a[2]
    local l = dx * dx + dy * dy
    if l == 0 then
        return false
    end
    local t = ((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / l
    return t >= -tol / sqrt(l)
        and t <= 1 + tol / sqrt(l)
        and abs(dx * (p[2] - a[2]) - dy * (p[1] - a[1])) <= tol * sqrt(l)
end
local function pairs_of(path)
    local out = {}
    for i = 2, #path.points do
        out[#out + 1] = { path.points[i - 1], path.points[i] }
    end
    if path.closed then
        out[#out + 1] = { path.points[#path.points], path.points[1] }
    end
    return out
end
--- Generate a shared constrained mesh, then classify elements independently of magnetics.
function F.mesh(m)
    local f = m.fe
    if f.nodes then
        if m.hybrid then
            require("luafemm-electric-hybrid").check_mesh(m)
        end
        return f
    end
    if m.hybrid then
        require("luafemm-electric-hybrid").validate(m)
    end
    local constraints = P.copy(f.regions)
    for _, p in ipairs(f.paths) do
        -- Each two-vertex contour contributes a segment twice; the mesher deduplicates
        -- the reversed copy. It does not close an open polyline across its endpoints.
        for _, s in ipairs(pairs_of(p)) do
            constraints[#constraints + 1] =
                { contours = { s }, fill_rule = "nonzero", mesh_size = p.mesh_size }
        end
    end
    local nodes, ids, stats
    if f.extended_mesh then
        local saved = f.extended_mesh
        nodes, ids, stats = P.copy(saved.nodes), P.copy(saved.elements), P.copy(saved.stats)
        f.periodic_pairs = P.copy(saved.periodic_pairs)
        f.periodic_nodes = {}
        for _, pair in ipairs(f.periodic_pairs) do
            f.periodic_nodes[pair[1]], f.periodic_nodes[pair[2]] = true, true
        end
    elseif f.native then
        nodes, ids, stats = require("luafemm-electric-native").generate(m, constraints)
    else
        nodes, ids, stats = require("luafemm-electric-periodic").generate(f, constraints)
    end
    local triangles, edges, edge_list = {}, {}, {}
    for ti, t in ipairs(ids) do
        local a, b, c = nodes[t[1]], nodes[t[2]], nodes[t[3]]
        local twice = (b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])
        check(twice > 0, "nonpositive element area")
        local tri = {
            ids = t,
            area = twice / 2,
            gx = { (b[2] - c[2]) / twice, (c[2] - a[2]) / twice, (a[2] - b[2]) / twice },
            gy = { (c[1] - b[1]) / twice, (a[1] - c[1]) / twice, (b[1] - a[1]) / twice },
            epsilon = m.epsilon,
            tensor = m.tensor,
            density = f.density,
            region = 0,
        }
        local x, y = (a[1] + b[1] + c[1]) / 3, (a[2] + b[2] + c[2]) / 3
        for ri, r in ipairs(f.regions) do
            if G.contains(r, x, y) then
                tri.epsilon, tri.density, tri.region = r.epsilon, r.density, ri
                tri.tensor = r.tensor
            end
        end
        if f.native then
            local label = f.native.labels[f.native_labels[ti]]
            local material = f.native.materials[label[3]]
            tri.tensor = { xx = material.ex, yy = material.ey, xy = 0 }
            tri.epsilon, tri.density, tri.region =
                E.epsilon0 * material.ex, material.qv, f.native_labels[ti]
        end
        for _, p in ipairs(f.paths) do
            local name = p.name
            if p.solid and G.contains({ contours = { p.points }, fill_rule = "nonzero" }, x, y) then
                check(not tri.electrode or tri.electrode == name, "overlapping solid electrodes")
                tri.electrode = name
            end
        end
        if f.axisymmetric and (a[1] == 0 or b[1] == 0 or c[1] == 0) then
            check(
                abs(tri.tensor.xy) < 1e-12 * max(tri.tensor.xx, tri.tensor.yy),
                "a rotated dielectric touching the regular axis must have zero radial-axial coupling"
            )
        end
        triangles[ti] = tri
        for j = 1, 3 do
            local u, v = t[j], t[j % 3 + 1]
            local k = key(u, v)
            if not edges[k] then
                edges[k] = { u, v, triangles = {} }
                edge_list[#edge_list + 1] = edges[k]
            end
            local e = edges[k]
            e.triangles[#e.triangles + 1] = ti
            check(#e.triangles <= 2, "nonmanifold solution mesh")
        end
    end
    local fixed, owners, explicit_boundary, floating = {}, {}, {}, {}
    -- Named electrodes may override the default box voltage, but not a user
    -- boundary voltage meeting them at only one endpoint.
    local function fix(id, value, name, explicit)
        local conductor = f.electrodes[name]
        if conductor and conductor.condition == "charge" then
            check(
                fixed[id] == nil and (not owners[id] or owners[id] == name),
                "charge-controlled electrode touches another electrode or voltage boundary"
            )
            floating[id], owners[id] = name, name
            return
        end
        check(not floating[id], "voltage boundary touches a charge-controlled electrode")
        if
            fixed[id] ~= nil
            and (explicit_boundary[id] or owners[id] ~= "outer boundary" or not explicit)
        then
            check(
                abs(fixed[id] - value) <= 1e-12 * max(1, abs(value), abs(fixed[id])),
                "incompatible potentials on touching electrodes/boundaries"
            )
            check(
                owners[id] == name or owners[id] == "outer boundary" or name == "outer boundary",
                "touching named electrodes must be one electrode"
            )
        end
        fixed[id], owners[id] = value, name
    end
    -- Default exterior conditions are overridden by explicitly tagged edge paths.
    local tol = max(
        stats.mesh_tolerance * 8,
        1e-10 * max(f.bounds[3] - f.bounds[1], f.bounds[4] - f.bounds[2])
    )
    for _, e in ipairs(edge_list) do
        e.outer = #e.triangles == 1
        e.loads = 0
        local a, b = nodes[e[1]], nodes[e[2]]
        e.axis = f.axisymmetric and a[1] == 0 and b[1] == 0
        if e.axis then
            e.condition = "insulated"
        elseif e.outer then
            e.condition = f.outer_boundary
            if f.outer_boundary == "voltage" then
                e.potential = f.outer_potential
            end
        end
        for pi, p in ipairs(f.paths) do
            local matched = false
            for _, s in ipairs(pairs_of(p)) do
                if onsegment(a, s[1], s[2], tol) and onsegment(b, s[1], s[2], tol) then
                    matched = true
                    break
                end
            end
            if matched and e.axis and p.kind == "boundary" then
                check(
                    p.condition == "insulated",
                    "the regular axis only accepts an insulated boundary"
                )
            end
            -- Axis closure segments have zero revolved area. They are not
            -- electrode surfaces or charged sheets; adjacent finite-area arcs
            -- and solid conductor interiors still constrain their own nodes.
            if matched and not e.axis then
                if p.kind == "boundary" and p.native_surface and not e.outer then
                    -- The adapter normalizes edition-specific native qs to
                    -- one total displacement jump before element assembly.
                    e.loads = e.loads + p.native_surface
                elseif p.kind == "boundary" then
                    check(e.outer, "electric boundary paths must lie on the active exterior")
                    check(not e.boundary_path, "overlapping electric boundary paths")
                    e.boundary_path = pi
                    e.condition = p.condition
                    e.potential = p.potential
                    e.displacement = p.displacement
                    e.coefficient = p.coefficient
                elseif p.kind == "electrode" then
                    check(not e.electrode or e.electrode == p.name, "overlapping electrode paths")
                    e.electrode = p.name
                else
                    e.loads = e.loads + p.density
                end
            end
        end
        if e.outer and e.condition == "voltage" then
            if e.boundary_path then
                explicit_boundary[e[1]], explicit_boundary[e[2]] = true, true
            end
            fix(e[1], e.potential, "outer boundary")
            fix(e[2], e.potential, "outer boundary")
        end
    end
    for _, e in ipairs(edge_list) do
        if e.electrode then
            local p = f.electrodes[e.electrode]
            fix(e[1], p.potential, p.name, true)
            fix(e[2], p.potential, p.name, true)
        end
    end
    for _, t in ipairs(triangles) do
        if t.electrode then
            for _, id in ipairs(t.ids) do
                fix(id, f.electrodes[t.electrode].potential, t.electrode, true)
            end
        end
    end
    require("luafemm-electric-native").points(m, nodes, fixed, owners, floating, fix)
    local anchored = next(fixed)
    for _, e in ipairs(edge_list) do
        anchored = anchored or (e.condition == "mixed" and e.coefficient > 0)
    end
    check(
        anchored or #f.periodic_pairs > 0,
        "at least one prescribed potential is required; pure Neumann models are not supported"
    )
    for _, e in ipairs(edge_list) do
        local active = false
        for _, ti in ipairs(e.triangles) do
            if not triangles[ti].electrode then
                active = true
            end
        end
        check(active or e.loads == 0, "surface charge inside a solid electrode")
    end
    f.interface_pairs = {}
    for _, e in ipairs(edge_list) do
        if e.loads ~= 0 then
            check(
                not e.outer and not e.electrode,
                "surface charge paths must be internal dielectric interfaces"
            )
            for _, id in ipairs(e.triangles) do
                check(
                    not triangles[id].electrode,
                    "surface charge on a solid electrode is over-specified"
                )
            end
        end
        if e.electrode and e.boundary_path then
            check(
                e.condition ~= "displacement"
                    and e.condition ~= "mixed"
                    and e.condition ~= "periodic"
                    and e.condition ~= "antiperiodic",
                "electrode overlaps a prescribed displacement boundary"
            )
            check(
                e.potential == nil or e.potential == f.electrodes[e.electrode].potential,
                "electrode conflicts with an explicit boundary potential"
            )
        end
        if #e.triangles == 2 then
            local a, b = e.triangles[1], e.triangles[2]
            if
                e.loads ~= 0
                or triangles[a].tensor.xx ~= triangles[b].tensor.xx
                or triangles[a].tensor.xy ~= triangles[b].tensor.xy
                or triangles[a].tensor.yy ~= triangles[b].tensor.yy
            then
                f.interface_pairs[key(a, b)] = true
            end
        end
    end
    f.nodes, f.triangles, f.edges, f.stats, f.fixed, f.owners =
        nodes, triangles, edge_list, stats, fixed, owners
    f.floating = floating
    if m.hybrid then
        require("luafemm-electric-hybrid").check_mesh(m)
    end
    return f
end
--- Solve the reduced SPD system. Keep full variational reactions for charges.
function F.solve(m)
    if m.hybrid then
        return require("luafemm-electric-hybrid").solve(m)
    end
    local f = F.mesh(m)
    if f.solution then
        return m
    end
    local K, b = {}, {}
    for i = 1, #f.nodes do
        K[i], b[i] = {}, 0
    end
    local source = 0
    for _, t in ipairs(f.triangles) do
        if not t.electrode then
            local weight, fractions = metric.triangle(m, t)
            local load = t.density * t.area * m.unit ^ 2 * weight / E.epsilon0
            source = source + load * E.epsilon0
            for i, u in ipairs(t.ids) do
                b[u] = b[u] + load * fractions[i]
                for j, v in ipairs(t.ids) do
                    K[u][v] = (K[u][v] or 0)
                        + t.area
                            * weight
                            * (t.gx[i] * (t.tensor.xx * t.gx[j] + t.tensor.xy * t.gy[j]) + t.gy[i] * (t.tensor.xy * t.gx[j] + t.tensor.yy * t.gy[j]))
                end
            end
        end
    end
    local exterior_flux = 0
    for _, e in ipairs(f.edges) do
        local a, c = f.nodes[e[1]], f.nodes[e[2]]
        local length, fractions = metric.edge(m, a, c)
        local outward = e.outer
                and (e.condition == "displacement" or e.condition == "mixed")
                and e.displacement
            or 0
        if e.condition == "mixed" and e.coefficient > 0 then
            local mass = metric.edge_mass(m, a, c)
            for i = 1, 2 do
                for j = 1, 2 do
                    K[e[i]][e[j]] = (K[e[i]][e[j]] or 0) + e.coefficient * mass[i][j] / E.epsilon0
                end
            end
        end
        local load = (e.loads - outward) * length / E.epsilon0
        for i = 1, 2 do
            b[e[i]] = b[e[i]] + load * fractions[i]
        end
        source = source + e.loads * length
        exterior_flux = exterior_flux + outward * length
    end
    for id, charge in pairs(f.point_loads or {}) do
        b[id] = b[id] + charge / E.epsilon0
        source = source + charge
    end
    local solution, iterations, stats
    if f.imported_solution then
        solution = P.copy(f.imported_solution)
        local floating_values = {}
        for id, value in pairs(f.fixed) do
            check(
                math.abs(solution[id] - value) <= 1e-7 * math.max(1, math.abs(value)),
                "native solution contradicts a prescribed potential"
            )
        end
        for id, name in pairs(f.floating) do
            local value = floating_values[name]
            check(
                not value or math.abs(solution[id] - value) <= 1e-7 * math.max(1, math.abs(value)),
                "native solution violates an equipotential conductor"
            )
            floating_values[name] = solution[id]
        end
        iterations, stats = 0, { relative_residual = 0, imported = true }
    else
        solution, iterations, stats =
            require("luafemm-electric-constraints").solve(f, K, b, E.epsilon0)
    end
    local charges, all_reactions = {}, 0
    for _, name in ipairs(f.electrodes) do
        charges[name] = 0
    end
    charges["outer boundary"] = 0
    local voltages = {}
    for i, owner in pairs(f.owners) do
        voltages[owner] = solution[i]
        local r = -b[i]
        for j, k in pairs(K[i]) do
            r = r + k * solution[j]
        end
        r = r * E.epsilon0
        charges[f.owners[i]] = charges[f.owners[i]] + r
        all_reactions = all_reactions + r
    end
    local periodic_reaction = 0
    for i in pairs(f.periodic_nodes) do
        if not f.owners[i] then
            local r = -b[i]
            for j, k in pairs(K[i]) do
                r = r + k * solution[j]
            end
            periodic_reaction = periodic_reaction + r * E.epsilon0
        end
    end
    f.periodic_flux = -periodic_reaction
    local boundary_energy = 0
    for _, e in ipairs(f.edges) do
        if e.condition == "mixed" then
            local a, c = f.nodes[e[1]], f.nodes[e[2]]
            local measure, fractions = metric.edge(m, a, c)
            local mass = metric.edge_mass(m, a, c)
            for i = 1, 2 do
                exterior_flux = exterior_flux
                    + e.coefficient * measure * fractions[i] * solution[e[i]]
                for j = 1, 2 do
                    boundary_energy = boundary_energy
                        + 0.5 * e.coefficient * mass[i][j] * solution[e[i]] * solution[e[j]]
                end
            end
        end
    end
    f.boundary_energy, f.voltages = boundary_energy, voltages
    local energy = 0
    for _, t in ipairs(f.triangles) do
        local ex, ey = 0, 0
        if not t.electrode then
            for j, id in ipairs(t.ids) do
                ex = ex - t.gx[j] * solution[id] / m.unit
                ey = ey - t.gy[j] * solution[id] / m.unit
            end
        end
        check(
            require("luafemm-linear").finite(ex) and require("luafemm-linear").finite(ey),
            "field evaluation overflow"
        )
        t.Ex, t.Ey = ex, ey
        t.E = sqrt(ex * ex + ey * ey)
        t.Dx = E.epsilon0 * (t.tensor.xx * ex + t.tensor.xy * ey)
        t.Dy = E.epsilon0 * (t.tensor.xy * ex + t.tensor.yy * ey)
        t.D = sqrt(t.Dx * t.Dx + t.Dy * t.Dy)
        energy = energy
            + 0.5 * (ex * t.Dx + ey * t.Dy) * t.area * m.unit ^ 2 * metric.triangle(m, t)
    end
    check(require("luafemm-linear").finite(energy), "field energy overflow")
    f.solution, f.charges, f.energy, f.source_charge = solution, charges, energy, source
    f.balance = source + all_reactions + periodic_reaction - exterior_flux
    f.iterations, f.linear_stats = iterations, stats
    m.frozen = true
    return m
end
-- Spatial acceleration is ephemeral and never copied into an owned snapshot.
local indices = setmetatable({}, { __mode = "k" })
function F.index(m)
    local f = m.fe
    local index = indices[m]
    if not index then
        index = require("luafemm-spatial").new(f.options.h)
        for i, t in ipairs(f.triangles) do
            local a, b, c = f.nodes[t.ids[1]], f.nodes[t.ids[2]], f.nodes[t.ids[3]]
            index.insert(
                i,
                min(a[1], b[1], c[1]),
                min(a[2], b[2], c[2]),
                max(a[1], b[1], c[1]),
                max(a[2], b[2], c[2])
            )
        end
        indices[m] = index
    end
    return index
end
function F.barycentric(m, t, x, y)
    local a = m.fe.nodes[t.ids[1]]
    return {
        1 + t.gx[1] * (x - a[1]) + t.gy[1] * (y - a[2]),
        t.gx[2] * (x - a[1]) + t.gy[2] * (y - a[2]),
        t.gx[3] * (x - a[1]) + t.gy[3] * (y - a[2]),
    }
end
function F.locate(m, x, y)
    local f = m.fe
    local b = f.bounds
    if x < b[1] or x > b[3] or y < b[2] or y > b[4] then
        return nil
    end
    for _, i in ipairs(F.index(m).query(x, y, x, y)) do
        local t = f.triangles[i]
        local w = F.barycentric(m, t, x, y)
        if min(w[1], w[2], w[3]) >= -1e-10 then
            return i, w
        end
    end
end
--- At a mesh face the lowest triangle index is the deterministic raw side.
-- An explicit element index supports one-sided profile limits without averaging.
function F.sample(m, x, y, element)
    if m.hybrid then
        return require("luafemm-electric-hybrid").sample(m, x, y, element)
    end
    F.solve(m)
    local i, w
    if type(element) == "number" then
        i = element
        w = F.barycentric(m, m.fe.triangles[i], x, y)
    else
        i, w = F.locate(m, x, y)
    end
    if not i then
        return nil, "outside electric FE domain"
    end
    for id, charge in pairs(m.fe.point_loads or {}) do
        local p = m.fe.nodes[id]
        if charge ~= 0 and x == p[1] and y == p[2] then
            return nil, "native line/ring source singularity"
        end
    end
    local t = m.fe.triangles[i]
    local v = 0
    for j, id in ipairs(t.ids) do
        v = v + w[j] * m.fe.solution[id]
    end
    return {
        V = v,
        Ex = t.Ex,
        Ey = t.Ey,
        Er = m.fe.axisymmetric and t.Ex or nil,
        Ez = m.fe.axisymmetric and t.Ey or nil,
        Dr = m.fe.axisymmetric and t.Dx or nil,
        Dz = m.fe.axisymmetric and t.Dy or nil,
        E = t.E,
        Dx = t.Dx,
        Dy = t.Dy,
        D = t.D,
        element = i,
        region = t.region,
        electrode = t.electrode,
    }
end
function F.properties(m, name)
    if m.hybrid then
        return require("luafemm-electric-hybrid").properties(m, name)
    end
    F.solve(m)
    local f = m.fe
    if name and name ~= "" then
        check(f.charges[name] ~= nil, "unknown electrode: " .. name)
        local v
        if name == "outer boundary" then
            local multiple = false
            for id, owner in pairs(f.owners) do
                if owner == name then
                    if v and v ~= f.fixed[id] then
                        multiple = true
                    end
                    v = f.fixed[id]
                end
            end
            if multiple then
                v = nil
            end
        else
            v = f.voltages[name]
        end
        return {
            ["charge per length"] = not f.axisymmetric and f.charges[name] or nil,
            charge = f.charges[name] * f.depth,
            voltage = v,
        }
    end
    return {
        nodes = #f.nodes,
        triangles = #f.triangles,
        iterations = f.iterations,
        ["relative residual"] = not f.linear_stats.imported and f.linear_stats.relative_residual
            or nil,
        ["energy per length"] = not f.axisymmetric and f.energy or nil,
        energy = f.energy * f.depth,
        ["periodic flux"] = f.periodic_flux * f.depth,
        ["boundary energy"] = f.boundary_energy * f.depth,
        ["total energy"] = (f.energy + f.boundary_energy) * f.depth,
        ["charge per length"] = not f.axisymmetric and f.source_charge or nil,
        charge = f.source_charge * f.depth,
        ["charge balance per length"] = not f.axisymmetric and f.balance or nil,
        ["charge balance"] = f.balance * f.depth,
    }
end
function F.mesh_paths(m)
    local f = F.mesh(m)
    local out = {}
    for _, e in ipairs(f.edges) do
        out[#out + 1] = { points = { f.nodes[e[1]], f.nodes[e[2]] } }
    end
    return out
end
--- Restore the flattened physical drawing from an owned model.
function F.geometry(m)
    local f = m.fe
    local b = f.bounds
    if f.native then
        local out = {}
        for _, s in ipairs(f.native.lines) do
            if not s.hidden then
                out[#out + 1] = {
                    points = { s.a, s.b },
                    style = s.conductor > 0 and "femm electric electrode"
                        or "femm electric boundary",
                }
            end
        end
        return out
    end
    local out = {
        {
            points = { { b[1], b[2] }, { b[3], b[2] }, { b[3], b[4] }, { b[1], b[4] } },
            closed = true,
            style = "femm electric boundary",
        },
    }
    for _, r in ipairs(f.regions) do
        -- Separate contours retain boundaries; compound fill stays in the user's
        -- original TikZ drawing. The snapshot view deliberately draws no solid fill.
        for _, c in ipairs(r.contours) do
            out[#out + 1] =
                { points = c, closed = true, style = "femm electric dielectric,fill=none" }
        end
    end
    for _, p in ipairs(f.paths) do
        out[#out + 1] = {
            points = p.points,
            closed = p.closed,
            style = p.kind == "electrode"
                    and (p.solid and "femm electric conductor" or "femm electric electrode")
                or (p.kind == "charge" and "femm electric charge layer" or "femm electric boundary"),
        }
    end
    return out
end
return F
