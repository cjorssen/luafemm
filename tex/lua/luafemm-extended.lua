-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Versioned, passive native problem/answer documents, independent of caches.
-- Adapters own evaluated geometry, constitutive data and ordered coefficients.
-- No imported string is executed or interpreted as a dependency filename.
-- Rebuilding element gradients/reactions is permitted; reading never solves.
-- Schema/provider revisions, not the release string, govern compatibility.
-- @module luafemm-extended
local X = {}
local C = require("luafemm-cache")
local copy = require("luafemm-physical").copy
local function check(v, s)
    assert(v, "luafemm extended: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
-- Sparse sequences and hidden dictionary entries must not escape ipairs.
local function sequence(t, name)
    check(type(t) == "table", "expected " .. name .. " sequence")
    local n, count = #t, 0
    for k in pairs(t) do
        check(
            type(k) == "number" and k % 1 == 0 and k >= 1 and k <= n,
            "invalid " .. name .. " sequence index"
        )
        count = count + 1
    end
    check(count == n, "sparse " .. name .. " sequence")
    return n
end
local function fields(t, names)
    local out = {}
    for name in names:gmatch("%S+") do
        out[name] = copy(t[name])
    end
    return out
end
local electric_fields = "descriptor sources material unit epsilon tensor bounds classic"
local fe_fields =
    "axisymmetric regions paths electrodes bounds density boundary_tolerance outer_boundary outer_potential depth options tolerance max_iterations"
local hybrid_fields =
    "quadrature_tolerance quadrature_depth quadrature_cells epsilon_r reference_radius"
local magnetic_fields =
    "materials regions boundaries frontiers domain pairs topology native_topology circuits point_constraints"
local adapters = {
    ["electric-analytic-1"] = true,
    ["electric-p1-1"] = true,
    ["electric-source-p1-1"] = true,
    ["magnetic-p1-1"] = true,
    ["magnetic-biot-1"] = true,
}
local function tensor(t)
    check(
        type(t) == "table"
            and finite(t.xx)
            and finite(t.xy)
            and finite(t.yy)
            and t.xx > 0
            and t.yy > 0
            and t.xx * t.yy > t.xy * t.xy,
        "invalid dielectric tensor"
    )
end
-- Pass declarations through the same validators as a TikZ builder. Catalogue
-- identifiers are provenance only: embedded values remain authoritative.
local function electric_model(p)
    local E, F = require("luafemm-electrostatic"), require("luafemm-electric-fe")
    sequence(p.sources, "sources")
    tensor(p.tensor)
    local f, h, c = p.fe, p.hybrid, p.classic
    local o = {
        length_unit = p.descriptor.length_unit,
        epsilon_r = p.epsilon / E.epsilon0,
        xmin = p.bounds[1],
        ymin = p.bounds[2],
        xmax = p.bounds[3],
        ymax = p.bounds[4],
    }
    if f then
        check(not c, "FE and analytical symmetry providers are mutually exclusive")
        o.model = h and "source correction" or "finite elements"
        o.problem_type = f.axisymmetric and "axisymmetric" or "planar"
        o.depth = not f.axisymmetric and f.depth / p.unit or nil
        o.volume_density, o.outer_boundary, o.outer_potential =
            f.density, f.outer_boundary, f.outer_potential
        o.boundary_tolerance, o.solve_tolerance, o.max_iterations =
            f.boundary_tolerance, f.tolerance, f.max_iterations
        o.mesh_size, o.mesh_tolerance, o.min_angle =
            f.options.h, f.options.mesh_tolerance, f.options.min_angle
        o.max_area, o.max_steiner = f.options.max_area, f.options.max_steiner
        if h then
            o.reference_radius = not f.axisymmetric and h.reference_radius or nil
            o.source_quadrature_tolerance = h.quadrature_tolerance
            o.source_quadrature_depth, o.source_quadrature_cells =
                h.quadrature_depth, h.quadrature_cells
        end
    elseif c then
        o.model, o.center_x, o.center_y = c.kind, c.x, c.y
        o.reference_potential, o.interface_value = c.reference, c.surface_value
        o.voltage, o.separation, o.normal_angle = c.voltage, c.separation, c.angle
        o.inner_radius, o.outer_radius, o.radius = c.inner, c.outer, c.radius
        o.reference_radius, o.distribution = c.reference_radius, c.distribution
        if c.distribution == "volume" then
            o.volume_density = c.density
        else
            o.surface_density = c.density
        end
    end
    local m = E.new(o)
    check(m.unit == p.unit, "length unit and scale disagree")
    for _, key in ipairs({ "geometry", "backend", "potential" }) do
        check(m.descriptor[key] == p.descriptor[key], "electric descriptor/provider mismatch")
    end
    m.descriptor, m.tensor, m.material = copy(p.descriptor), copy(p.tensor), p.material
    if f then
        sequence(f.regions, "dielectric regions")
        sequence(f.paths, "electric paths")
        for _, r in ipairs(f.regions) do
            tensor(r.tensor)
            F.region(m, r.contours, {
                fill_rule = r.fill_rule,
                mesh_size = r.mesh_size,
                epsilon_r = r.epsilon / E.epsilon0,
                density = r.density,
            })
            local q = m.fe.regions[#m.fe.regions]
            q.tensor, q.material = copy(r.tensor), r.material
        end
        for _, path in ipairs(f.paths) do
            local q = copy(path)
            q.interior = q.solid and "conductor" or "dielectric"
            F.path(m, path.points, path.closed, path.kind, q)
        end
    end
    check(not (c or f and not h) or #p.sources == 0, "this provider does not accept point sources")
    for _, s in ipairs(p.sources) do
        if h then
            require("luafemm-electric-hybrid").add(m, s.x, s.y, {
                kind = s.kind,
                charge = f.axisymmetric and s.q or nil,
                charge_per_length = not f.axisymmetric and s.q or nil,
                guard_radius = s.guard,
                radius = s.radius,
            })
        else
            E.add(m, s.x, s.y, s.q, s.radius)
        end
    end
    return m
end
local function magnetic_model(p)
    local M = require("luafemm")
    check(
        not p.options.source and not p.options.boundary and not p.options.log,
        "callbacks cannot be stored"
    )
    check(
        p.options.field_model == "open" and p.options.problem_type ~= "axial",
        "unsupported magnetic provider"
    )
    check(
        p.options.material_edition == "femm5" and p.options.material_file == "",
        "native materials must be embedded, without an external catalogue selector"
    )
    check(
        p.options.cache and p.options.cache.mode == "off",
        "native documents cannot select an execution cache"
    )
    local m = M.new(p.options)
    for name, prepared in pairs(p.materials) do
        local raw = copy(prepared)
        check(
            raw.library == nil or type(raw.library) == "table",
            "material provenance must be embedded, not a catalogue lookup"
        )
        raw.mur_x, raw.mur_y = raw.intrinsic_mur_x, raw.intrinsic_mur_y
        raw.curve, raw.tensor = nil, nil
        M.material(m, name, raw)
        check(
            C.encode(m.materials[name].tensor) == C.encode(prepared.tensor),
            "inconsistent prepared material tensor"
        )
    end
    for k, v in pairs(fields(p, magnetic_fields)) do
        if k ~= "materials" then
            m[k] = v
        end
    end
    for _, r in ipairs(m.regions) do
        check(m.materials[r.material], "unknown region material")
    end
    return m
end
local function meshcheck(mesh, values, imaginary)
    check(
        type(mesh) == "table" and type(mesh.nodes) == "table" and type(mesh.elements) == "table",
        "missing mesh"
    )
    local n = sequence(mesh.nodes, "nodes")
    sequence(mesh.elements, "elements")
    sequence(values, "coefficients")
    if imaginary then
        sequence(imaginary, "imaginary coefficients")
    end
    check(
        n >= 3 and n <= 1000000 and #mesh.elements > 0 and #mesh.elements <= 2000000,
        "mesh size outside supported range"
    )
    check(type(values) == "table" and #values == n, "coefficient count mismatch")
    check(not imaginary or #imaginary == n, "imaginary coefficient count mismatch")
    for i, p in ipairs(mesh.nodes) do
        sequence(p, "coordinate")
        check(
            type(p) == "table"
                and #p == 2
                and finite(p[1])
                and finite(p[2])
                and finite(values[i])
                and (not imaginary or finite(imaginary[i])),
            "invalid nodal data"
        )
    end
    for _, ids in ipairs(mesh.elements) do
        sequence(ids, "element")
        check(type(ids) == "table" and #ids == 3, "expected triangular P1 element")
        for _, id in ipairs(ids) do
            check(finite(id) and id % 1 == 0 and mesh.nodes[id], "invalid node index")
        end
        local a, b, c = mesh.nodes[ids[1]], mesh.nodes[ids[2]], mesh.nodes[ids[3]]
        check(
            (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]) > 0,
            "nonpositive element orientation"
        )
    end
    sequence(mesh.periodic_pairs or {}, "periodic pairs")
    for _, p in ipairs(mesh.periodic_pairs or {}) do
        sequence(p, "periodic pair")
        check(
            #p == 3 and mesh.nodes[p[1]] and mesh.nodes[p[2]] and (p[3] == 1 or p[3] == -1),
            "invalid periodic constraint"
        )
    end
end
local function meshrecord(nodes, triangles, stats, boundary, pairs)
    local mesh = {
        nodes = {},
        elements = {},
        stats = copy(stats),
        boundary = copy(boundary),
        periodic_pairs = copy(pairs or {}),
    }
    for i, p in ipairs(nodes) do
        mesh.nodes[i] = { p[1], p[2] }
    end
    for i, t in ipairs(triangles) do
        mesh.elements[i] = copy(t.ids)
    end
    return mesh
end
local function electric(m, answer)
    check(not m.fe or not m.fe.native, "native electric topology adapter is unavailable")
    local p = fields(m, electric_fields)
    p.fe = m.fe and fields(m.fe, fe_fields)
    p.hybrid = m.hybrid and fields(m.hybrid, hybrid_fields)
    local adapter = m.hybrid and "electric-source-p1-1"
        or m.fe and "electric-p1-1"
        or "electric-analytic-1"
    local s
    if answer then
        require("luafemm-electrostatic").freeze(m)
        s = {
            representation = m.hybrid and "analytic-plus-p1-correction"
                or m.fe and "scalar-p1"
                or "registered-analytic",
            frozen = true,
        }
        if m.fe then
            s.mesh = meshrecord(m.fe.nodes, m.fe.triangles, m.fe.stats, nil, m.fe.periodic_pairs)
            s.values = copy(m.hybrid and m.hybrid.solution or m.fe.solution)
            s.residual = m.hybrid and m.hybrid.residual or m.fe.linear_stats.relative_residual
        end
    end
    return adapter, p, s
end
local function magnetic(m, answer)
    check(
        m.options.field_model == "open"
            and m.options.problem_type ~= "axial"
            and not m.operating
            and not m.sliding
            and not m.transient
            and not m.import_answer
            and not m.phasor_model
            and not m.network
            and not require("luafemm-flux-drive").active(m),
        "magnetic adapter requires an ordinary static or harmonic planar/axisymmetric model"
    )
    check(not m.options.source and not m.options.boundary, "callbacks cannot be stored")
    -- Solving compiles face topology and circuit-derived sources. Capture the
    -- declarations after that preparation, so a first answer save is identical
    -- to saving an already drawn model. Problem-only capture still never solves.
    if answer then
        require("luafemm").solve(m)
    end
    local p = fields(m, magnetic_fields)
    p.options = copy(m.options)
    p.options.log, p.options.cache = nil, { mode = "off" }
    -- Material records already own their values and edition provenance. Do not
    -- reopen a user-supplied FEMM 4 catalogue through M.new during restoration.
    p.options.material_edition, p.options.material_file = "femm5", ""
    local s
    if answer then
        s = {
            representation = m.options.frequency > 0 and "magnetic-p1-phasor" or "magnetic-p1",
            mesh = meshrecord(m.nodes, m.triangles, m.mesh_stats, m.exterior_nodes),
            values = copy(m.A),
            imaginary = copy(m.A_im),
            stats = copy(m.stats),
            circuit_state = copy(m.circuit_state),
        }
    end
    return "magnetic-p1-1", p, s
end
--- Capture owned physical data. Answer capture may explicitly solve a builder.
-- kind is problem or answer; saving a problem never performs a calculation.
function X.capture(m, kind)
    kind = kind or "problem"
    check(kind == "problem" or kind == "answer", "unknown document kind")
    if m.biot then
        local r = {
            schema = 1,
            kind = kind,
            adapter = "magnetic-biot-1",
            required = { "magnetic-biot-1" },
            physics = "magnetism",
            restart = false,
            descriptor = copy(m.descriptor),
            problem = { options = copy(m.options), sources = copy(m.sources) },
            answer = kind == "answer" and { representation = "registered-biot-1" } or nil,
            provenance = { package = "luafemm", version = require("luafemm").version },
        }
        X.validate(r)
        return r
    end
    local is_electric = m.descriptor and m.descriptor.physics == "electrostatics"
    local adapter, p, s = (is_electric and electric or magnetic)(m, kind == "answer")
    local r = {
        schema = 1,
        kind = kind,
        adapter = adapter,
        required = { adapter },
        physics = is_electric and "electrostatics" or "magnetism",
        descriptor = require("luafemm-physical").describe(m),
        problem = p,
        answer = s,
        provenance = { package = "luafemm", version = require("luafemm").version },
        restart = false,
    }
    X.validate(r)
    return r
end
--- Validate identification, dimensions, representation and every mesh index.
-- The codec additionally bounds bytes, nesting and record counts and refuses
-- non-finite values, functions, duplicate keys and trailing payloads.
function X.validate(r)
    check(type(r) == "table" and r.schema == 1, "unknown schema version")
    check(r.kind == "problem" or r.kind == "answer", "unknown document kind")
    check(adapters[r.adapter], "unsupported adapter")
    sequence(r.required, "required capability")
    check(
        type(r.required) == "table" and #r.required == 1 and r.required[1] == r.adapter,
        "unknown required capability"
    )
    check(r.restart == false, "restart capability is unavailable; use .lft")
    check(type(r.problem) == "table" and type(r.descriptor) == "table", "missing physical data")
    local d, p = r.descriptor, r.problem
    if r.adapter == "magnetic-biot-1" then
        local m = require("luafemm-biot").restore(p)
        check(
            r.physics == "magnetism" and C.encode(d) == C.encode(m.descriptor),
            "prescribed-current descriptor mismatch"
        )
        if r.kind == "answer" then
            check(
                type(r.answer) == "table"
                    and r.answer.representation == "registered-biot-1"
                    and not r.answer.mesh
                    and not r.answer.values,
                "prescribed-current solution basis mismatch"
            )
            for key in pairs(r.answer) do
                check(key == "representation", "unknown prescribed-current answer field")
            end
        else
            check(not r.answer, "problem must not contain an answer")
        end
        C.encode(r)
        return r
    end
    check(
        d.physics == r.physics
            and d.physical_dimension == 3
            and d.evaluation_dimension == 2
            and finite(d.length_scale)
            and d.length_scale > 0,
        "invalid physical descriptor"
    )
    check((r.physics == "magnetism") == (r.adapter == "magnetic-p1-1"), "physics/adapter mismatch")
    check(r.physics == "magnetism" or r.physics == "electrostatics", "unknown physics")
    if r.physics == "electrostatics" then
        check(p.descriptor and C.encode(p.descriptor) == C.encode(d), "descriptor mismatch")
        check(
            finite(p.unit) and p.unit == d.length_scale and finite(p.epsilon) and p.epsilon > 0,
            "invalid electric units or permittivity"
        )
        check(
            (not not p.fe) == (r.adapter ~= "electric-analytic-1")
                and (not not p.hybrid) == (r.adapter == "electric-source-p1-1"),
            "electric provider mismatch"
        )
        electric_model(p)
    else
        check(p.options and p.options.unit == d.length_scale, "magnetic unit mismatch")
        check(
            d.geometry == p.options.problem_type
                and d.backend == "finite elements"
                and d.potential == "magnetic vector potential",
            "magnetic descriptor/provider mismatch"
        )
        magnetic_model(p)
    end
    if r.kind == "answer" then
        local a = r.answer
        check(type(a) == "table", "missing answer")
        if r.physics == "magnetism" then
            check(type(a.stats) == "table", "missing magnetic answer diagnostics")
        end
        local repr = r.adapter == "electric-analytic-1" and "registered-analytic"
            or r.adapter == "electric-source-p1-1" and "analytic-plus-p1-correction"
            or r.adapter == "electric-p1-1" and "scalar-p1"
            or p.options.frequency > 0 and "magnetic-p1-phasor"
            or "magnetic-p1"
        check(a.representation == repr, "solution basis mismatch")
        if repr == "registered-analytic" then
            check(not a.mesh and not a.values, "analytical answer cannot contain FE coefficients")
        else
            meshcheck(a.mesh, a.values, a.imaginary)
            check(
                (not not a.imaginary) == (repr == "magnetic-p1-phasor"),
                "phasor coefficient mismatch"
            )
        end
    else
        check(not r.answer, "problem must not contain an answer")
    end
    C.encode(r)
    return r
end
--- Restore independent data; never solve, execute input or access dependencies.
function X.restore(record)
    X.validate(record)
    local r = copy(record)
    local p, a = r.problem, r.answer
    if r.adapter == "magnetic-biot-1" then
        local m = require("luafemm-biot").restore(p)
        m.frozen = r.kind == "answer"
        m.extended_provenance = copy(r.provenance)
        return m
    end
    if r.physics == "electrostatics" then
        p = electric_model(p)
        if a and p.fe then
            p.fe.extended_mesh = a.mesh
            if p.hybrid then
                p.hybrid.imported_solution = a.values
                p.hybrid.imported_residual = a.residual
            else
                p.fe.imported_solution = a.values
            end
            require("luafemm-electric-fe").solve(p)
            local f = p.fe
            local values = {}
            for i, point in ipairs(f.nodes) do
                if f.floating[i] or f.periodic_nodes[i] then
                    values[i] = p.hybrid
                            and p.hybrid.solution[i] + require("luafemm-electric-hybrid").base(
                                p,
                                point[1],
                                point[2]
                            ).V
                        or f.solution[i]
                end
            end
            for _, pair in ipairs(f.periodic_pairs) do
                local u, v = values[pair[1]], pair[3] * values[pair[2]]
                check(
                    math.abs(u - v) <= 1e-8 * math.max(1, math.abs(u), math.abs(v)),
                    "restored field violates a periodic constraint"
                )
            end
            local floating = {}
            for id, name in pairs(f.floating) do
                local u, v = values[id], floating[name]
                check(
                    not v or math.abs(u - v) <= 1e-8 * math.max(1, math.abs(u), math.abs(v)),
                    "restored field violates an equipotential conductor"
                )
                floating[name] = u
            end
        elseif a then
            p.frozen = true
        end
        p.extended_provenance = copy(r.provenance)
        return p
    end
    local M = require("luafemm")
    local m = magnetic_model(p)
    if a then
        if m.topology then
            require("luafemm-static").excite(m, m.topology)
        end
        m.mesh_stats = a.mesh.stats
        M.prepare_mesh(m, a.mesh.nodes, a.mesh.elements, a.mesh.boundary)
        local real, imaginary = {}, {}
        for i, node in ipairs(m.nodes) do
            local function close(x, y)
                return math.abs(x - y) <= 1e-8 * math.max(1e-12, math.abs(x), math.abs(y))
            end
            if node.dof then
                local sign = node.dof_sign or 1
                local value = a.values[i] / sign
                check(
                    not real[node.dof] or close(real[node.dof], value),
                    "restored potential violates a periodic constraint"
                )
                real[node.dof] = value
                if a.imaginary then
                    value = a.imaginary[i] / sign
                    check(
                        not imaginary[node.dof] or close(imaginary[node.dof], value),
                        "restored phasor violates a periodic constraint"
                    )
                    imaginary[node.dof] = value
                end
            else
                check(close(m.A[i], a.values[i]), "restored potential violates a fixed constraint")
                check(
                    not a.imaginary or close(0, a.imaginary[i]),
                    "restored phasor violates a fixed constraint"
                )
            end
        end
        m.A, m.A_im, m.stats, m.circuit_state = a.values, a.imaginary, a.stats, a.circuit_state
        M.update_fields(m)
    end
    m.extended_provenance = copy(r.provenance)
    return m
end
function X.write(m, file, kind)
    local record = X.capture(m, kind)
    C.write_data(file, record, record.kind == "answer" and "LUAFEMM-ANSX" or "LUAFEMM-FEMX")
    return record
end
function X.read(file, kind)
    check(kind == "problem" or kind == "answer", "specify problem or answer")
    local r = C.read_data(file, kind == "answer" and "LUAFEMM-ANSX" or "LUAFEMM-FEMX")
    check(r.kind == kind, "document kind mismatch")
    return X.restore(r)
end
return X
