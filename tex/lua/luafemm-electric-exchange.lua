-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Bounded, explicit native electrostatic exchange with polygonal geometry.
-- Imported native lengths are converted to metres. Export partitions the active
-- P1 mesh into labelled triangular cells: arbitrary PGF paths remain available,
-- while native consumers receive straight constraints and diagonal dielectrics.
-- The original TikZ drawing is never replaced by this derived representation.
-- @module luafemm-electric-exchange
local X = {}
local E = require("luafemm-electrostatic")
local F = require("luafemm-electric-fe")
local I = require("luafemm-fee")
local N = require("luafemm-electric-native")
local copy = require("luafemm-physical").copy
local function check(v, s)
    assert(v, "luafemm electric exchange: " .. s)
end
local function ref(v, list, name, zero)
    return E.integer(v, name, zero and 0 or 1, #list)
end
local function number(t, k, default)
    t[k] = E.number(t[k] or default, k)
    return t[k]
end
local function angle(x, y)
    return math.atan(y, x)
end
--- Validate native data before assigning any model state. Unknown physical
-- extensions, external-region transforms and contradictory ownership fail.
function X.import(d, o)
    d = copy(d)
    o = o or {}
    local h = d.header
    d.edition = o.edition or "femm5"
    check(d.edition == "femm4" or d.edition == "femm5", "unknown native edition")
    local unit = require("luafemm-fem").units[h.lengthunits]
    check(unit, "unsupported native length unit")
    check(
        (h.extro or 0) == 0 and (h.extri or 0) == 0 and (h.extzo or 0) == 0,
        "external-region transformation is unsupported"
    )
    for _, p in ipairs(d.pointprops) do
        number(p, "qp", 0)
        number(p, "vp", 0)
    end
    for _, p in ipairs(d.boundaries) do
        p.bdrytype = E.integer(p.bdrytype or 0, "boundary type", 0, 4)
        for _, k in ipairs({ "vs", "qs", "c0", "c1" }) do
            number(p, k, 0)
        end
        check(p.c0 >= 0 or p.bdrytype ~= 1, "negative Robin coefficient is unsupported")
        check(
            not d.solution or p.bdrytype < 3,
            "periodic native solution import is unsupported; import the .fee problem"
        )
    end
    for _, p in ipairs(d.materials) do
        check(number(p, "ex", 1) > 0 and number(p, "ey", 1) > 0, "permittivity must be positive")
        number(p, "qv")
    end
    local names = {}
    for _, p in ipairs(d.conductors) do
        check(
            type(p.conductorname) == "string"
                and p.conductorname:match("%S")
                and p.conductorname ~= "outer boundary",
            "invalid conductor name"
        )
        check(not names[p.conductorname], "duplicate conductor name")
        names[p.conductorname] = true
        number(p, "vc", 0)
        number(p, "qc", 0)
        p.conductortype = E.integer(p.conductortype or 0, "conductor type", 0, 1)
    end
    for _, p in ipairs(d.points) do
        p[1], p[2] = p[1] * unit, p[2] * unit
        if p[6] then
            p[6] = p[6] * unit
        end
        ref(p[3], d.pointprops, "point property", true)
        ref(p[5], d.conductors, "point conductor", true)
    end
    for _, p in ipairs(d.labels) do
        p[1], p[2], p[4] = p[1] * unit, p[2] * unit, p[4] * unit
        ref(p[3], d.materials, "dielectric label")
        check(p[6] == 0 or p[6] == 2, "external material labels are unsupported")
    end
    for _, p in ipairs(d.holes) do
        p[1], p[2] = p[1] * unit, p[2] * unit
    end
    if d.solution then
        for _, p in ipairs(d.solution.nodes) do
            p[1], p[2] = p[1] * unit, p[2] * unit
            E.integer(p[4], "solution node marker", -2, #d.conductors - 1)
        end
    end
    local lines = {}
    d.lines = lines
    local function segment(a, b, boundary, conductor, size, hidden)
        ref(boundary, d.boundaries, "segment boundary", true)
        ref(conductor, d.conductors, "segment conductor", true)
        check(a[1] ~= b[1] or a[2] ~= b[2], "zero-length native segment")
        check(not (boundary > 0 and conductor > 0), "segment combines boundary and conductor")
        lines[#lines + 1] = {
            a = { a[1], a[2] },
            b = { b[1], b[2] },
            boundary = boundary,
            conductor = conductor,
            h = size,
            hidden = hidden == 1,
        }
    end
    for _, s in ipairs(d.segments) do
        local a = d.points[E.integer(s[1], "segment vertex", 0, #d.points - 1) + 1]
        local b = d.points[E.integer(s[2], "segment vertex", 0, #d.points - 1) + 1]
        segment(a, b, s[4], s[7], s[3] > 0 and s[3] * unit or 0, s[5])
    end
    local tolerance = E.positive(o.curve_tolerance or 1e-5, "native curve tolerance (m)")
    for _, s in ipairs(d.arcs) do
        local a = d.points[E.integer(s[1], "arc vertex", 0, #d.points - 1) + 1]
        local b = d.points[E.integer(s[2], "arc vertex", 0, #d.points - 1) + 1]
        local sweep = E.positive(s[3], "arc angle") * math.pi / 180
        check(sweep < 2 * math.pi, "arc must be less than 360 degrees")
        local dx, dy = b[1] - a[1], b[2] - a[2]
        local chord = math.sqrt(dx * dx + dy * dy)
        check(chord > 0, "zero arc chord")
        local cx, cy =
            (a[1] + b[1]) / 2 - dy / (2 * math.tan(sweep / 2)),
            (a[2] + b[2]) / 2 + dx / (2 * math.tan(sweep / 2))
        local radius = chord / (2 * math.sin(sweep / 2))
        local start = angle(a[1] - cx, a[2] - cy)
        local step = 2 * math.acos(math.max(-1, 1 - tolerance / radius))
        step =
            math.min(step, E.positive(s[4], "arc subdivision angle") * math.pi / 180, math.pi / 4)
        local n = math.ceil(sweep / step)
        check(n <= 100000, "arc tolerance requires too many chords")
        local points = {}
        if d.solution then
            for _, p in ipairs(d.solution.nodes) do
                if
                    math.abs(math.sqrt((p[1] - cx) ^ 2 + (p[2] - cy) ^ 2) - radius)
                    <= 1e-8 * radius
                then
                    local t = (angle(p[1] - cx, p[2] - cy) - start) % (2 * math.pi)
                    if t > 2 * math.pi - 1e-8 then
                        t = 0
                    end
                    if t <= sweep + 1e-8 then
                        points[#points + 1] = { t, p }
                    end
                end
            end
            table.sort(points, function(u, v)
                return u[1] < v[1]
            end)
            check(#points >= 2, "solution arc lacks mesh vertices")
        else
            for j = 0, n do
                local t = sweep * j / n
                points[#points + 1] =
                    { t, { cx + radius * math.cos(start + t), cy + radius * math.sin(start + t) } }
            end
            points[1][2], points[#points][2] = a, b
        end
        for j = 2, #points do
            segment(points[j - 1][2], points[j][2], s[5], s[8], 0, s[6])
        end
    end
    local bounds = { math.huge, math.huge, -math.huge, -math.huge }
    local function include(p)
        bounds[1] = math.min(bounds[1], p[1])
        bounds[2] = math.min(bounds[2], p[2])
        bounds[3] = math.max(bounds[3], p[1])
        bounds[4] = math.max(bounds[4], p[2])
    end
    for _, p in ipairs(d.points) do
        include(p)
    end
    for _, s in ipairs(lines) do
        include(s.a)
        include(s.b)
    end
    if d.solution then
        for _, p in ipairs(d.solution.nodes) do
            include(p)
        end
    end
    local span = math.max(bounds[3] - bounds[1], bounds[4] - bounds[2])
    check(span > 0, "empty native geometry")
    if h.problemtype == "axisymmetric" then
        check(bounds[1] >= -span * 1e-12, "negative radius")
        bounds[1] = math.max(0, bounds[1])
    end
    local size = E.positive(o.mesh_size or span / 16, "native mesh size (m)")
    if not o.mesh_size then
        for _, p in ipairs(d.points) do
            if p[6] and p[6] > 0 then
                size = math.min(size, p[6])
            end
        end
        for _, p in ipairs(d.labels) do
            if p[4] > 0 then
                size = math.min(size, p[4])
            end
        end
    end
    local units = { m = 1, cm = 0.01, mm = 0.001, um = 1e-6 }
    local display = o.length_unit or "m"
    local scale = units[display]
    check(scale, "unknown display length unit")
    local function rescale(p)
        p[1], p[2] = p[1] / scale, p[2] / scale
    end
    for _, p in ipairs(d.points) do
        rescale(p)
    end
    for _, p in ipairs(d.labels) do
        rescale(p)
        p[4] = p[4] / scale
    end
    for _, p in ipairs(d.holes) do
        rescale(p)
    end
    for _, p in ipairs(lines) do
        rescale(p.a)
        rescale(p.b)
        p.h = p.h / scale
    end
    if d.solution then
        for _, p in ipairs(d.solution.nodes) do
            rescale(p)
        end
    end
    for i = 1, 4 do
        bounds[i] = bounds[i] / scale
    end
    local m = E.new({
        model = "finite elements",
        length_unit = display,
        problem_type = h.problemtype,
        xmin = bounds[1],
        ymin = bounds[2],
        xmax = bounds[3],
        ymax = bounds[4],
        depth = h.problemtype == "planar"
                and E.positive(h.depth or 1, "native depth") * unit / scale
            or nil,
        mesh_size = size / scale,
        outer_boundary = "insulated",
        boundary_tolerance = 0,
        solve_tolerance = o.solve_tolerance or h.precision or 1e-10,
    })
    m.fe.native = d
    m.descriptor.native_format = "FEMM electrostatic format 1"
    m.descriptor.validity = "native labelled dielectric cells; raw P1 gradients"
    for _, p in ipairs(d.conductors) do
        local c = {
            name = p.conductorname,
            condition = p.conductortype == 1 and "voltage" or "charge",
            potential = p.conductortype == 1 and p.vc or nil,
            charge = p.conductortype == 0 and p.qc or nil,
        }
        m.fe.electrodes[c.name] = c
        table.insert(m.fe.electrodes, c.name)
    end
    local periodic = {}
    for _, s in ipairs(lines) do
        if m.fe.axisymmetric and s.a[1] == 0 and s.b[1] == 0 then
            -- Native segment voltage/conductor data on the axis are nodal
            -- constraints, despite the zero surface measure. Bind after meshing.
            if s.boundary > 0 then
                check(
                    d.boundaries[s.boundary].bdrytype == 0,
                    "only a native voltage condition may constrain the regular axis"
                )
            end
        elseif s.conductor > 0 then
            local c = m.fe.electrodes[d.conductors[s.conductor].conductorname]
            F.path(m, { s.a, s.b }, false, "electrode", c)
        elseif s.boundary > 0 then
            local p = d.boundaries[s.boundary]
            if p.bdrytype == 0 then
                -- A native prescribed-potential segment may be internal. Represent it as
                -- a named, zero-thickness electrode rather than an exterior condition.
                local name = "native voltage " .. s.boundary
                check(not names[name], "reserved native voltage name collision")
                F.path(m, { s.a, s.b }, false, "electrode", { name = name, potential = p.vs })
            elseif p.bdrytype == 1 then
                F.path(
                    m,
                    { s.a, s.b },
                    false,
                    "boundary",
                    { condition = "mixed", coefficient = p.c0, displacement = p.c1 }
                )
            elseif p.bdrytype == 2 then
                F.path(
                    m,
                    { s.a, s.b },
                    false,
                    "boundary",
                    { condition = "displacement", displacement = -p.qs }
                )
                m.fe.paths[#m.fe.paths].native_surface = p.qs * (d.edition == "femm5" and 2 or 1)
            else
                periodic[s.boundary] = periodic[s.boundary] or {}
                table.insert(periodic[s.boundary], s)
            end
        end
    end
    for id, list in pairs(periodic) do
        check(#list == 2, "native periodic property requires exactly two straight segments")
        local a, b = list[1], list[2]
        if (a.b[1] - a.a[1]) * (b.b[1] - b.a[1]) + (a.b[2] - a.a[2]) * (b.b[2] - b.a[2]) < 0 then
            b.a, b.b = b.b, b.a
        end
        for _, s in ipairs(list) do
            F.path(m, { s.a, s.b }, false, "boundary", {
                condition = d.boundaries[id].bdrytype == 3 and "periodic" or "antiperiodic",
                pair = "native pair " .. id,
            })
        end
    end
    if d.solution then
        m.fe.imported_solution = {}
        for i, p in ipairs(d.solution.nodes) do
            m.fe.imported_solution[i] = p[3]
        end
    end
    m.fe.native_locked = true
    return m
end
function X.read(path, o)
    return X.import(I.read(path), o)
end

--- Export each triangle as a native material cell. This deliberately favours
-- exact region ownership over compact CAD geometry. Solid interiors become holes.
-- Off-diagonal tensors and periodic answer meshes are rejected, never flattened
-- into an isotropic surrogate. All validation precedes opening the output file.
function X.document(model, answer, edition)
    check(not model.hybrid, "native exports are unavailable with source correction")
    edition = edition or "femm5"
    check(edition == "femm4" or edition == "femm5", "unknown native edition")
    check(model.fe, "native exchange requires a finite-element electrostatic model")
    local m = copy(model)
    local f = F.mesh(m)
    check(#f.periodic_pairs == 0, "periodic native export is unsupported")
    if answer then
        E.freeze(m)
    end
    local d = {
        header = {
            format = 1,
            precision = f.tolerance,
            minangle = 0,
            depth = f.axisymmetric and 1 or f.depth,
            lengthunits = "meters",
            problemtype = f.axisymmetric and "axisymmetric" or "planar",
            coordinates = "cartesian",
            comment = "luafemm polygonal finite-element exchange; " .. edition,
        },
        pointprops = {},
        boundaries = {},
        materials = {},
        conductors = {},
        points = {},
        segments = {},
        arcs = {},
        holes = {},
        labels = {},
    }
    local conductors = {}
    for _, name in ipairs(f.electrodes) do
        local c = f.electrodes[name]
        conductors[name] = #d.conductors + 1
        d.conductors[#d.conductors + 1] = {
            conductorname = name,
            conductortype = c.condition == "voltage" and 1 or 0,
            vc = c.potential or 0,
            qc = c.charge or 0,
        }
    end
    local active = {}
    for _, t in ipairs(f.triangles) do
        if not t.electrode then
            for _, id in ipairs(t.ids) do
                active[id] = true
            end
        end
    end
    for i, p in ipairs(f.nodes) do
        local property = 0
        local owner = f.owners[i]
        local conductor = active[i] and owner and conductors[owner] or 0
        local load = (f.point_loads or {})[i]
        if active[i] and ((f.fixed[i] ~= nil and conductor == 0) or load) then
            local qp = load
                    and load ~= 0
                    and load / (f.axisymmetric and 2 * math.pi * p[1] * m.unit or 1)
                or 0
            property = #d.pointprops + 1
            d.pointprops[property] = { pointname = "node " .. i, vp = f.fixed[i] or 0, qp = qp }
        end
        d.points[i] = { p[1] * m.unit, p[2] * m.unit, property, 0, conductor }
    end
    local element_labels = {}
    for i, t in ipairs(f.triangles) do
        local x, y = 0, 0
        for _, id in ipairs(t.ids) do
            x = x + f.nodes[id][1] * m.unit / 3
            y = y + f.nodes[id][2] * m.unit / 3
        end
        if t.electrode then
            d.holes[#d.holes + 1] = { x, y }
        else
            check(
                math.abs(t.tensor.xy) <= 1e-12 * math.max(t.tensor.xx, t.tensor.yy),
                "off-diagonal permittivity cannot be represented in FEMM"
            )
            local id = #d.materials + 1
            d.materials[id] =
                { blockname = "cell " .. i, ex = t.tensor.xx, ey = t.tensor.yy, qv = t.density }
            d.labels[#d.labels + 1] = { x, y, id, 0, 0, 0 }
            element_labels[i] = #d.labels - 1
        end
    end
    for _, e in ipairs(f.edges) do
        local boundary = 0

        if
            not e.electrode
            and (
                e.condition == "mixed"
                or e.condition == "displacement"
                or e.condition == "voltage"
            )
        then
            boundary = #d.boundaries + 1
            d.boundaries[boundary] = {
                bdryname = "edge " .. boundary,
                bdrytype = e.condition == "voltage" and 0 or 1,
                vs = e.potential or 0,
                qs = 0,
                c0 = e.coefficient or 0,
                c1 = e.displacement or 0,
            }
        end
        if e.loads ~= 0 then
            check(boundary == 0, "sheet overlaps an exterior boundary")
            boundary = #d.boundaries + 1
            d.boundaries[boundary] = {
                bdryname = "sheet " .. boundary,
                bdrytype = 2,
                vs = 0,
                qs = e.loads / (edition == "femm5" and 2 or 1),
                c0 = 0,
                c1 = 0,
            }
        end
        local visible = e.outer or e.electrode or boundary > 0
        if #e.triangles == 2 then
            local a, b = f.triangles[e.triangles[1]], f.triangles[e.triangles[2]]
            visible = visible
                or a.region ~= b.region
                or a.tensor.xx ~= b.tensor.xx
                or a.tensor.yy ~= b.tensor.yy
                or a.density ~= b.density
        end
        d.segments[#d.segments + 1] = {
            e[1] - 1,
            e[2] - 1,
            -1,
            boundary,
            visible and 0 or 1,
            0,
            e.electrode and conductors[e.electrode] or 0,
        }
    end
    if answer then
        local s = { nodes = {}, triangles = {}, conductors = {} }
        d.solution = s
        local used, map = {}, {}
        for i, t in ipairs(f.triangles) do
            if element_labels[i] then
                for _, id in ipairs(t.ids) do
                    used[id] = true
                end
            end
        end
        for i, p in ipairs(f.nodes) do
            if used[i] then
                local owner = f.owners[i]
                local marker = owner and conductors[owner] and conductors[owner] - 1
                    or (f.fixed[i] ~= nil and -1 or -2)
                map[i] = #s.nodes
                s.nodes[#s.nodes + 1] = { p[1] * m.unit, p[2] * m.unit, f.solution[i], marker }
            end
        end
        for i, t in ipairs(f.triangles) do
            if element_labels[i] then
                s.triangles[#s.triangles + 1] =
                    { map[t.ids[1]], map[t.ids[2]], map[t.ids[3]], element_labels[i] }
            end
        end
        for _, name in ipairs(f.electrodes) do
            s.conductors[#s.conductors + 1] = { f.voltages[name], f.charges[name] * f.depth }
        end
    end
    return d
end
function X.write(m, path, format, edition)
    check(format == "fee" or format == "res", "format must be fee or res")
    check(
        type(path) == "string" and path:lower():sub(-4) == "." .. format,
        "output filename must end in ." .. format
    )
    local text = I.serialize(X.document(m, format == "res", edition))
    require("luafemm-fem").write_text(path, text, m.fe.native and m.fe.native.source)
end
return X
