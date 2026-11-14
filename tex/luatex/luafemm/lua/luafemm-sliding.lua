-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Circular moving-band meshes for independent planar static rotor positions.
-- Bulk vertices/connectivity are owned and reused. Only source-free air-band
-- triangles are rebuilt; this is a P1 moving band, not FEMM's special AGE.
-- @module luafemm-sliding
local S = {}
local pred = require("luafemm-predicates")
local geometry = require("luafemm-geometry")
local function check(v, message)
    assert(v, "luafemm sliding band: " .. message)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
local function copy(v)
    if type(v) ~= "table" then
        return v
    end
    local out = {}
    for k, x in pairs(v) do
        out[k] = copy(x)
    end
    return out
end
local function edge(a, b)
    return math.min(a, b) .. ":" .. math.max(a, b)
end
local function radius(p, o)
    return math.sqrt((p[1] - o.center_x) ^ 2 + (p[2] - o.center_y) ^ 2)
end
local function rotate(p, o, angle)
    local c, s = math.cos(math.rad(angle)), math.sin(math.rad(angle))
    local x, y = p[1] - o.center_x, p[2] - o.center_y
    return { o.center_x + c * x - s * y, o.center_y + s * x + c * y }
end
local function circle(o, r)
    local points = {}
    for i = 0, o.segments - 1 do
        local a = 2 * math.pi * i / o.segments
        points[#points + 1] = { o.center_x + r * math.cos(a), o.center_y + r * math.sin(a) }
    end
    return points
end
local function eligible(m)
    check(m and not m.operating and not m.sliding, "use an ordinary native template")
    local o = m.options
    check(
        o.problem_type == "planar" and o.frequency == 0 and o.field_model == "open",
        "only planar static open-field models are supported"
    )
    check(o.mesher == "delaunay", "use mesher=delaunay")
    check(
        not m.import_document and not m.domain and not m.frontiers and not m.pairs,
        "imported/custom domains and tagged/periodic frontiers are unsupported"
    )
    check(not o.source and not o.boundary, "callbacks cannot be captured")
    for _, c in pairs(m.circuits or {}) do
        check(
            c.type == "series" and c.drive == "current" and c.phase == 0,
            "use prescribed DC series currents"
        )
    end
    for _, r in ipairs(m.regions) do
        check((r.current_phase or 0) == 0, "use DC sources")
    end
    local p = m.materials.air
    check(
        p and not p.bh and p.hc == 0 and p.tensor[1] == 1 and p.tensor[2] == 0 and p.tensor[3] == 1,
        "the air material must retain vacuum permeability and zero coercivity"
    )
end
local function ringsort(nodes, ids, o)
    table.sort(ids, function(a, b)
        local p, q = nodes[a], nodes[b]
        return math.atan(p[2] - o.center_y, p[1] - o.center_x)
            < math.atan(q[2] - o.center_y, q[1] - o.center_x)
    end)
end
-- Greedy visible bridges between strictly nested convex polygons. Choosing the
-- shorter next bridge avoids winding the strip around the hole. Positive exact
-- orientations and the subsequent incidence/area audit certify the result.
local function bridge(nodes, inner, outer, tri)
    local start, best = 1, math.huge
    local function length(a, b)
        return (nodes[a][1] - nodes[b][1]) ^ 2 + (nodes[a][2] - nodes[b][2]) ^ 2
    end
    for i, id in ipairs(inner) do
        local d = length(id, outer[1])
        if d < best then
            start, best = i, d
        end
    end
    local i, j = 0, 0
    local function ia(k)
        return inner[(start + k - 1) % #inner + 1]
    end
    local function ob(k)
        return outer[k % #outer + 1]
    end
    while i < #inner or j < #outer do
        local a, b, c, d = ia(i), ob(j), ia(i + 1), ob(j + 1)
        local left = i < #inner and pred.orient(nodes[a], nodes[b], nodes[c]) > 0
        local right = j < #outer and pred.orient(nodes[a], nodes[b], nodes[d]) > 0
        check(left or right, "cannot triangulate the band; increase segments or widen the band")
        if left and (not right or length(c, b) <= length(a, d)) then
            tri[#tri + 1] = { a, b, c }
            i = i + 1
        else
            tri[#tri + 1] = { a, b, d }
            j = j + 1
        end
    end
end
-- Recheck the full assembled disk, including every protected bulk edge. No
-- partial mesh or inverted band element can reach the field solver/cache.
local function audit(m, nodes, tri, boundary, stats)
    local edges, used, area, minangle, maxarea = {}, {}, 0, 180, 0
    for _, t in ipairs(tri) do
        local twice = pred.orient(nodes[t[1]], nodes[t[2]], nodes[t[3]])
        check(twice > 0, "nonpositive moving-band triangle")
        area = area + twice / 2
        maxarea = math.max(maxarea, twice / 2)
        for j = 1, 3 do
            local a, b, c = t[j], t[j % 3 + 1], t[(j + 1) % 3 + 1]
            local key = edge(a, b)
            local e = edges[key] or { a, b, count = 0 }
            edges[key] = e
            e.count = e.count + 1
            check(e.count <= 2, "non-manifold band mesh")
            used[a] = true
            local p, q, r = nodes[a], nodes[b], nodes[c]
            local dot = (q[1] - p[1]) * (r[1] - p[1]) + (q[2] - p[2]) * (r[2] - p[2])
            minangle = math.min(minangle, math.deg(math.atan(twice, dot)))
        end
    end
    local ne = 0
    for _, e in pairs(edges) do
        ne = ne + 1
        if e.count == 1 then
            check(boundary[e[1]] and boundary[e[2]], "unclosed air band")
        end
    end
    for i in ipairs(nodes) do
        check(used[i], "unused band vertex")
    end
    check(#nodes - ne + #tri == 1, "invalid band mesh topology")
    local o = m.options
    local expected = (o.xmax - o.xmin) * (o.ymax - o.ymin) * o.unit ^ 2
    check(math.abs(area / expected - 1) < 1e-9, "band mesh coverage mismatch")
    for _, e in ipairs(stats.segments) do
        check(edges[edge(e[1], e[2])], "lost bulk constraint")
    end
    stats.min_angle, stats.max_area = minangle, maxarea / o.unit ^ 2
    stats.boundary_nodes = copy(boundary)
    -- The initial mesh's quality targets do not certify a newly connected band.
    if o.min_angle > 0 then
        check(minangle + 1e-8 >= o.min_angle, "moving band misses minimum angle target")
    end
    if o.max_area > 0 then
        check(stats.max_area <= o.max_area * (1 + 1e-8), "moving band misses maximum area target")
    end
end
--- Mesh and capture a template after all declarations, before any field solve.
-- An air collar must separate every physical path from both circular traces.
-- The caller's model is meshed once; returned data are independently owned.
-- @tparam table m Unmeshed native planar DC model.
-- @tparam table o inner_radius, outer_radius; optional center_x/y (model units),
-- segments (default 96), layers (default 2).
-- @treturn table Owned template for repeated new(template,options) calls.
function S.capture(m, o)
    eligible(m)
    check(not m.nodes, "capture the sliding band before meshing or solving")
    o = copy(o or {})
    for k in pairs(o) do
        check(
            k == "inner_radius"
                or k == "outer_radius"
                or k == "center_x"
                or k == "center_y"
                or k == "segments"
                or k == "layers",
            "unknown option " .. tostring(k)
        )
    end
    o.center_x, o.center_y = o.center_x or 0, o.center_y or 0
    o.segments, o.layers = o.segments or 96, o.layers or 2
    for _, k in ipairs({
        "center_x",
        "center_y",
        "inner_radius",
        "outer_radius",
        "segments",
        "layers",
    }) do
        check(finite(o[k]), "invalid " .. k)
    end
    check(
        o.inner_radius > 0 and o.outer_radius > o.inner_radius,
        "require 0 < inner radius < outer radius"
    )
    check(
        o.segments >= 16 and o.segments <= 1024 and o.segments % 1 == 0,
        "segments must be an integer from 16 to 1024"
    )
    check(
        o.layers >= 1 and o.layers <= 16 and o.layers % 1 == 0,
        "layers must be an integer from 1 to 16"
    )
    local dr = (o.outer_radius - o.inner_radius) / o.layers
    check(
        (o.inner_radius + dr) * math.cos(math.pi / o.segments) > o.inner_radius,
        "polygonal layers overlap; increase segments or reduce layers"
    )
    local bounds = m.options
    check(
        o.center_x - o.outer_radius > bounds.xmin
            and o.center_x + o.outer_radius < bounds.xmax
            and o.center_y - o.outer_radius > bounds.ymin
            and o.center_y + o.outer_radius < bounds.ymax,
        "the complete band must lie strictly inside the problem box"
    )
    local rotating = {}
    local inner_limit = o.inner_radius * math.cos(math.pi / o.segments)
    for i, r in ipairs(m.regions) do
        local maximum, minimum = 0, math.huge
        for _, contour in ipairs(r.contours) do
            for j, p in ipairs(contour) do
                maximum = math.max(maximum, radius(p, o))
                local distance = require("luafemm-topology").distance(
                    { o.center_x, o.center_y },
                    p,
                    contour[j % #contour + 1]
                )
                minimum = math.min(minimum, distance)
            end
        end
        if maximum < inner_limit * (1 - 1e-8) then
            rotating[i] = true
        else
            check(
                minimum > o.outer_radius * (1 + 1e-8)
                    and not geometry.contains(
                        r,
                        o.center_x + (o.inner_radius + o.outer_radius) / 2,
                        o.center_y
                    ),
                "a region crosses or fills the air band; keep both traces inside source-free air"
            )
        end
    end
    local M = require("luafemm")
    M.region_contours(
        m,
        "air",
        { circle(o, o.outer_radius), circle(o, o.inner_radius) },
        0,
        0,
        0,
        "even odd"
    )
    local marker = #m.regions
    m.frontiers = {} -- Request a native arrangement even without named circuits.
    M.mesh(m)
    local bulk, used, edges = {}, {}, {}
    for _, e in ipairs(m.triangles) do
        if m.topology.faces[e.face].region ~= marker then
            bulk[#bulk + 1] = copy(e.ids)
            for j, id in ipairs(e.ids) do
                used[id] = true
                local b = e.ids[j % 3 + 1]
                local key = edge(id, b)
                local v = edges[key] or { id, b, count = 0 }
                edges[key] = v
                v.count = v.count + 1
            end
        end
    end
    local in_set, out_set = {}, {}
    for _, e in pairs(edges) do
        if e.count == 1 and not (m.exterior_nodes[e[1]] and m.exterior_nodes[e[2]]) then
            for _, id in ipairs({ e[1], e[2] }) do
                local p = { m.nodes[id][1] / bounds.unit, m.nodes[id][2] / bounds.unit }
                local set = radius(p, o) < (o.inner_radius + o.outer_radius) / 2 and in_set
                    or out_set
                set[id] = true
            end
        end
    end
    local nodes, map, inner, outer, boundary, moving = {}, {}, {}, {}, {}, {}
    for id, p in ipairs(m.nodes) do
        if used[id] then
            local n = #nodes + 1
            map[id] = n
            nodes[n] = { p[1], p[2] }
            if in_set[id] then
                inner[#inner + 1] = n
            end
            if out_set[id] then
                outer[#outer + 1] = n
            end
            if m.exterior_nodes[id] then
                boundary[n] = true
            end
            moving[n] = radius({ p[1] / bounds.unit, p[2] / bounds.unit }, o)
                    <= o.inner_radius * (1 + 1e-8)
                or nil
        end
    end
    check(#inner >= o.segments and #outer >= o.segments, "incomplete circular traces")
    for _, e in ipairs(bulk) do
        for j, id in ipairs(e) do
            e[j] = map[id]
        end
    end
    local stats = copy(m.mesh_stats)
    stats.segments = {}
    for _, e in ipairs(m.mesh_stats.segments) do
        if map[e[1]] and map[e[2]] and edges[edge(e[1], e[2])] then
            stats.segments[#stats.segments + 1] = { map[e[1]], map[e[2]] }
        end
    end
    local model = {}
    for _, k in ipairs({
        "options",
        "materials",
        "regions",
        "boundaries",
        "circuits",
        "frontiers",
        "frame",
        "component_basis",
        "machine_nodes",
    }) do
        model[k] = copy(m[k])
    end
    return {
        model = model,
        options = o,
        nodes = nodes,
        bulk = bulk,
        inner = inner,
        outer = outer,
        boundary = boundary,
        moving = moving,
        rotating = rotating,
        marker = marker,
        stats = stats,
        provenance = require("luafemm-cache").signature(m),
    }
end
--- Construct an independent DC position, retaining every bulk triangle.
-- Angle is mechanical, counterclockwise in degrees, relative to the template.
-- Currents replace named prescribed currents; omitted currents retain the bias.
-- @tparam table template Owned capture.
-- @tparam[opt] table o angle (default 0), currents (A), cache (usual options).
-- @treturn table Meshed, unsolved physical-position model in world coordinates.
function S.new(template, o)
    o = o or {}
    for k in pairs(o) do
        check(
            k == "angle" or k == "currents" or k == "cache",
            "unknown position option " .. tostring(k)
        )
    end
    local angle = o.angle or 0
    check(finite(angle), "invalid angle")
    angle = angle % 360
    local m = copy(template.model)
    local spec = copy(template.options)
    m.options.cache = copy(o.cache or { mode = "off" })
    require("luafemm-cache").configure(m)
    local cs, sn = math.cos(math.rad(angle)), math.sin(math.rad(angle))
    local materials = {}
    for i, r in ipairs(m.regions) do
        if template.rotating[i] or i == template.marker then
            for j, contour in ipairs(r.contours) do
                if template.rotating[i] or j == 2 then
                    for k, p in ipairs(contour) do
                        contour[k] = rotate(p, spec, angle)
                    end
                end
            end
            if template.rotating[i] then
                local p = m.materials[r.material]
                if not p.bh and (p.tensor[1] ~= p.tensor[3] or p.tensor[2] ~= 0) then
                    local name = materials[r.material]
                    if not name then
                        name = "sliding-rotor-" .. r.material
                        check(
                            not m.materials[name],
                            "reserved rotated material name already exists"
                        )
                        local q = copy(p)
                        local a, b, c = table.unpack(p.tensor)
                        q.tensor = {
                            cs * cs * a - 2 * cs * sn * b + sn * sn * c,
                            cs * sn * (a - c) + (cs * cs - sn * sn) * b,
                            sn * sn * a + 2 * cs * sn * b + cs * cs * c,
                        }
                        q.library = nil
                        m.materials[name] = q
                        materials[r.material] = name
                    end
                    r.material = name
                end
                r.mx, r.my = cs * r.mx - sn * r.my, sn * r.mx + cs * r.my
                local law = r.magnetization
                local order = law
                        and law.type ~= "uniform"
                        and (law.type == "angular" and law.order or 1)
                    or 0
                r.angle = r.angle + (1 - order) * angle
                if law then
                    local p = rotate({ law.x, law.y }, spec, angle)
                    law.x, law.y = p[1], p[2]
                end
                r.source = nil
            end
        end
    end
    for _, node in pairs(m.machine_nodes or {}) do
        if
            radius(node.east, spec) < spec.inner_radius
            and radius(node.north, spec) < spec.inner_radius
        then
            for _, k in ipairs({ "center", "east", "north" }) do
                node[k] = rotate(node[k], spec, angle)
            end
        end
    end
    check(type(o.currents or {}) == "table", "currents must be a table")
    for name, value in pairs(o.currents or {}) do
        check(m.circuits and m.circuits[name], "unknown circuit " .. tostring(name))
        check(finite(value), "invalid circuit current")
        m.circuits[name].current = value
    end
    require("luafemm-frontiers").compile(m)
    local nodes, tri, boundary = copy(template.nodes), copy(template.bulk), copy(template.boundary)
    local unit = m.options.unit
    local si = { center_x = spec.center_x * unit, center_y = spec.center_y * unit }
    for id in pairs(template.moving) do
        nodes[id] = rotate(nodes[id], si, angle)
    end
    local weights = {}
    for id in ipairs(nodes) do
        weights[id] = template.moving[id] and 1 or 0
    end
    local inner, outer = copy(template.inner), copy(template.outer)
    ringsort(nodes, inner, si)
    ringsort(nodes, outer, si)
    for layer = 1, spec.layers do
        local next_ring = outer
        if layer < spec.layers then
            next_ring = {}
            local r = spec.inner_radius
                + (spec.outer_radius - spec.inner_radius) * layer / spec.layers
            for _, p in ipairs(circle(spec, r)) do
                nodes[#nodes + 1] = { p[1] * unit, p[2] * unit }
                weights[#nodes] = 1 - layer / spec.layers
                next_ring[#next_ring + 1] = #nodes
            end
            ringsort(nodes, next_ring, si)
        end
        bridge(nodes, inner, next_ring, tri)
        inner = next_ring
    end
    m.mesh_stats = copy(template.stats)
    audit(m, nodes, tri, boundary, m.mesh_stats)
    m.mesh_stats.band_triangles = #tri - #template.bulk
    m.mesh_stats.bulk_triangles = #template.bulk
    require("luafemm").prepare_mesh(m, nodes, tri, boundary)
    for i = #template.bulk + 1, #m.triangles do
        local e = m.triangles[i]
        check(
            not e.p.bh
                and e.p.hc == 0
                and e.current == 0
                and e.p.tensor[1] == 1
                and e.p.tensor[2] == 0
                and e.p.tensor[3] == 1,
            "moving-band element is not source-free air"
        )
    end
    m.sliding = {
        method = "moving-band-p1-v1",
        aliases = materials,
        weights = weights,
        angle = angle,
        options = spec,
        provenance = template.provenance,
        nodes = copy(nodes),
        elements = copy(tri),
    }
    -- Store coordinates only, never mesh-owned DOF metadata in the contract.
    for i, p in ipairs(nodes) do
        m.sliding.nodes[i] = { p[1], p[2] }
    end
    m.sliding_contract = require("luafemm-cache").signature(m)
    return m
end
--- Reject edits to an inherited position before solve/cache reuse.
-- @tparam table m Position model; validation never changes its state.
function S.validate(m)
    local s = m.sliding
    check(m.options.frequency == 0 and not m.operating, "sliding positions are static only")
    check(#m.nodes == #s.nodes and #m.triangles == #s.elements, "position mesh changed")
    for i, p in ipairs(m.nodes) do
        check(p[1] == s.nodes[i][1] and p[2] == s.nodes[i][2], "position mesh changed")
    end
    for i, e in ipairs(m.triangles) do
        for j, id in ipairs(e.ids) do
            check(id == s.elements[i][j], "position mesh changed")
        end
    end
    check(
        m.sliding_contract == require("luafemm-cache").signature(m),
        "position declarations changed; create a new position from the template"
    )
end
return S
