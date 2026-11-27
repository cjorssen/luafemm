-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Reduced-sector analytical interfaces in separate rotor/stator frames.
-- The bulk mesh stays in each body's reference frame. Relative angle enters
-- only the Fourier trace map; it is not a rotation of the clipped sector mesh.
local Q = {}
local copy = require("luafemm-fem").copy
local function check(v, s)
    assert(v, "luafemm sector gap: " .. s)
end
local function edge(a, b)
    return math.min(a, b) .. ":" .. math.max(a, b)
end
local function theta(p, o, u)
    local t = (math.atan(p[2] / u - o.center_y, p[1] / u - o.center_x) - math.rad(o.start_angle))
        % (2 * math.pi)
    if math.min(t, 2 * math.pi - t) < 1e-9 then
        return 0
    end
    return t % (2 * math.pi)
end
function Q.capture(m, o)
    local M = require("luafemm")
    check(
        not m.nodes and not m.sliding and not m.operating and not m.import_document,
        "use an unmeshed native template"
    )
    check(
        m.options.problem_type == "planar"
            and m.options.field_model == "open"
            and m.options.mesher == "delaunay",
        "use a planar open/delaunay model"
    )
    check(
        m.domain and m.pairs and #m.pairs > 0,
        "declare a sector domain and paired radial frontiers"
    )
    check(not m.network, "network-driven moving sectors are not supported")
    check(not m.options.source and not m.options.boundary, "callbacks cannot be captured")
    local air = m.materials.air
    check(
        air
            and not air.bh
            and air.hc == 0
            and air.conductivity == 0
            and air.tensor[1] == 1
            and air.tensor[2] == 0
            and air.tensor[3] == 1,
        "air must retain vacuum permeability and zero conductivity/coercivity"
    )
    o = copy(o)
    for k in pairs(o) do
        check(
            ({
                inner_radius = true,
                outer_radius = true,
                center_x = true,
                center_y = true,
                segments = true,
                layers = true,
                method = true,
                modes = true,
                sectors = true,
                periodicity = true,
                start_angle = true,
            })[k],
            "unknown option " .. k
        )
    end
    o.center_x, o.center_y, o.start_angle = o.center_x or 0, o.center_y or 0, o.start_angle or 0
    o.segments, o.modes, o.layers = o.segments or 24, o.modes or 0, 1
    o.periodicity = o.periodicity or "periodic"
    for _, k in ipairs({
        "center_x",
        "center_y",
        "start_angle",
        "sectors",
        "segments",
        "modes",
        "inner_radius",
        "outer_radius",
    }) do
        local v = o[k]
        check(type(v) == "number" and v == v and math.abs(v) < math.huge, "invalid " .. k)
    end
    check(
        o.modes % 1 == 0 and o.modes >= 0 and o.modes <= 256,
        "modes must be an integer from 0 to 256"
    )
    check(o.method == "analytical", "reduced sectors require method=analytical")
    check(
        o.sectors % 1 == 0 and o.sectors >= 2 and o.sectors <= 32,
        "sectors must be an integer from 2 to 32"
    )
    check(o.periodicity == "periodic" or o.periodicity == "antiperiodic", "invalid periodicity")
    check(
        o.periodicity == "periodic" or o.sectors % 2 == 0,
        "antiperiodic full machines need an even sector count"
    )
    check(
        o.segments % 1 == 0 and o.segments >= 4 and o.segments <= 128,
        "segments must be an integer from 4 to 128 per sector"
    )
    check(o.inner_radius > 0 and o.outer_radius > o.inner_radius, "require positive nested radii")
    local span = 2 * math.pi / o.sectors
    -- A sector's wedge is convex (N >= 2), so checking every polygon vertex
    -- also checks its straight edges. Extra bulk outside that wedge would form
    -- an unintended magnetic bypass around the analytical interface.
    for _, contour in ipairs(m.domain.contours) do
        for _, p in ipairs(contour) do
            if math.abs(p[1] - o.center_x) + math.abs(p[2] - o.center_y) > 1e-10 then
                check(theta(p, o, 1) <= span + 1e-9, "domain extends outside the declared sector")
            end
        end
    end
    for _, p in ipairs(m.pairs) do
        check(
            p.x == 0
                and p.y == 0
                and math.abs(p.cx - o.center_x) < 1e-10
                and math.abs(p.cy - o.center_y) < 1e-10,
            "pair must rotate about the gap centre without translation"
        )
        check(
            math.abs(math.abs(p.angle) - 360 / o.sectors) < 1e-8,
            "pair rotation must match the sector angle"
        )
        check(
            (p.sign == -1) == (o.periodicity == "antiperiodic"),
            "pair sign must match gap periodicity"
        )
    end
    local paths = { {}, {} }
    for j, r in ipairs({ o.inner_radius, o.outer_radius }) do
        for k = 0, o.segments do
            local t = math.rad(o.start_angle) + span * k / o.segments
            paths[j][#paths[j] + 1] = { o.center_x + r * math.cos(t), o.center_y + r * math.sin(t) }
        end
    end
    -- Require an air collar: no physical region can straddle either trace.
    for _, r in ipairs(m.regions) do
        local lo, hi = math.huge, 0
        for _, c in ipairs(r.contours) do
            for j, p in ipairs(c) do
                local d = math.sqrt((p[1] - o.center_x) ^ 2 + (p[2] - o.center_y) ^ 2)
                local distance = require("luafemm-topology").distance(
                    { o.center_x, o.center_y },
                    p,
                    c[j % #c + 1]
                )
                lo = math.min(lo, distance)
                hi = math.max(hi, d)
            end
        end
        check(
            hi < o.inner_radius * math.cos(span / o.segments / 2) or lo > o.outer_radius,
            "keep physical regions away from the source-free gap"
        )
    end
    local contour = {}
    for _, p in ipairs(paths[2]) do
        contour[#contour + 1] = p
    end
    for j = #paths[1], 1, -1 do
        contour[#contour + 1] = paths[1][j]
    end
    M.region(m, "air", contour)
    local marker = #m.regions
    M.mesh(m)
    local bulk, used, edges = {}, {}, {}
    for _, e in ipairs(m.triangles) do
        if m.topology.faces[e.face].region ~= marker then
            bulk[#bulk + 1] = copy(e.ids)
            for j, a in ipairs(e.ids) do
                used[a] = true
                local b = e.ids[j % 3 + 1]
                local key = edge(a, b)
                local v = edges[key] or { a, b, count = 0 }
                edges[key] = v
                v.count = v.count + 1
            end
        end
    end
    local sets = { {}, {} }
    local distance = require("luafemm-topology").distance
    for _, e in pairs(edges) do
        if e.count == 1 then
            for _, id in ipairs({ e[1], e[2] }) do
                local p = { m.nodes[id][1] / m.options.unit, m.nodes[id][2] / m.options.unit }
                for j, path in ipairs(paths) do
                    for k = 1, #path - 1 do
                        if distance(p, path[k], path[k + 1]) < 1e-7 then
                            sets[j][id] = true
                        end
                    end
                end
            end
        end
    end
    local nodes, map, boundary, rings = {}, {}, {}, { {}, {} }
    for id in ipairs(m.nodes) do
        if used[id] then
            map[id] = #nodes + 1
            nodes[#nodes + 1] = copy(m.nodes[id])
            if m.exterior_nodes[id] then
                boundary[map[id]] = true
            end
            for j = 1, 2 do
                if sets[j][id] then
                    rings[j][#rings[j] + 1] = map[id]
                end
            end
        end
    end
    for j, ring in ipairs(rings) do
        table.sort(ring, function(a, b)
            return theta(nodes[a], o, m.options.unit) < theta(nodes[b], o, m.options.unit)
        end)
        check(#ring >= o.segments + 1, "incomplete sector trace")
        check(
            math.abs(theta(nodes[ring[#ring]], o, m.options.unit) - span) < 1e-7,
            "trace does not end at the sector boundary"
        )
        local radius = (j == 1 and o.inner_radius or o.outer_radius) * m.options.unit
        for _, id in ipairs(ring) do
            local p = nodes[id]
            local x, y = p[1] - o.center_x * m.options.unit, p[2] - o.center_y * m.options.unit
            local length = math.sqrt(x * x + y * y)
            p[1], p[2] =
                o.center_x * m.options.unit + radius * x / length,
                o.center_y * m.options.unit + radius * y / length
        end
    end
    for _, e in ipairs(bulk) do
        for j, id in ipairs(e) do
            e[j] = map[id]
        end
    end
    local tri = copy(bulk)
    local i, j = 1, 1
    local pred = require("luafemm-predicates")
    while i < #rings[1] or j < #rings[2] do
        local a, b = rings[1][i], rings[2][j]
        local c, d = rings[1][i + 1], rings[2][j + 1]
        local left = c and pred.orient(nodes[a], nodes[b], nodes[c]) > 0
        local right = d and pred.orient(nodes[a], nodes[b], nodes[d]) > 0
        check(left or right, "invalid sector strip")
        if
            left
            and (
                not right
                or theta(nodes[c], o, m.options.unit) <= theta(nodes[d], o, m.options.unit)
            )
        then
            tri[#tri + 1] = { a, b, c }
            i = i + 1
        else
            tri[#tri + 1] = { a, b, d }
            j = j + 1
        end
    end
    local stats = copy(m.mesh_stats)
    stats.min_angle, stats.max_area = 180, 0
    for _, e in ipairs(tri) do
        local twice = pred.orient(nodes[e[1]], nodes[e[2]], nodes[e[3]])
        check(twice > 0, "projected sector trace inverts a bulk element; refine the air collar")
        stats.max_area = math.max(stats.max_area, twice / (2 * m.options.unit ^ 2))
        for j = 1, 3 do
            local a, b, c = nodes[e[j]], nodes[e[j % 3 + 1]], nodes[e[(j + 1) % 3 + 1]]
            local dot = (b[1] - a[1]) * (c[1] - a[1]) + (b[2] - a[2]) * (c[2] - a[2])
            stats.min_angle = math.min(stats.min_angle, math.deg(math.atan(twice, dot)))
        end
    end
    check(
        m.options.min_angle == 0 or stats.min_angle + 1e-8 >= m.options.min_angle,
        "sector strip misses minimum angle"
    )
    check(
        m.options.max_area == 0 or stats.max_area <= m.options.max_area * (1 + 1e-8),
        "sector strip misses maximum area"
    )
    stats.segments = {}
    stats.boundary_nodes = copy(boundary)
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
        "domain",
        "pairs",
        "frame",
        "component_basis",
    }) do
        model[k] = copy(m[k])
    end
    return {
        sector = true,
        model = model,
        options = o,
        nodes = nodes,
        elements = tri,
        boundary = boundary,
        inner = rings[1],
        outer = rings[2],
        bulk_count = #bulk,
        stats = stats,
        provenance = require("luafemm-cache").signature(m),
    }
end
function Q.new(t, o)
    o = o or {}
    for k in pairs(o) do
        check(
            k == "angle" or k == "cache" or k == "currents" or k == "phases" or k == "numerical",
            "unknown position option " .. k
        )
    end
    local m = copy(t.model)
    local angle = o.angle or 0
    check(
        type(angle) == "number" and angle == angle and math.abs(angle) < math.huge,
        "invalid angle"
    )
    m.options.cache = copy(o.cache or { mode = "off" })
    require("luafemm-cache").configure(m)
    require("luafemm-linear").override(m, o.numerical)
    for name, v in pairs(o.currents or {}) do
        check(
            m.circuits and m.circuits[name] and m.circuits[name].drive == "current",
            "unknown prescribed current"
        )
        check(type(v) == "number" and v == v and math.abs(v) < math.huge, "invalid current")
        m.circuits[name].current = v
    end
    for name, v in pairs(o.phases or {}) do
        check(m.circuits and m.circuits[name], "unknown circuit phase")
        m.circuits[name].phase = v
    end
    require("luafemm-frontiers").compile(m)
    m.mesh_stats = copy(t.stats)
    require("luafemm").prepare_mesh(m, copy(t.nodes), copy(t.elements), copy(t.boundary))
    m.sliding = {
        sector = true,
        method = "sector-fourier-dtn-v1",
        options = copy(t.options),
        inner = copy(t.inner),
        outer = copy(t.outer),
        angle = angle,
        bulk_count = t.bulk_count,
        nodes = copy(t.nodes),
        elements = copy(t.elements),
        aliases = {},
        provenance = t.provenance,
    }
    require("luafemm-airgap").prepare(m)
    m.sliding_contract = require("luafemm-cache").signature(m)
    return m
end
return Q
