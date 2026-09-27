-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Confined planar magnetostatics with floating, equipotential boundary loops.
-- Each hole contributes one unknown Az and one prescribed enclosed current.
-- Condensing all vertices of a loop into one degree of freedom enforces zero
-- normal flux exactly. Interior P1 continuity conserves flux across elements.
-- @module luafemm-ideal
local I = {}
local topology = require("luafemm-topology")
local geometry = require("luafemm-geometry")
local function check(v, s)
    assert(v, "luafemm ideal: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end

--- Add an enclosed-current seed before meshing, in model coordinates and amperes.
-- The seed must lie strictly in a bounded excluded window, not in iron or air gaps.
-- Positive current points along +z. Several seeds in a window add algebraically.
function I.excitation(m, x, y, turns)
    check(m.options.field_model == "ideal", "excitation seeds require field model=ideal")
    check(not m.nodes and not m.topology, "declare excitations before meshing")
    check(finite(x) and finite(y) and finite(turns), "invalid excitation seed")
    m.ideal_sources = m.ideal_sources or {}
    m.ideal_sources[#m.ideal_sources + 1] = { x, y, turns }
end

--- Compile the paint-order ideal domain before cache lookup or mesh generation.
-- Current-bearing excluded faces remain available to determine enclosed currents.
function I.compile(m)
    if m.options.field_model ~= "ideal" or m.ideal_compiled then
        return
    end
    check(not m.topology, "FEM import cannot be combined with field model=ideal")
    check(m.options.mesher == "delaunay", "the ideal domain requires mesher=delaunay")
    check(
        not m.options.source and not m.options.boundary,
        "callbacks are unavailable in ideal mode"
    )
    for _, b in pairs(m.boundaries) do
        check(
            b.type == "dirichlet" and b.value == 0 and b.dx == 0 and b.dy == 0,
            "rectangular boundary overrides are unavailable in ideal mode"
        )
    end
    local t = topology.native(m)
    local active = 0
    for _, f in ipairs(t.faces) do
        local included = false
        for _, r in ipairs(m.regions) do
            if geometry.contains(r, f.seed[1], f.seed[2]) then
                included = r.ideal_domain == true
            end
        end
        f.hole = not included
        if included then
            active = active + 1
            check(f.current == 0, "windings must lie outside the confined domain")
        end
    end
    check(active > 0, "declare at least one region with ideal domain=true")
    m.topology, m.ideal_compiled = t, true
end

--- Recover oriented boundary cycles from counterclockwise triangles.
-- Positive area is an exterior boundary; negative area encloses a window.
-- Reject touching components and open/nonmanifold boundaries explicitly.
function I.loops(nodes, elements)
    local edges = {}
    for _, e in ipairs(elements) do
        for j = 1, 3 do
            local a, b = e[j], e[j % 3 + 1]
            local key = math.min(a, b) .. ":" .. math.max(a, b)
            local edge = edges[key] or { a, b, count = 0 }
            edge.count = edge.count + 1
            check(edge.count <= 2, "nonmanifold mesh edge")
            edges[key] = edge
        end
    end
    local successor, incoming = {}, {}
    for _, e in pairs(edges) do
        if e.count == 1 then
            check(not successor[e[1]] and not incoming[e[2]], "touching or nonmanifold domain")
            successor[e[1]], incoming[e[2]] = e[2], e[1]
        end
    end
    local loops, visited = {}, {}
    for start = 1, #nodes do
        if successor[start] and not visited[start] then
            local ids, points, id = {}, {}, start
            repeat
                check(id and not visited[id], "open boundary")
                visited[id] = true
                ids[#ids + 1], points[#points + 1] = id, nodes[id]
                id = successor[id]
            until id == start
            local area = topology.area(points)
            check(area ~= 0, "zero-area boundary")
            loops[#loops + 1] = { ids = ids, contours = { points }, area = area }
        end
    end
    return loops
end

--- Return fixed outer values, boundary representatives and enclosed-current loads.
-- Nodes use SI coordinates. Loads use amperes, before multiplication by mu0.
function I.prepare(m, nodes, elements)
    local loops = I.loops(nodes, elements)
    local fixed, representative, loads = {}, {}, {}
    local unit = m.options.unit
    for _, loop in ipairs(loops) do
        local first = loop.ids[1]
        if loop.area > 0 then
            for _, id in ipairs(loop.ids) do
                fixed[id] = 0
            end
        else
            loads[first] = 0
            for _, id in ipairs(loop.ids) do
                representative[id] = first
            end
            for _, face in ipairs(m.topology.faces) do
                if
                    face.hole
                    and face.current ~= 0
                    and geometry.contains(loop, face.seed[1] * unit, face.seed[2] * unit)
                then
                    loads[first] = loads[first] + face.current * face.area * unit ^ 2
                end
            end
        end
    end
    for _, source in ipairs(m.ideal_sources or {}) do
        local face = topology.locate(m.topology, source[1], source[2], true)
        check(face and face.hole, "excitation seed must be strictly inside an excluded window")
        local found = false
        for _, loop in ipairs(loops) do
            if loop.area < 0 and geometry.contains(loop, source[1] * unit, source[2] * unit) then
                loads[loop.ids[1]] = loads[loop.ids[1]] + source[3]
                found = true
            end
        end
        check(found, "excitation seed lies in exterior air, not a bounded window")
    end
    return fixed, representative, loads, loops
end

--- Check cached boundary constants independently of the nonlinear residual.
function I.validate_potential(m, nodes, elements, A)
    local fixed, representatives = I.prepare(m, nodes, elements)
    for id, value in ipairs(A) do
        check(fixed[id] == nil or value == fixed[id], "invalid exterior potential")
        check(
            not representatives[id] or value == A[representatives[id]],
            "cached potential violates flux confinement"
        )
    end
end

--- Oriented flux through a path, positive along its right normal, in webers.
-- Constant depth is in model units. Endpoint Az differences integrate the P1
-- normal field exactly; intermediate vertices do not change the result.
-- @tparam table m Source model, solved on demand.
-- @tparam table points Section polyline in model units, copied before projection.
-- @tparam[opt=0] number snap_tolerance Maximum endpoint-to-wall distance in model units.
-- @treturn number Signed right-normal flux in Wb.
function I.flux(m, points, snap_tolerance)
    local f = require("luafemm")
    check(m.options.depth and m.options.depth > 0, "set depth before requesting flux in Wb")
    check(#points >= 2, "a section needs two endpoints")
    snap_tolerance = snap_tolerance or 0
    check(finite(snap_tolerance) and snap_tolerance >= 0, "invalid endpoint snap tolerance")
    local owned = {}
    for _, p in ipairs(points) do
        check(type(p) == "table" and finite(p[1]) and finite(p[2]), "invalid section vertex")
        owned[#owned + 1] = { p[1], p[2] }
    end
    points = owned
    f.solve(m)
    if snap_tolerance > 0 then
        local loops = m.ideal_loops
        if not loops then
            local elements = {}
            for _, e in ipairs(m.triangles) do
                elements[#elements + 1] = e.ids
            end
            loops = I.loops(m.nodes, elements)
        end
        for _, index in ipairs({ 1, #points }) do
            local p, best, nearest = points[index], snap_tolerance
            local si = { p[1] * m.options.unit, p[2] * m.options.unit }
            for _, loop in ipairs(loops) do
                for j, id in ipairs(loop.ids) do
                    local a, b = m.nodes[id], m.nodes[loop.ids[j % #loop.ids + 1]]
                    local distance, t = topology.distance(si, a, b)
                    distance = distance / m.options.unit
                    if distance < best then
                        t = math.max(0, math.min(1, t))
                        best = distance
                        nearest = {
                            (a[1] + t * (b[1] - a[1])) / m.options.unit,
                            (a[2] + t * (b[2] - a[2])) / m.options.unit,
                        }
                    end
                end
            end
            if nearest then
                points[index] = nearest
            end
        end
    end
    for j = 2, #points do
        require("luafemm-integrals").segment(m, points[j - 1], points[j])
    end
    local first, last = points[1], points[#points]
    local _, _, a = f.sample(m, first[1], first[2])
    local _, _, b = f.sample(m, last[1], last[2])
    return (b - a) * m.options.depth * m.options.unit
end

--- Select the geometric or isolated half-flux cycle of a named component.
-- Geometric contours are oriented with their net B circulation. A flux cycle
-- must enclose exactly one window and may not accidentally enclose another.
function I.mean(m, component, cycle, kind)
    local f = require("luafemm")
    local cycles = m.component_cycles and m.component_cycles[component]
    check(cycles and cycles[cycle], "unknown component cycle " .. component .. "/" .. cycle)
    local points = {}
    for _, p in ipairs(cycles[cycle]) do
        points[#points + 1] = { p[1], p[2] }
    end
    f.solve(m)
    if kind == "flux" then
        check(m.options.field_model == "ideal", "half-flux contours require ideal mode")
        local selected = {}
        local region = { contours = { points } }
        for _, loop in ipairs(m.ideal_loops) do
            local p = m.nodes[loop.ids[1]]
            if
                loop.area < 0
                and geometry.contains(region, p[1] / m.options.unit, p[2] / m.options.unit)
            then
                selected[#selected + 1] = loop
            end
        end
        check(#selected == 1, "a half-flux cycle must enclose exactly one window")
        local loop = selected[1]
        local value = m.A[loop.ids[1]] / 2
        check(math.abs(value) > 1e-20, "no half-flux contour in a zero-flux window")
        local candidates = {}
        for _, path in ipairs(f.contours(m, 1, { value })) do
            local poly = {}
            for _, p in ipairs(path) do
                poly[#poly + 1] = { p.x, p.y }
            end
            local p, q = poly[1], poly[#poly]
            if (p[1] - q[1]) ^ 2 + (p[2] - q[2]) ^ 2 < 1e-16 then
                local contains, others = false, 0
                for _, other in ipairs(m.ideal_loops) do
                    local v = m.nodes[other.ids[1]]
                    if
                        other.area < 0
                        and geometry.contains(
                            { contours = { poly } },
                            v[1] / m.options.unit,
                            v[2] / m.options.unit
                        )
                    then
                        if other == loop then
                            contains = true
                        else
                            others = others + 1
                        end
                    end
                end
                if contains and others == 0 then
                    candidates[#candidates + 1] = poly
                end
            end
        end
        check(#candidates == 1, "no isolated half-flux line for this cycle; use kind=geometric")
        return candidates[1]
    end
    check(kind == "geometric", "unknown mean contour kind")
    local sign = 0
    for i = 2, #points do
        local a, b = points[i - 1], points[i]
        local bx, by = f.sample(m, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2)
        sign = sign + bx * (b[1] - a[1]) + by * (b[2] - a[2])
    end
    if sign < 0 then
        local reverse = {}
        for i = #points, 1, -1 do
            reverse[#reverse + 1] = points[i]
        end
        points = reverse
    end
    return points
end
return I
