-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Planar arrangements shared by TikZ geometry and magnetic file interchange.
-- Original implementation: split constraints, walk oriented half-edges, nest
-- bounded cycles, then assign labels. All coordinates are in model units.
-- @module luafemm-topology
local T = {}
local geom = require("luafemm-geometry")
local pred = require("luafemm-predicates")
local function check(v, s)
    assert(v, "luafemm topology: " .. s)
end
local function area(p)
    local sum = 0
    for i, a in ipairs(p) do
        local b = p[i % #p + 1]
        sum = sum + pred.orient(p[1], a, b) / 2
    end
    return sum
end
T.area = area
local function distance(p, a, b)
    local x, y = b[1] - a[1], b[2] - a[2]
    local t = ((p[1] - a[1]) * x + (p[2] - a[2]) * y) / (x * x + y * y)
    local u = math.max(0, math.min(1, t))
    return math.sqrt((p[1] - a[1] - u * x) ^ 2 + (p[2] - a[2] - u * y) ^ 2), t
end
T.distance = distance

--- Compile a list of {point,point,condition=...} into bounded polygonal faces.
-- Tags are copied to every split segment; conflicting tags are rejected.
-- @tparam table raw Segments with model-coordinate endpoints and optional tags.
-- @tparam[opt] number tolerance Positive vertex merging distance in model units.
-- @treturn table Owned points, split segments and bounded faces.
function T.arrange(raw, tolerance)
    local eps = tolerance or 1e-9
    local points, segments, seen = {}, {}, {}
    local function vertex(p)
        for i, q in ipairs(points) do
            if (p[1] - q[1]) ^ 2 + (p[2] - q[2]) ^ 2 <= eps ^ 2 then
                return i
            end
        end
        points[#points + 1] = { p[1], p[2] }
        return #points
    end
    local work = {}
    for _, s in ipairs(raw) do
        check(distance(s[1], s[1], s[2]) == 0, "zero-length or invalid segment")
        local a, b = vertex(s[1]), vertex(s[2])
        check(a ~= b, "segment collapsed at geometry tolerance")
        work[#work + 1] = { a, b, condition = s.condition, group = s.group, hidden = s.hidden }
    end
    for i, s in ipairs(work) do
        local a, b = points[s[1]], points[s[2]]
        for j = i + 1, #work do
            local t = work[j]
            local c, d = points[t[1]], points[t[2]]
            local u, v = pred.orient(a, b, c), pred.orient(a, b, d)
            local w, z = pred.orient(c, d, a), pred.orient(c, d, b)
            if u * v < 0 and w * z < 0 then
                local f = w / (w - z)
                vertex({ a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2]) })
            end
        end
    end
    local function same(a, b)
        if not a or not b then
            return true
        end
        for k, v in pairs(a) do
            if b[k] ~= v then
                return false
            end
        end
        for k, v in pairs(b) do
            if a[k] ~= v then
                return false
            end
        end
        return true
    end
    for _, s in ipairs(work) do
        local a, b = points[s[1]], points[s[2]]
        local cuts = {}
        for i, p in ipairs(points) do
            local d, t = distance(p, a, b)
            if d <= eps and t >= -eps and t <= 1 + eps then
                cuts[#cuts + 1] = { t, i }
            end
        end
        table.sort(cuts, function(x, y)
            return x[1] < y[1]
        end)
        for j = 2, #cuts do
            local u, v = cuts[j - 1][2], cuts[j][2]
            local k = math.min(u, v) .. ":" .. math.max(u, v)
            if not seen[k] then
                local e = {
                    u,
                    v,
                    condition = s.condition,
                    group = s.group or 0,
                    hidden = s.hidden or 0,
                }
                segments[#segments + 1] = e
                seen[k] = e
            else
                check(
                    same(seen[k].condition, s.condition),
                    "conflicting conditions on coincident edges"
                )
                seen[k].condition = seen[k].condition or s.condition
            end
        end
    end
    local adj = {}
    for i in ipairs(points) do
        adj[i] = {}
    end
    for _, s in ipairs(segments) do
        adj[s[1]][#adj[s[1]] + 1] = s[2]
        adj[s[2]][#adj[s[2]] + 1] = s[1]
    end
    for i, neighbors in ipairs(adj) do
        check(#neighbors >= 2, "open or isolated constraint at vertex " .. i)
        table.sort(neighbors, function(a, b)
            return math.atan(points[a][2] - points[i][2], points[a][1] - points[i][1])
                < math.atan(points[b][2] - points[i][2], points[b][1] - points[i][1])
        end)
    end
    local visited, faces = {}, {}
    for _, s in ipairs(segments) do
        for direction = 1, 2 do
            local u, v = s[direction], s[3 - direction]
            local first = u .. ":" .. v
            if not visited[first] then
                local polygon, ids = {}, {}
                repeat
                    local k = u .. ":" .. v
                    check(not visited[k], "invalid face walk")
                    visited[k] = true
                    polygon[#polygon + 1] = points[u]
                    ids[#ids + 1] = u
                    local n = adj[v]
                    local pos
                    for j, id in ipairs(n) do
                        if id == u then
                            pos = j
                            break
                        end
                    end
                    u, v = v, n[(pos - 2) % #n + 1]
                until u .. ":" .. v == first
                local a = area(polygon)
                if a > eps * eps then
                    faces[#faces + 1] = {
                        points = polygon,
                        contours = { polygon },
                        fill_rule = "even odd",
                        area = a,
                        ids = ids,
                    }
                end
            end
        end
    end
    -- A dangling/bridging constraint is not a region contour. Reject it
    -- explicitly rather than dropping it when adapting faces to the mesher.
    local face_edges = {}
    for _, f in ipairs(faces) do
        local used = {}
        for i, a in ipairs(f.ids) do
            check(not used[a], "bridged or self-touching face boundary is unsupported")
            used[a] = true
            local b = f.ids[i % #f.ids + 1]
            face_edges[math.min(a, b) .. ":" .. math.max(a, b)] = true
        end
    end
    for _, s in ipairs(segments) do
        check(
            face_edges[math.min(s[1], s[2]) .. ":" .. math.max(s[1], s[2])],
            "open or bridging constraint does not bound a valid face"
        )
    end
    -- An interior point close to an edge avoids centroid failures on concavities.
    local function seed(face)
        for i, a in ipairs(face.points) do
            local b = face.points[i % #face.points + 1]
            local x, y = (a[1] + b[1]) / 2, (a[2] + b[2]) / 2
            local dx, dy = b[1] - a[1], b[2] - a[2]
            local length = math.sqrt(dx * dx + dy * dy)
            local clearance = length / 4
            for _, e in ipairs(segments) do
                local d = distance({ x, y }, points[e[1]], points[e[2]])
                if d > eps then
                    clearance = math.min(clearance, d / 4)
                end
            end
            local p = { x - dy / length * clearance, y + dx / length * clearance }
            if geom.contains(face, p[1], p[2]) then
                return p
            end
        end
        error("luafemm topology: cannot locate an interior label")
    end
    for _, f in ipairs(faces) do
        f.seed = seed(f)
    end
    for i, f in ipairs(faces) do
        local parent
        for j, g in ipairs(faces) do
            if
                i ~= j
                and g.area > f.area
                and geom.contains(g, f.seed[1], f.seed[2])
                and (not parent or g.area < faces[parent].area)
            then
                parent = j
            end
        end
        f.parent = parent
    end
    for _, f in ipairs(faces) do
        if f.parent then
            local g = faces[f.parent]
            g.contours[#g.contours + 1] = f.points
        end
    end
    for i, f in ipairs(faces) do
        f.id = i
        f.seed = seed(f)
        local a = f.area
        for _, g in ipairs(faces) do
            if g.parent == i then
                a = a - g.area
            end
        end
        f.area = a
        check(a > eps * eps, "collapsed face")
    end
    return { points = points, segments = segments, faces = faces, tolerance = eps }
end

--- Find one bounded face, with explicit rejection for ambiguous boundary seeds.
-- @tparam table t Compiled arrangement.
-- @tparam number x Model-coordinate abscissa.
-- @tparam number y Model-coordinate ordinate.
-- @tparam[opt] boolean strict Reject seeds on an interface.
-- @treturn table Face, or nil outside all bounded faces.
function T.locate(t, x, y, strict)
    if strict then
        for _, s in ipairs(t.segments) do
            check(
                distance({ x, y }, t.points[s[1]], t.points[s[2]]) > t.tolerance,
                "label lies on a boundary"
            )
        end
    end
    for _, f in ipairs(t.faces) do
        if geom.contains(f, x, y) then
            return f
        end
    end
end

--- Compile native paint-order contours without meshing or solving.
-- @tparam table m Native model with ordered material contours.
-- @treturn table Computational arrangement; the model is not mutated.
function T.native(m)
    local o = m.options
    local raw = {}
    local function loop(p, conditions)
        for i, a in ipairs(p) do
            raw[#raw + 1] = { a, p[i % #p + 1], condition = conditions and conditions[i] }
        end
    end
    loop(
        { { o.xmin, o.ymin }, { o.xmax, o.ymin }, { o.xmax, o.ymax }, { o.xmin, o.ymax } },
        { m.boundaries.bottom, m.boundaries.right, m.boundaries.top, m.boundaries.left }
    )
    for _, r in ipairs(m.regions) do
        for _, p in ipairs(r.contours) do
            loop(p)
        end
    end
    local t = T.arrange(
        raw,
        o.mesh_tolerance and o.mesh_tolerance > 0 and o.mesh_tolerance
            or math.max(o.xmax - o.xmin, o.ymax - o.ymin) * 1e-12
    )
    for _, f in ipairs(t.faces) do
        f.material, f.current, f.angle, f.mx, f.my, f.mesh_size = "air", 0, 0, 1, 0, 0
        for _, r in ipairs(m.regions) do
            if geom.contains(r, f.seed[1], f.seed[2]) then
                for _, k in ipairs({ "material", "current", "angle", "mx", "my", "mesh_size" }) do
                    f[k] = r[k]
                end
            end
        end
    end
    return t
end

--- Remove excluded/outside triangles from a constrained rectangular staging mesh.
-- The generator's own audit precedes carving; exposed edges are checked again.
-- @tparam table t Tagged arrangement with explicit hole flags.
-- @tparam table points Staging vertices in model units.
-- @tparam table elements Staging triangle triples; consumed by this operation.
-- @tparam table stats Staging mesh diagnostics, updated in place.
-- @return Compacted vertices, retained elements and updated diagnostics.
function T.carve(t, points, elements, stats)
    local used, kept = {}, {}
    for _, e in ipairs(elements) do
        local a, b, c = points[e[1]], points[e[2]], points[e[3]]
        local f = T.locate(t, (a[1] + b[1] + c[1]) / 3, (a[2] + b[2] + c[2]) / 3)
        if f and not f.hole then
            kept[#kept + 1] = e
            for _, id in ipairs(e) do
                used[id] = true
            end
        end
    end
    local out, map = {}, {}
    for id, p in ipairs(points) do
        if used[id] then
            map[id] = #out + 1
            out[#out + 1] = p
        end
    end
    local edges = {}
    for _, e in ipairs(kept) do
        for j = 1, 3 do
            e[j] = map[e[j]]
        end
        for j = 1, 3 do
            local a, b = e[j], e[j % 3 + 1]
            local key = math.min(a, b) .. ":" .. math.max(a, b)
            local edge = edges[key] or { a, b, count = 0 }
            edge.count = edge.count + 1
            edges[key] = edge
        end
    end
    local constraints, boundary = {}, {}
    for _, s in ipairs(stats.segments) do
        local a, b = map[s[1]], map[s[2]]
        if a and b and edges[math.min(a, b) .. ":" .. math.max(a, b)] then
            constraints[#constraints + 1] = { a, b }
        end
    end
    for _, e in pairs(edges) do
        if e.count == 1 then
            local p, q = out[e[1]], out[e[2]]
            local found = false
            for _, s in ipairs(t.segments) do
                if
                    distance(p, t.points[s[1]], t.points[s[2]]) < t.tolerance * 20
                    and distance(q, t.points[s[1]], t.points[s[2]]) < t.tolerance * 20
                then
                    found = true
                    break
                end
            end
            check(found, "unprotected carved boundary")
            boundary[e[1]], boundary[e[2]] = true, true
        end
    end
    check(#kept > 0, "empty computational domain")
    local actual, edge_count, minangle, maxarea = 0, 0, 180, 0
    for _ in pairs(edges) do
        edge_count = edge_count + 1
    end
    for _, e in ipairs(kept) do
        local a, b, c = out[e[1]], out[e[2]], out[e[3]]
        local twice = pred.orient(a, b, c)
        actual = actual + twice / 2
        maxarea = math.max(maxarea, twice / 2)
        for j = 1, 3 do
            local p, q, r = out[e[j]], out[e[j % 3 + 1]], out[e[(j + 1) % 3 + 1]]
            local dot = (q[1] - p[1]) * (r[1] - p[1]) + (q[2] - p[2]) * (r[2] - p[2])
            minangle = math.min(minangle, math.deg(math.atan(twice, dot)))
        end
    end
    local expected = 0
    for _, f in ipairs(t.faces) do
        if not f.hole then
            expected = expected + f.area
        end
    end
    check(math.abs(actual / expected - 1) < 1e-8, "carved domain coverage mismatch")
    check(#out - edge_count + #kept == T.characteristic(t), "carved domain Euler check failed")
    stats.min_angle, stats.max_area = minangle, maxarea
    stats.segments, stats.boundary_nodes = constraints, boundary
    return out, kept, stats
end
--- Euler characteristic of the active polygonal cell complex.
-- @tparam table t Compiled arrangement with active/hole faces.
-- @treturn number Euler characteristic C-H of the active domain.
function T.characteristic(t)
    local vertices, edges, chi = {}, {}, 0
    local function key(p)
        return string.format("%.17g/%.17g", p[1], p[2])
    end
    local function add(poly, v, e)
        for i, a in ipairs(poly) do
            local ka, kb = key(a), key(poly[i % #poly + 1])
            v[ka] = true
            e[ka < kb and ka .. ":" .. kb or kb .. ":" .. ka] = true
        end
    end
    local function difference(v, e)
        local n = 0
        for _ in pairs(v) do
            n = n + 1
        end
        for _ in pairs(e) do
            n = n - 1
        end
        return n
    end
    for _, f in ipairs(t.faces) do
        if not f.hole then
            -- Child cycles may share edges. Their union, not their raw count,
            -- determines the number of holes in this face's boundary complex.
            local hv, he = {}, {}
            for i, p in ipairs(f.contours) do
                add(p, vertices, edges)
                if i > 1 then
                    add(p, hv, he)
                end
            end
            chi = chi + 1 - (difference(hv, he) + #f.contours - 1)
        end
    end
    return difference(vertices, edges) + chi
end

return T
