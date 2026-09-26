-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- luafemm-mesh implementation.
-- @module luafemm-mesh
-- Original Lua constrained triangulator. No Triangle source is embedded.
-- Incremental Delaunay, segment cavity recovery, constrained Lawson flips.
local G = {}
local abs, sqrt, min, max = math.abs, math.sqrt, math.min, math.max
local predicates = require("luafemm-predicates")
local function check(v, s)
    assert(v, "luafemm mesh: " .. s)
end
local orient = predicates.orient
G.orient = orient
local function neighbor(e, id)
    if e[3] == id then
        return e[4]
    else
        return e[3]
    end
end
local function key(a, b)
    if a > b then
        a, b = b, a
    end
    return a .. ":" .. b
end
local function sorted_keys(t)
    local out = {}
    for k in pairs(t) do
        out[#out + 1] = k
    end
    table.sort(out)
    return out
end
local function distance2(p, a, b)
    local dx, dy = b[1] - a[1], b[2] - a[2]
    local t = max(0, min(1, ((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / (dx * dx + dy * dy)))
    return (p[1] - a[1] - t * dx) ^ 2 + (p[2] - a[2] - t * dy) ^ 2, t
end
local function proper(a, b, c, d)
    local x, y = orient(a, b, c), orient(a, b, d)
    local u, v = orient(c, d, a), orient(c, d, b)
    return ((x < 0 and y > 0) or (x > 0 and y < 0)) and ((u < 0 and v > 0) or (u > 0 and v < 0))
end
--- Validate a simple polygon, rejecting self-contact and backtracking.
-- @tparam table points Array of finite {x,y} model coordinates.
-- @raise Geometry error for an invalid contour.
function G.validate(points)
    check(#points >= 3, "a closed region needs at least three vertices")
    local area = 0
    for i, p in ipairs(points) do
        local q = points[i % #points + 1]
        check(p[1] ~= q[1] or p[2] ~= q[2], "zero-length edge")
        area = area + orient(points[1], p, q)
        local prev = points[(i - 2) % #points + 1]
        check(
            not (
                    orient(prev, p, q) == 0
                    and (predicates.on_segment(q, prev, p) or predicates.on_segment(prev, p, q))
                ),
            "overlapping adjacent edges"
        )
        for j = i + 2, #points do
            if not (i == 1 and j == #points) then
                local a, b = points[j], points[j % #points + 1]
                check(
                    not proper(p, q, a, b)
                        and not predicates.on_segment(p, a, b)
                        and not predicates.on_segment(q, a, b)
                        and not predicates.on_segment(a, p, q)
                        and not predicates.on_segment(b, p, q),
                    "self-intersecting or self-touching region"
                )
            end
        end
    end
    check(area ~= 0, "zero-area region")
end

local incircle = predicates.incircle

--- Construct and audit a constrained Delaunay triangulation.
-- Coordinates are normalised before topological decisions. Intersections
-- remain floating-point constructions; mesh_tolerance controls snapping.
-- @tparam table o Validated model options.
-- @tparam table regions Material polygons with optional local mesh_size.
-- @return Points in model units, counterclockwise index triples, diagnostics.
function G.generate(o, regions)
    local start = os.clock()
    local pred0, pred1 = predicates.stats.orientation, predicates.stats.incircle
    local scale = max(o.xmax - o.xmin, o.ymax - o.ymin)
    local ox, oy = (o.xmin + o.xmax) / 2, (o.ymin + o.ymax) / 2
    local function normalized(p)
        return { (p[1] - ox) / scale, (p[2] - oy) / scale }
    end
    local h = o.h / scale
    local tol = (o.mesh_tolerance and o.mesh_tolerance > 0) and o.mesh_tolerance / scale or 1e-12
    check(tol > 0 and tol < h * 0.01, "mesh tolerance must be less than 1% of mesh size")
    local raw = {}
    local local_regions = {}
    local function loop(points, boundary, size)
        for i, a in ipairs(points) do
            raw[#raw + 1] = {
                normalized(a),
                normalized(points[i % #points + 1]),
                cuts = { 0, 1 },
                boundary = boundary,
                h = size or h,
            }
        end
    end
    loop({ { o.xmin, o.ymin }, { o.xmax, o.ymin }, { o.xmax, o.ymax }, { o.xmin, o.ymax } }, true)
    for _, r in ipairs(regions) do
        local size = r.mesh_size and r.mesh_size > 0 and min(h, r.mesh_size / scale) or h
        loop(r.points, false, size)
        if size < h then
            local poly = {}
            local box = { math.huge, math.huge, -math.huge, -math.huge }
            for _, p in ipairs(r.points) do
                local q = normalized(p)
                poly[#poly + 1] = q
                box[1] = min(box[1], q[1])
                box[2] = min(box[2], q[2])
                box[3] = max(box[3], q[1])
                box[4] = max(box[4], q[2])
            end
            local_regions[#local_regions + 1] = { points = poly, h = size, box = box }
        end
    end
    -- Split proper crossings, T junctions and collinear overlaps before CDT.
    local function cut(s, t)
        if t > 0 and t < 1 then
            s.cuts[#s.cuts + 1] = t
        end
    end
    for i, s in ipairs(raw) do
        for j = i + 1, #raw do
            local t = raw[j]
            local a, b, c, d = s[1], s[2], t[1], t[2]
            if proper(a, b, c, d) then
                local u, v = orient(c, d, a), orient(c, d, b)
                cut(s, u / (u - v))
                u, v = orient(a, b, c), orient(a, b, d)
                cut(t, u / (u - v))
            else
                local ds, u = distance2(c, a, b)
                if ds <= tol * tol then
                    cut(s, u)
                end
                ds, u = distance2(d, a, b)
                if ds <= tol * tol then
                    cut(s, u)
                end
                ds, u = distance2(a, c, d)
                if ds <= tol * tol then
                    cut(t, u)
                end
                ds, u = distance2(b, c, d)
                if ds <= tol * tol then
                    cut(t, u)
                end
            end
        end
    end
    local points, lookup, segments, seen = {}, {}, {}, {}
    local function point(p)
        local ix, iy = math.floor(p[1] / tol), math.floor(p[2] / tol)
        for x = ix - 1, ix + 1 do
            for y = iy - 1, iy + 1 do
                for _, id in ipairs(lookup[x .. "/" .. y] or {}) do
                    local q = points[id]
                    if (p[1] - q[1]) ^ 2 + (p[2] - q[2]) ^ 2 <= tol * tol then
                        return id
                    end
                end
            end
        end
        local k = ix .. "/" .. iy
        check(#points < 40000, "more than 40000 vertices; increase mesh size or curve tolerance")
        local id = #points + 1
        points[id] = p
        lookup[k] = lookup[k] or {}
        table.insert(lookup[k], id)
        return id
    end
    for _, s in ipairs(raw) do
        table.sort(s.cuts)
        local prev
        for _, t in ipairs(s.cuts) do
            if prev and t > prev then
                local a, b = s[1], s[2]
                local length = sqrt((b[1] - a[1]) ^ 2 + (b[2] - a[2]) ^ 2) * (t - prev)
                local n = math.ceil(length / s.h)
                local last = point({ a[1] + prev * (b[1] - a[1]), a[2] + prev * (b[2] - a[2]) })
                for j = 1, n do
                    local u = prev + (t - prev) * j / n
                    local id = point({ a[1] + u * (b[1] - a[1]), a[2] + u * (b[2] - a[2]) })
                    local k = key(last, id)
                    if last ~= id and not seen[k] then
                        segments[#segments + 1] = { last, id, boundary = s.boundary }
                        seen[k] = #segments
                    elseif seen[k] and s.boundary then
                        segments[seen[k]].boundary = true
                    end
                    last = id
                end
            end
            prev = t
        end
    end
    -- Overlapping segments may have different subdivisions. Split them at
    -- every vertex lying on the segment, then deduplicate the atomic edges.
    local atomic = {}
    seen = {}
    for _, s in ipairs(segments) do
        local a, b = points[s[1]], points[s[2]]
        local ids = { { 0, s[1] }, { 1, s[2] } }
        for id, p in ipairs(points) do
            if id ~= s[1] and id ~= s[2] then
                local ds, t = distance2(p, a, b)
                if ds <= tol * tol and t > 0 and t < 1 then
                    ids[#ids + 1] = { t, id }
                end
            end
        end
        table.sort(ids, function(x, y)
            return x[1] < y[1]
        end)
        for i = 2, #ids do
            local a, b = ids[i - 1][2], ids[i][2]
            local k = key(a, b)
            if a ~= b and not seen[k] then
                atomic[#atomic + 1] = { a, b, boundary = s.boundary }
                seen[k] = #atomic
            elseif seen[k] and s.boundary then
                atomic[seen[k]].boundary = true
            end
        end
    end
    segments = atomic
    local lo = normalized({ o.xmin, o.ymin })
    local hi = normalized({ o.xmax, o.ymax })
    local function inside(p, poly)
        local yes = false
        for i, a in ipairs(poly) do
            local b = poly[i % #poly + 1]
            if
                (a[2] > p[2]) ~= (b[2] > p[2])
                and p[1] < (b[1] - a[1]) * (p[2] - a[2]) / (b[2] - a[2]) + a[1]
            then
                yes = not yes
            end
        end
        return yes
    end
    local function target(p)
        local size = h
        for _, r in ipairs(local_regions) do
            if inside(p, r.points) then
                size = min(size, r.h)
            end
        end
        return size
    end
    local function seed(size, box, poly)
        local dy = size * sqrt(3) / 2
        -- A common lattice origin makes overlapping equal-size regions agree.
        local j0 = math.ceil((box[2] - lo[2] - size * 0.5) / dy)
        local j1 = math.floor((box[4] - lo[2] - size * 0.5) / dy)
        for j = j0, j1 do
            local y = lo[2] + size * 0.5 + j * dy
            local x0 = lo[1] + size * (j % 2 == 0 and 1 or 0.5)
            for i = math.ceil((box[1] - x0) / size), math.floor((box[3] - x0) / size) do
                local x = x0 + i * size
                local p = { x, y }
                if (not poly or inside(p, poly)) and target(p) == size then
                    local away = true
                    for _, s in ipairs(raw) do
                        if distance2(p, s[1], s[2]) < (size * 0.18) ^ 2 then
                            away = false
                            break
                        end
                    end
                    if away then
                        point(p)
                    end
                end
            end
        end
    end
    seed(h, { lo[1] + h * 0.2, lo[2] + h * 0.2, hi[1] - h * 0.2, hi[2] - h * 0.2 })
    for _, r in ipairs(local_regions) do
        seed(r.h, r.box, r.points)
    end
    check(#points <= 40000, "more than 40000 vertices; increase mesh size or curve tolerance")
    local count = #points
    points[count + 1] = { -16, -8 }
    points[count + 2] = { 16, -8 }
    points[count + 3] = { 0, 16 }
    local triangles, edges = {}, {}
    local function add(a, b, c)
        local area = orient(points[a], points[b], points[c])
        check(area ~= 0, "degenerate triangle; simplify coincident geometry")
        if area < 0 then
            b, c = c, b
        end
        local id = #triangles + 1
        local t = { a, b, c, alive = true }
        triangles[id] = t
        for i = 1, 3 do
            local u, v = t[i], t[i % 3 + 1]
            local k = key(u, v)
            local e = edges[k]
            if not e then
                e = { u, v }
                edges[k] = e
            end
            check(not e[4], "non-manifold edge")
            if not e[3] then
                e[3] = id
            else
                e[4] = id
            end
        end
        return id
    end
    local function remove(id)
        local t = triangles[id]
        t.alive = false
        for i = 1, 3 do
            local k = key(t[i], t[i % 3 + 1])
            local e = edges[k]
            if e[3] == id then
                e[3] = e[4]
                e[4] = nil
            else
                e[4] = nil
            end
            if not e[3] then
                edges[k] = nil
            end
        end
    end
    local last = add(count + 1, count + 2, count + 3)
    for id = 1, count do
        local p = points[id]
        local at = last
        local located = false
        for walk = 1, #triangles do
            local t = triangles[at]
            local nextid
            for j = 1, 3 do
                local a, b = t[j], t[j % 3 + 1]
                if orient(points[a], points[b], p) < 0 then
                    local e = edges[key(a, b)]
                    nextid = neighbor(e, at)
                    break
                end
            end
            if not nextid then
                located = true
                break
            end
            at = nextid
        end
        check(located, "point location failed")
        local cavity, queue, visited = { [at] = true }, { at }, { [at] = true }
        local q = 1
        while q <= #queue do
            local tid = queue[q]
            q = q + 1
            local t = triangles[tid]
            for j = 1, 3 do
                local e = edges[key(t[j], t[j % 3 + 1])]
                local other = neighbor(e, tid)
                if other and not visited[other] then
                    visited[other] = true
                    local v = triangles[other]
                    local d = incircle(points[v[1]], points[v[2]], points[v[3]], p)
                    if d > 0 then
                        cavity[other] = true
                        queue[#queue + 1] = other
                    end
                end
            end
        end
        local border = {}
        for _, tid in ipairs(queue) do
            local t = triangles[tid]
            for j = 1, 3 do
                local a, b = t[j], t[j % 3 + 1]
                local e = edges[key(a, b)]
                local other = neighbor(e, tid)
                if not other or not cavity[other] then
                    border[#border + 1] = { a, b }
                end
            end
        end
        for _, tid in ipairs(queue) do
            remove(tid)
        end
        for _, e in ipairs(border) do
            last = add(e[1], e[2], id)
        end
    end
    -- Recover each PSLG segment by removing its crossed cavity and
    -- triangulating the two polygons on either side of the segment.
    local function earclip(poly)
        if #poly < 3 then
            return
        end
        local area = 0
        for i, a in ipairs(poly) do
            local b = poly[i % #poly + 1]
            area = area + orient(points[poly[1]], points[a], points[b])
        end
        if area < 0 then
            local rev = {}
            for i = #poly, 1, -1 do
                rev[#rev + 1] = poly[i]
            end
            poly = rev
        end
        while #poly > 3 do
            local found = false
            for i, b in ipairs(poly) do
                local before = (i - 2) % #poly + 1
                local after = i % #poly + 1
                local a, c = poly[before], poly[after]
                if orient(points[a], points[b], points[c]) > 0 then
                    local empty = true
                    for j, p in ipairs(poly) do
                        if j ~= before and j ~= i and j ~= after then
                            if
                                orient(points[a], points[b], points[p]) >= 0
                                and orient(points[b], points[c], points[p]) >= 0
                                and orient(points[c], points[a], points[p]) >= 0
                            then
                                empty = false
                                break
                            end
                        end
                    end
                    if empty then
                        add(a, b, c)
                        table.remove(poly, i)
                        found = true
                        break
                    end
                end
            end
            check(found, "cannot triangulate constraint cavity")
        end
        if #poly == 3 then
            add(poly[1], poly[2], poly[3])
        end
    end
    local fixed = {}
    for _, s in ipairs(segments) do
        local sk = key(s[1], s[2])
        if not edges[sk] then
            local crossed = {}
            for k, e in pairs(edges) do
                if proper(points[s[1]], points[s[2]], points[e[1]], points[e[2]]) then
                    check(not fixed[k], "intersecting constraints were not split")
                    crossed[e[3]] = true
                    if e[4] then
                        crossed[e[4]] = true
                    end
                end
            end
            local adjacency = {}
            local crossed_ids = sorted_keys(crossed)
            for _, tid in ipairs(crossed_ids) do
                local t = triangles[tid]
                for j = 1, 3 do
                    local a, b = t[j], t[j % 3 + 1]
                    local e = edges[key(a, b)]
                    local other = neighbor(e, tid)
                    if not other or not crossed[other] then
                        adjacency[a] = adjacency[a] or {}
                        adjacency[b] = adjacency[b] or {}
                        table.insert(adjacency[a], b)
                        table.insert(adjacency[b], a)
                    end
                end
            end
            check(
                adjacency[s[1]] and adjacency[s[2]],
                "unrecoverable segment; near-coincident vertices"
            )
            local chains = {}
            for _, first in ipairs(adjacency[s[1]]) do
                local chain = { s[1] }
                local prev, at = s[1], first
                for step = 1, count do
                    chain[#chain + 1] = at
                    if at == s[2] then
                        break
                    end
                    local links = adjacency[at]
                    check(links and #links == 2, "invalid cavity boundary")
                    prev, at = at, links[1] == prev and links[2] or links[1]
                end
                check(chain[#chain] == s[2], "open constraint cavity")
                chains[#chains + 1] = chain
            end
            for _, tid in ipairs(crossed_ids) do
                remove(tid)
            end
            for _, chain in ipairs(chains) do
                earclip(chain)
            end
        end
        check(edges[sk], "missing recovered segment")
        fixed[sk] = true
    end
    -- Lawson legalization, keeping all recovered segments fixed.
    local queue = sorted_keys(edges)
    local head, flips = 1, 0
    while head <= #queue do
        local k = queue[head]
        head = head + 1
        local e = edges[k]
        if e and e[4] and not fixed[k] then
            local a, b = e[1], e[2]
            local t, u = triangles[e[3]], triangles[e[4]]
            local c, d
            for j = 1, 3 do
                if t[j] ~= a and t[j] ~= b then
                    c = t[j]
                end
                if u[j] ~= a and u[j] ~= b then
                    d = u[j]
                end
            end
            if proper(points[a], points[b], points[c], points[d]) then
                local val = incircle(points[t[1]], points[t[2]], points[t[3]], points[d])
                if val > 0 then
                    local first, second = e[3], e[4]
                    remove(first)
                    remove(second)
                    add(c, d, a)
                    add(d, c, b)
                    flips = flips + 1
                    check(flips < 20 * count, "edge legalization did not converge")
                    for _, pair in ipairs({ { a, c }, { a, d }, { b, c }, { b, d } }) do
                        queue[#queue + 1] = key(pair[1], pair[2])
                    end
                end
            end
        end
    end
    local output = {}
    local area = 0
    local minangle, maxarea = 180, 0
    for _, t in ipairs(triangles) do
        if t.alive and t[1] <= count and t[2] <= count and t[3] <= count then
            output[#output + 1] = { t[1], t[2], t[3] }
            local a = orient(points[t[1]], points[t[2]], points[t[3]]) / 2
            area = area + a
            maxarea = max(maxarea, a)
            for j = 1, 3 do
                local p, q, r = points[t[j]], points[t[j % 3 + 1]], points[t[(j + 1) % 3 + 1]]
                local dot = (q[1] - p[1]) * (r[1] - p[1]) + (q[2] - p[2]) * (r[2] - p[2])
                minangle = min(minangle, math.deg(math.atan(2 * a, dot)))
            end
        end
    end
    check(
        abs(area - (o.xmax - o.xmin) * (o.ymax - o.ymin) / (scale * scale)) < 1e-8,
        "mesh does not cover the domain"
    )
    for _, s in ipairs(segments) do
        check(edges[key(s[1], s[2])], "lost constrained segment")
    end
    -- Independent audit on the final mesh after removal of the supertriangle.
    local audit, boundary, boundary_edges, used = {}, {}, {}, {}
    for _, s in ipairs(segments) do
        if s.boundary then
            boundary[s[1]] = true
            boundary[s[2]] = true
            boundary_edges[key(s[1], s[2])] = true
        end
    end
    for tid, t in ipairs(output) do
        for j = 1, 3 do
            local a, b = t[j], t[j % 3 + 1]
            local k = key(a, b)
            used[a] = true
            audit[k] = audit[k] or {}
            table.insert(audit[k], tid)
            check(#audit[k] <= 2, "non-manifold final mesh")
        end
    end
    local edge_count, vertex_count = 0, 0
    for _ in pairs(used) do
        vertex_count = vertex_count + 1
    end
    check(vertex_count == count, "isolated or lost vertex")
    for k, incident in pairs(audit) do
        edge_count = edge_count + 1
        local e = edges[k]
        if #incident == 1 then
            check(boundary_edges[k], "unprotected hole in final mesh")
        elseif not fixed[k] then
            local t, u = output[incident[1]], output[incident[2]]
            local c, d
            for j = 1, 3 do
                if t[j] ~= e[1] and t[j] ~= e[2] then
                    c = t[j]
                end
                if u[j] ~= e[1] and u[j] ~= e[2] then
                    d = u[j]
                end
            end
            if proper(points[e[1]], points[e[2]], points[c], points[d]) then
                check(
                    incircle(points[t[1]], points[t[2]], points[t[3]], points[d]) <= 0,
                    "non-Delaunay unconstrained edge"
                )
            end
        end
    end
    check(count - edge_count + #output == 1, "Euler topology check failed")
    for _, s in ipairs(segments) do
        check(audit[key(s[1], s[2])], "constraint not in final mesh")
    end
    points[count + 1] = nil
    points[count + 2] = nil
    points[count + 3] = nil
    for _, p in ipairs(points) do
        p[1] = p[1] * scale + ox
        p[2] = p[2] * scale + oy
    end
    return points,
        output,
        {
            segments = segments,
            seconds = os.clock() - start,
            flips = flips,
            boundary_nodes = boundary,
            min_angle = minangle,
            max_area = maxarea * scale * scale,
            mesh_tolerance = tol * scale,
            exact_orientation = predicates.stats.orientation - pred0,
            exact_incircle = predicates.stats.incircle - pred1,
        }
end
return G
