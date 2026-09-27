-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Bounded quality refinement of an existing constrained Delaunay mesh.
-- Original Lua implementation of segment-first Delaunay refinement, inspired
-- by the algorithms studied in Triangle. No Triangle source is embedded.
-- Acute input corners are not specially protected: unattainable requests fail
-- at the insertion budget or geometric tolerance, never silently relax quality.
-- @module luafemm-refine
local R = {}
local predicates = require("luafemm-predicates")
local orient = predicates.orient

local function check(value, message)
    assert(value, "luafemm mesh quality: " .. message)
end
local function key(a, b)
    return math.min(a, b) .. ":" .. math.max(a, b)
end
local function opposite(t, a, b)
    for j = 1, 3 do
        if t[j] ~= a and t[j] ~= b then
            return t[j]
        end
    end
end

--- Refine a mesh in place; called only when a quality target is nonzero.
-- @tparam table mesh Normalised points, alive CCW triangles, edge incidences,
--   fixed constraint keys and segments, plus point/add/remove/legalize callbacks.
-- @tparam table options min_angle (degrees), max_area (normalised area),
--   max_steiner (additional vertices) and tolerance (normalised distance).
-- @treturn table Counts of additional vertices, segment splits and circumcentres.
function R.run(mesh, options)
    local points, triangles, edges = mesh.points, mesh.triangles, mesh.edges
    local segments, fixed = mesh.segments, mesh.fixed
    local constraints, pending, segment_queue = {}, {}, {}
    local segment_head, heap = 1, {}
    local stats = { steiner_points = 0, segment_splits = 0, circumcenters = 0 }

    -- Short-edge priority is deterministic and avoids repeatedly refining large
    -- triangles while a much smaller geometric feature still needs resolution.
    local function less(a, b)
        return a.size < b.size or (a.size == b.size and a.id < b.id)
    end
    local function push(item)
        local i = #heap + 1
        while i > 1 do
            local parent = math.floor(i / 2)
            if not less(item, heap[parent]) then
                break
            end
            heap[i], i = heap[parent], parent
        end
        heap[i] = item
    end
    local function pop()
        local first, last = heap[1], table.remove(heap)
        if #heap > 0 then
            local i = 1
            while 2 * i <= #heap do
                local child = 2 * i
                if child < #heap and less(heap[child + 1], heap[child]) then
                    child = child + 1
                end
                if not less(heap[child], last) then
                    break
                end
                heap[i], i = heap[child], child
            end
            heap[i] = last
        end
        return first
    end
    local function assess(id)
        local t = triangles[id]
        if not t.alive then
            return
        end
        local area2 = orient(points[t[1]], points[t[2]], points[t[3]])
        local angle, shortest = 180, math.huge
        for j = 1, 3 do
            local a, b, c = points[t[j]], points[t[j % 3 + 1]], points[t[(j + 1) % 3 + 1]]
            local x, y, u, v = b[1] - a[1], b[2] - a[2], c[1] - a[1], c[2] - a[2]
            shortest = math.min(shortest, x * x + y * y)
            angle = math.min(angle, math.deg(math.atan(area2, x * u + y * v)))
        end
        if
            angle + 1e-10 < options.min_angle
            or (options.max_area > 0 and area2 > 2 * options.max_area * (1 + 1e-12))
        then
            push({ id = id, size = shortest })
        end
    end
    local function encroached(s)
        local e = edges[key(s[1], s[2])]
        check(e, "lost constraint during refinement")
        for j = 3, 4 do
            if e[j] then
                local id = opposite(triangles[e[j]], s[1], s[2])
                if predicates.diametral(points[s[1]], points[s[2]], points[id]) < 0 then
                    return true
                end
            end
        end
        return false
    end
    local function enqueue(s)
        if not pending[s] and encroached(s) then
            pending[s] = true
            segment_queue[#segment_queue + 1] = s
        end
    end
    for i, s in ipairs(segments) do
        constraints[key(s[1], s[2])] = { segment = s, index = i }
        enqueue(s)
    end
    for id in ipairs(triangles) do
        assess(id)
    end
    local function vertex(p)
        check(
            stats.steiner_points < options.max_steiner,
            "Steiner point budget exhausted; reduce the quality targets, increase max Steiner points, "
                .. "or simplify acute/narrow geometry"
        )
        check(
            p[1] == p[1]
                and p[2] == p[2]
                and math.abs(p[1]) < math.huge
                and math.abs(p[2]) < math.huge,
            "non-finite circumcentre"
        )
        local count = #points
        local id = mesh.point(p)
        check(
            id == count + 1,
            "refinement reached an existing vertex within mesh tolerance; "
                .. "reduce the quality targets or review the geometry/tolerance"
        )
        stats.steiner_points = stats.steiner_points + 1
        return id
    end
    -- After insertion and local Lawson flips, only newly created triangles can
    -- have acquired a defect or changed the visible apex of a constrained edge.
    local function finish(first)
        local queue = {}
        for id = first, #triangles do
            local t = triangles[id]
            if t.alive then
                for j = 1, 3 do
                    queue[#queue + 1] = key(t[j], t[j % 3 + 1])
                end
            end
        end
        mesh.legalize(queue)
        for id = first, #triangles do
            local t = triangles[id]
            if t.alive then
                assess(id)
                for j = 1, 3 do
                    local record = constraints[key(t[j], t[j % 3 + 1])]
                    if record then
                        enqueue(record.segment)
                    end
                end
            end
        end
    end
    local function split(s)
        local a, b = s[1], s[2]
        local pa, pb = points[a], points[b]
        check(
            (pb[1] - pa[1]) ^ 2 + (pb[2] - pa[2]) ^ 2 > (4 * options.tolerance) ^ 2,
            "segment splitting reached mesh tolerance; reduce the quality targets "
                .. "or review the geometry/tolerance"
        )
        local id = vertex({ pa[1] + (pb[1] - pa[1]) / 2, pa[2] + (pb[2] - pa[2]) / 2 })
        local k = key(a, b)
        local e, record = edges[k], constraints[k]
        local adjacent = { e[3], e[4] }
        local first = #triangles + 1
        for _, tid in ipairs(adjacent) do
            local c = opposite(triangles[tid], a, b)
            mesh.remove(tid)
            mesh.add(a, id, c)
            mesh.add(id, b, c)
        end
        constraints[k], fixed[k] = nil, nil
        local left, right = { a, id, boundary = s.boundary }, { id, b, boundary = s.boundary }
        segments[record.index], segments[#segments + 1] = left, right
        constraints[key(a, id)] = { segment = left, index = record.index }
        constraints[key(id, b)] = { segment = right, index = #segments }
        fixed[key(a, id)], fixed[key(id, b)] = true, true
        stats.segment_splits = stats.segment_splits + 1
        finish(first)
    end
    local function locate(p, at)
        for _ = 1, #triangles do
            local t = triangles[at]
            local nextid
            for j = 1, 3 do
                local a, b = t[j], t[j % 3 + 1]
                if orient(points[a], points[b], p) < 0 then
                    local k = key(a, b)
                    if constraints[k] then
                        return nil, constraints[k].segment
                    end
                    local e = edges[k]
                    nextid = e[3] == at and e[4] or e[3]
                    check(nextid, "circumcentre escaped the protected domain")
                    break
                end
            end
            if not nextid then
                return at
            end
            at = nextid
        end
        error("luafemm mesh quality: circumcentre point location failed")
    end
    local function insert(p, at)
        local t = triangles[at]
        local cavity = { at }
        for j = 1, 3 do
            local a, b = t[j], t[j % 3 + 1]
            if orient(points[a], points[b], p) == 0 then
                local k = key(a, b)
                check(not fixed[k], "attempt to insert on a protected edge")
                local e = edges[k]
                cavity[2] = e[3] == at and e[4] or e[3]
                break
            end
        end
        local border = {}
        for _, tid in ipairs(cavity) do
            local tri = triangles[tid]
            for j = 1, 3 do
                local a, b = tri[j], tri[j % 3 + 1]
                local e = edges[key(a, b)]
                if
                    not (
                        cavity[2]
                        and e[3]
                        and e[4]
                        and (e[3] == cavity[1] or e[3] == cavity[2])
                        and (e[4] == cavity[1] or e[4] == cavity[2])
                    )
                then
                    border[#border + 1] = { a, b }
                end
            end
        end
        local id, first = vertex(p), #triangles + 1
        for _, tid in ipairs(cavity) do
            mesh.remove(tid)
        end
        for _, e in ipairs(border) do
            mesh.add(e[1], e[2], id)
        end
        stats.circumcenters = stats.circumcenters + 1
        finish(first)
    end
    while segment_head <= #segment_queue or #heap > 0 do
        if segment_head <= #segment_queue then
            local s = segment_queue[segment_head]
            segment_head = segment_head + 1
            pending[s] = nil
            local record = constraints[key(s[1], s[2])]
            if record and record.segment == s and encroached(s) then
                split(s)
            end
        else
            local item = pop()
            local t = triangles[item.id]
            if t.alive then
                local a, b, c = points[t[1]], points[t[2]], points[t[3]]
                local x, y, u, v = b[1] - a[1], b[2] - a[2], c[1] - a[1], c[2] - a[2]
                local bb, cc = x * x + y * y, u * u + v * v
                local denominator = 2 * orient(a, b, c)
                local p = {
                    a[1] + (bb * v - cc * y) / denominator,
                    a[2] + (cc * x - bb * u) / denominator,
                }
                local blocking
                -- Conservative global disk check; visibility is not required.
                -- This can split more segments than Triangle's visible test.
                for _, s in ipairs(segments) do
                    if predicates.diametral(points[s[1]], points[s[2]], p) < 0 then
                        blocking = s
                        break
                    end
                end
                local at
                if not blocking then
                    at, blocking = locate(p, item.id)
                end
                if blocking then
                    push(item) -- It may survive the segment split; dead entries are discarded.
                    split(blocking)
                else
                    insert(p, at)
                end
            end
        end
    end
    return stats
end
return R
