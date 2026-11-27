-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Bounded vector visibility for projected polylines and opaque convex planes.
-- All three coordinates use drawing centimetres here; larger depth is nearer.
-- BSP splits surfaces, analytical interval clipping splits lines. No field
-- provider is called and no average-depth object ordering is used. Labels and
-- axes are overlays. Line crossings use a schematic gap, not a conductor radius.
local M = {}
local copy = require("luafemm-physical").copy
local function check(ok, message)
    assert(ok, "luafemm 3D visibility: " .. message)
end
local function lerp(a, b, t)
    return { a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]), a[3] + t * (b[3] - a[3]) }
end
local function cross(a, b, c)
    return (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
end
local function distance(a, b)
    return math.sqrt((b[1] - a[1]) ^ 2 + (b[2] - a[2]) ^ 2)
end
local function equation(points)
    local a, b, c = points[1], points[2], points[3]
    local u, v = {}, {}
    for i = 1, 3 do
        u[i], v[i] = b[i] - a[i], c[i] - a[i]
    end
    local n = { u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1] }
    local length = math.sqrt(n[1] ^ 2 + n[2] ^ 2 + n[3] ^ 2)
    if length == 0 then
        return nil
    end
    local sign = n[3] < 0 and -1 or 1
    for i = 1, 3 do
        n[i] = n[i] * sign / length
    end
    return { n[1], n[2], n[3], -(n[1] * a[1] + n[2] * a[2] + n[3] * a[3]) }
end
local function signed(p, e)
    return p[1] * e[1] + p[2] * e[2] + p[3] * e[3] + e[4]
end
local function union(intervals, tolerance)
    table.sort(intervals, function(a, b)
        return a[1] < b[1]
    end)
    local out = {}
    for _, v in ipairs(intervals) do
        local last = out[#out]
        if last and v[1] <= last[2] + tolerance then
            last[2] = math.max(last[2], v[2])
        else
            out[#out + 1] = { v[1], v[2] }
        end
    end
    return out
end
-- Disjoint interval ends are monotone after union. Binary searches keep
-- crossing and mark queries independent of the total number of remote cuts.
local function first_interval(intervals, value)
    local lo, hi = 1, #intervals + 1
    while lo < hi do
        local mid = math.floor((lo + hi) / 2)
        if intervals[mid][2] < value then
            lo = mid + 1
        else
            hi = mid
        end
    end
    return lo
end
local function contains(intervals, s, tolerance)
    local v = intervals[first_interval(intervals, s - tolerance)]
    return v ~= nil and v[1] <= s + tolerance
end

--- The projected convex polygon and linear depth give an exact hidden interval
-- on a straight segment, up to the explicitly chosen geometric tolerance.
local function behind(a, b, face, eps)
    local lo, hi = 0, 1
    local points = face.points
    local orientation = cross(points[1], points[2], points[3]) >= 0 and 1 or -1
    for i, p in ipairs(points) do
        local q = points[i % #points + 1]
        local fa, fb = orientation * cross(p, q, a), orientation * cross(p, q, b)
        local edge = distance(p, q)
        local limit = eps * edge
        if fa < -limit and fb < -limit then
            return nil
        end
        if fa < 0 and fb >= 0 then
            lo = math.max(lo, fa / (fa - fb))
        elseif fb < 0 and fa >= 0 then
            hi = math.min(hi, fa / (fa - fb))
        end
    end
    -- Normal faces toward the viewer. Coplanar lines remain visible.
    local da, db = signed(a, face.equation) + eps, signed(b, face.equation) + eps
    if da >= 0 and db >= 0 then
        return nil
    end
    if da >= 0 then
        lo = math.max(lo, da / (da - db))
    elseif db >= 0 then
        hi = math.min(hi, da / (da - db))
    end
    if hi > lo then
        return lo, hi
    end
end

--- Painter traversal of a BSP of finite planar polygons. Splitting polygons
-- preserves each original face's style; outlines are emitted separately so
-- artificial BSP edges are never stroked.
local function surfaces(faces, eps, work, emit, budget)
    local generated = #faces
    local function split(face, e)
        local front, back = {}, {}
        local points = face.points
        for i, a in ipairs(points) do
            local b = points[i % #points + 1]
            local da, db = signed(a, e), signed(b, e)
            if da >= -eps then
                front[#front + 1] = a
            end
            if da <= eps then
                back[#back + 1] = a
            end
            if (da > eps and db < -eps) or (da < -eps and db > eps) then
                local p = lerp(a, b, da / (da - db))
                front[#front + 1] = p
                back[#back + 1] = p
            end
        end
        return front, back
    end
    local function visit(list, level)
        if #list == 0 then
            return
        end
        check(level <= 128, "surface partition depth exceeded")
        local root = list[1]
        local front, back, same = {}, {}, { root }
        for i = 2, #list do
            work()
            local p = list[i]
            local positive, negative = false, false
            for _, q in ipairs(p.points) do
                local d = signed(q, root.equation)
                positive = positive or d > eps
                negative = negative or d < -eps
            end
            if positive and negative then
                generated = generated + 2
                check(generated <= budget, "surface fragment budget exceeded")
                local a, b = split(p, root.equation)
                local f, g = copy(p), copy(p)
                f.points = a
                g.points = b
                if #a >= 3 then
                    front[#front + 1] = f
                end
                if #b >= 3 then
                    back[#back + 1] = g
                end
            elseif positive then
                front[#front + 1] = p
            elseif negative then
                back[#back + 1] = p
            else
                same[#same + 1] = p
            end
        end
        visit(back, level + 1)
        for _, p in ipairs(same) do
            local q = copy(p)
            q.kind = "path"
            q.closed = true
            q.style = q.style
                .. ",luafemm@fragment,fill,opacity=1,fill opacity="
                .. (q.occlusion == "translucent" and tostring(q.alpha) or "1")
                .. ",draw=none"
            emit(q)
        end
        visit(front, level + 1)
    end
    visit(faces, 1)
end

function M.render(layout)
    local o = layout.options
    local work_count, fragment_count = 0, 0
    local function work()
        work_count = work_count + 1
        check(
            work_count <= o.visibility_budget,
            "work budget exceeded; simplify the scene or raise visibility budget"
        )
    end
    local segment_count = 0
    for _, p in ipairs(layout.primitives) do
        if p.role and p.points then
            segment_count = segment_count + #p.points - 1
        end
    end
    check(
        segment_count <= o.fragment_budget,
        "input segment budget exceeded; raise fragment budget"
    )
    local result = {}
    local function emit(p)
        fragment_count = fragment_count + 1
        check(fragment_count <= o.fragment_budget, "fragment budget exceeded")
        result[#result + 1] = p
    end
    local faces, paths, overlays, background, translucent = {}, {}, {}, {}, {}
    local extent = 0
    -- Convert retained model-unit depth to the same units as screen coordinates.
    local primitives = copy(layout.primitives)
    local low, high = { math.huge, math.huge, math.huge }, { -math.huge, -math.huge, -math.huge }
    for _, p in ipairs(primitives) do
        for _, q in ipairs(p.points or {}) do
            q[3] = q[3] and q[3] * o.coordinate_scale or 0
            check(math.abs(q[3]) < 1e100, "projected depth outside supported numerical range")
            if p.role then
                for i = 1, 3 do
                    low[i] = math.min(low[i], q[i])
                    high[i] = math.max(high[i], q[i])
                end
            end
        end
    end
    for i = 1, 3 do
        extent = math.max(extent, high[i] - low[i])
    end
    local eps = o.visibility_tolerance * (extent > 0 and extent or 1)
    for _, p in ipairs(primitives) do
        if p.role == "surface" and p.occlusion == "opaque" then
            p.equation = equation(p.points)
            -- An edge-on surface has no projected interior and cannot occlude.
            if p.equation and math.abs(p.equation[3]) > o.visibility_tolerance then
                faces[#faces + 1] = p
            end
            if not p.tube then
                local edge = copy(p)
                edge.closed = nil
                edge.role = "plane"
                edge.points[#edge.points + 1] = copy(edge.points[1])
                edge.style = edge.style .. ",fill=none,opacity=1,draw opacity=1"
                paths[#paths + 1] = edge
            end
        elseif p.role == "surface" and p.occlusion == "translucent" then
            p.equation = equation(p.points)
            if p.equation and math.abs(p.equation[3]) > o.visibility_tolerance then
                translucent[#translucent + 1] = p
            end
        elseif p.role == "surface" or p.role == "plane" and p.occlusion == "none" then
            background[#background + 1] = p
        elseif p.role then
            paths[#paths + 1] = p
        else
            overlays[#overlays + 1] = p
        end
    end
    -- Bounded alpha composition: every translucent plane must be strictly
    -- in front of all opaque geometry and curves. Mutual plane intersections
    -- are split by the same BSP; arbitrary interleaved alpha is rejected.
    local front = -math.huge
    for _, list in ipairs({ faces, paths }) do
        for _, p in ipairs(list) do
            for _, q in ipairs(p.points) do
                front = math.max(front, q[3])
            end
        end
    end
    for _, p in ipairs(translucent) do
        for _, q in ipairs(p.points) do
            check(
                q[3] > front + eps,
                "translucent plane must be wholly in front of opaque geometry and curves"
            )
        end
    end
    -- Non-occluding annotations are deliberately background layers. Opaque
    -- surfaces get exact BSP painter ordering, independent of declared layers.
    for _, p in ipairs(background) do
        emit(p)
    end
    surfaces(faces, eps, work, emit, o.fragment_budget)
    local segments = {}
    for id, path in ipairs(paths) do
        path.id = id
        path.intervals = {}
        path.segments = {}
        path.length = 0
        for i = 2, #path.points do
            local a, b = path.points[i - 1], path.points[i]
            local length = distance(a, b)
            if length > eps then
                local s = {
                    a = a,
                    b = b,
                    length = length,
                    start = path.length,
                    path = path,
                    index = i - 1,
                    xmin = math.min(a[1], b[1]),
                    xmax = math.max(a[1], b[1]),
                    ymin = math.min(a[2], b[2]),
                    ymax = math.max(a[2], b[2]),
                }
                path.segments[#path.segments + 1] = s
                segments[#segments + 1] = s
                for _, face in ipairs(faces) do
                    work()
                    if path.plane_id ~= face.plane_id then
                        local lo, hi = behind(a, b, face, eps)
                        if lo then
                            path.intervals[#path.intervals + 1] =
                                { s.start + lo * length, s.start + hi * length }
                        end
                    end
                end
                path.length = path.length + length
            end
        end
        path.intervals = union(path.intervals, eps)
        path.plane_intervals = copy(path.intervals)
    end
    local function hidden(s, lo, hi)
        s.path.intervals[#s.path.intervals + 1] = {
            math.max(0, s.start + lo * s.length),
            math.min(s.path.length, s.start + hi * s.length),
        }
    end
    -- Sweep projected x intervals; bounded candidate comparisons avoid an
    -- unbounded quadratic pass on a dense view. Exact crossings interpolate
    -- both depths. Collinear overlaps split where their depth order changes.
    if o.crossing_gap > 0 then
        table.sort(segments, function(a, b)
            if a.xmin ~= b.xmin then
                return a.xmin < b.xmin
            end
            if a.path.id ~= b.path.id then
                return a.path.id < b.path.id
            end
            return a.index < b.index
        end)
        for i, s in ipairs(segments) do
            if s.path.role ~= "plane" then
                for j = i + 1, #segments do
                    local t = segments[j]
                    if t.xmin > s.xmax + eps then
                        break
                    end
                    work()
                    if
                        t.path.role ~= "plane"
                        and t.ymin <= s.ymax + eps
                        and t.ymax >= s.ymin - eps
                        and not (s.path == t.path and math.abs(s.index - t.index) <= 1)
                    then
                        local ax, ay = s.b[1] - s.a[1], s.b[2] - s.a[2]
                        local bx, by = t.b[1] - t.a[1], t.b[2] - t.a[2]
                        local dx, dy = t.a[1] - s.a[1], t.a[2] - s.a[2]
                        local det = ax * by - ay * bx
                        if math.abs(det) > o.visibility_tolerance * s.length * t.length then
                            local u, v = (dx * by - dy * bx) / det, (dx * ay - dy * ax) / det
                            if u >= 0 and u <= 1 and v >= 0 and v <= 1 then
                                local za, zb = lerp(s.a, s.b, u)[3], lerp(t.a, t.b, v)[3]
                                if
                                    not contains(
                                        s.path.plane_intervals,
                                        s.start + u * s.length,
                                        eps
                                    )
                                    and not contains(
                                        t.path.plane_intervals,
                                        t.start + v * t.length,
                                        eps
                                    )
                                then
                                    if za < zb - eps then
                                        hidden(
                                            s,
                                            u - o.crossing_gap / s.length,
                                            u + o.crossing_gap / s.length
                                        )
                                    elseif zb < za - eps then
                                        hidden(
                                            t,
                                            v - o.crossing_gap / t.length,
                                            v + o.crossing_gap / t.length
                                        )
                                    end
                                end
                            end
                        elseif math.abs(dx * ay - dy * ax) <= eps * s.length then
                            local den = s.length ^ 2
                            local u0 = (dx * ax + dy * ay) / den
                            local u1 = u0 + (bx * ax + by * ay) / den
                            local lo, hi =
                                math.max(0, math.min(u0, u1)), math.min(1, math.max(u0, u1))
                            if hi > lo and math.abs(u1 - u0) > o.visibility_tolerance then
                                local function other(u)
                                    return (u - u0) / (u1 - u0)
                                end
                                local function delta(u)
                                    return lerp(s.a, s.b, u)[3] - lerp(t.a, t.b, other(u))[3]
                                end
                                local cuts = { lo, hi }
                                local dl, dh = delta(lo), delta(hi)
                                if dl * dh < 0 then
                                    cuts[#cuts + 1] = lo + (hi - lo) * dl / (dl - dh)
                                end
                                -- Plane visibility can change within a collinear overlap.
                                -- Split at those boundaries before testing either line.
                                local function boundaries(seg, map)
                                    local ranges = seg.path.plane_intervals
                                    local index = first_interval(ranges, seg.start - eps)
                                    while
                                        ranges[index]
                                        and ranges[index][1] <= seg.start + seg.length + eps
                                    do
                                        work()
                                        for _, v in ipairs(ranges[index]) do
                                            local u = map(v)
                                            if u > lo and u < hi then
                                                cuts[#cuts + 1] = u
                                            end
                                        end
                                        index = index + 1
                                    end
                                end
                                boundaries(s, function(v)
                                    return (v - s.start) / s.length
                                end)
                                boundaries(t, function(v)
                                    return u0 + (v - t.start) / t.length * (u1 - u0)
                                end)
                                table.sort(cuts)
                                for k = 2, #cuts do
                                    local a, b = cuts[k - 1], cuts[k]
                                    local mid = (a + b) / 2
                                    local v = other(mid)
                                    if
                                        not contains(
                                            s.path.plane_intervals,
                                            s.start + mid * s.length,
                                            eps
                                        )
                                        and not contains(
                                            t.path.plane_intervals,
                                            t.start + v * t.length,
                                            eps
                                        )
                                    then
                                        if delta(mid) < -eps then
                                            hidden(s, a, b)
                                        elseif delta(mid) > eps then
                                            hidden(
                                                t,
                                                math.min(other(a), other(b)),
                                                math.max(other(a), other(b))
                                            )
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    local hidden_parts, visible_parts, marks = {}, {}, {}
    local pending = 0
    local function reserve()
        pending = pending + 1
        check(fragment_count + pending <= o.fragment_budget, "fragment budget exceeded")
    end
    local function segment_at(path, s)
        local lo, hi = 1, #path.segments + 1
        while lo < hi do
            local mid = math.floor((lo + hi) / 2)
            local seg = path.segments[mid]
            if seg.start + seg.length + eps < s then
                lo = mid + 1
            else
                hi = mid
            end
        end
        return lo
    end
    local function at(path, s)
        local seg = path.segments[segment_at(path, s)]
        if seg then
            local t = math.max(0, math.min(1, (s - seg.start) / seg.length))
            return lerp(seg.a, seg.b, t),
                { (seg.b[1] - seg.a[1]) / seg.length, (seg.b[2] - seg.a[2]) / seg.length }
        end
    end
    local function piece(path, lo, hi, is_hidden)
        if hi - lo <= eps then
            return
        end
        reserve()
        local first = at(path, lo)
        local points = { first }
        for i = segment_at(path, lo), #path.segments do
            local s = path.segments[i]
            local finish = s.start + s.length
            if finish >= hi - eps then
                break
            end
            if finish > lo + eps and finish < hi - eps then
                points[#points + 1] = s.b
            end
        end
        points[#points + 1] = at(path, hi)
        local p = {
            kind = "path",
            points = points,
            clip = path.clip,
            role = path.role,
            source_id = path.id,
            interval = { lo, hi },
            hidden = is_hidden,
            style = path.style .. ",luafemm@fragment",
        }
        if is_hidden then
            p.style = p.style .. ",femm 3d hidden"
            hidden_parts[#hidden_parts + 1] = p
        else
            visible_parts[#visible_parts + 1] = p
        end
    end
    for _, p in ipairs(paths) do
        p.intervals = union(p.intervals, eps)
        local start = 0
        for _, r in ipairs(p.intervals) do
            piece(p, start, r[1], false)
            if o.hidden == "dashed" then
                piece(p, r[1], r[2], true)
            end
            start = r[2]
        end
        piece(p, start, p.length, false)
        -- Generate marks incrementally: budget checks precede allocation,
        -- including when a tiny spacing would otherwise collect huge arrays.
        local curve = require("luafemm-arrows").path(p.points)
        local function mark(s)
            work()
            if not contains(p.intervals, s, o.mark_clearance) then
                local available = s
                local previous = p.intervals[first_interval(p.intervals, s) - 1]
                if previous then
                    available = math.min(available, s - previous[2] - o.mark_clearance)
                end
                local carrier = require("luafemm-arrows").mark(curve, s, {
                    length = o.mark_length,
                    tolerance = o.mark_tolerance,
                    minimum = o.mark_minimum,
                    mode = o.mark_mode,
                }, available)
                if carrier then
                    reserve()
                    marks[#marks + 1] = {
                        kind = "path",
                        points = carrier.points,
                        controls = carrier.controls,
                        clip = p.clip,
                        mark_at = s,
                        source_id = p.id,
                        style = p.style
                            .. ",luafemm@fragment,luafemm@curve tip={"
                            .. carrier.factor
                            .. "},femm 3d mark",
                    }
                end
            end
        end
        if p.role == "field" and o.field_marks == "true" then
            if o.mark_spacing == 0 then
                mark(o.mark_position * p.length)
            else
                check(
                    p.length / o.mark_spacing <= o.fragment_budget,
                    "fragment budget exceeded by requested marks"
                )
                for s = o.mark_position * o.mark_spacing, p.length, o.mark_spacing do
                    mark(s)
                end
            end
        elseif p.role == "vectors" then
            mark(p.length)
        end
    end
    for _, list in ipairs({ hidden_parts, visible_parts, marks }) do
        for _, p in ipairs(list) do
            emit(p)
        end
    end
    surfaces(translucent, eps, work, emit, o.fragment_budget)
    for _, p in ipairs(overlays) do
        emit(p)
    end
    layout.primitives = result
    layout.visibility = { work = work_count, fragments = fragment_count, tolerance = eps }
    return layout
end
return M
