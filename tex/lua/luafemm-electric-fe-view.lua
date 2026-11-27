-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Views and measurements of raw planar electric P1 fields.
-- Lines walk elements exactly; profiles retain one-sided gradients. No recovered
-- field is substituted across a dielectric or an ordinary discretization edge.
-- @module luafemm-electric-fe-view
local V = {}
local F = require("luafemm-electric-fe")
local E = require("luafemm-electrostatic")
local P = require("luafemm-physical")
local min, max, abs, sqrt = math.min, math.max, math.abs, math.sqrt
local function check(v, s)
    assert(v, "luafemm electric FE view: " .. tostring(s))
end
local function distance(a, b)
    return sqrt((a[1] - b[1]) ^ 2 + (a[2] - b[2]) ^ 2)
end
local function at(a, b, t)
    return { a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]) }
end
-- Partition one measurement segment by the actual triangulation. Candidates are
-- conservative boxes; barycentric clipping makes the final geometric decision.
function V.cuts(m, a, b)
    local cuts = { 0, 1 }
    for _, id in
        ipairs(F.index(m).query(min(a[1], b[1]), min(a[2], b[2]), max(a[1], b[1]), max(a[2], b[2])))
    do
        local tri = m.fe.triangles[id]
        local u, w = F.barycentric(m, tri, a[1], a[2]), F.barycentric(m, tri, b[1], b[2])
        local lo, hi = 0, 1
        for j = 1, 3 do
            local d = w[j] - u[j]
            if abs(d) < 1e-15 then
                if u[j] < -1e-10 then
                    lo, hi = 1, 0
                end
            elseif d > 0 then
                lo = max(lo, -u[j] / d)
            else
                hi = min(hi, -u[j] / d)
            end
        end
        if lo <= hi + 1e-12 then
            cuts[#cuts + 1] = max(0, min(1, lo))
            cuts[#cuts + 1] = max(0, min(1, hi))
        end
    end
    if m.hybrid then
        -- Include both disk boundaries and an interior mark so an exclusion
        -- narrower than the profile spacing still produces an explicit gap.
        local dx, dy = b[1] - a[1], b[2] - a[2]
        local aa = dx * dx + dy * dy
        for _, source in ipairs(m.sources) do
            local x, y = a[1] - source.x, a[2] - source.y
            local bb = x * dx + y * dy
            local disc = bb * bb - aa * (x * x + y * y - source.radius * source.radius)
            if aa > 0 and disc >= 0 then
                local lo, hi = max(0, (-bb - sqrt(disc)) / aa), min(1, (-bb + sqrt(disc)) / aa)
                if lo <= hi then
                    cuts[#cuts + 1] = lo
                    cuts[#cuts + 1] = (lo + hi) / 2
                    cuts[#cuts + 1] = hi
                end
            end
        end
    end
    table.sort(cuts)
    local unique = {}
    for _, t in ipairs(cuts) do
        if #unique == 0 or t - unique[#unique] > 1e-12 then
            unique[#unique + 1] = t
        end
    end
    return unique
end
--- Raw profiles have finite potential and explicit E/D gaps between distinct
-- element gradients. samples refines arclength sampling; mesh faces are always
-- represented, including interfaces narrower than the uniform sample spacing.
function V.profile(m, points, o)
    E.freeze(m)
    o = o or {}
    local samples = E.integer(o.samples or 201, "profile samples", 2, 100000)
    local outside = o.outside or "error"
    check(outside == "error" or outside == "nan", "unknown outside policy")
    local lengths = { 0 }
    for i = 2, #points do
        lengths[i] = lengths[i - 1] + distance(points[i - 1], points[i])
    end
    local total = lengths[#points]
    check(total and total > 0, "empty profile")
    local marks = { 0, total }
    for i = 2, #points do
        for _, t in ipairs(V.cuts(m, points[i - 1], points[i])) do
            marks[#marks + 1] = lengths[i - 1] + t * (lengths[i] - lengths[i - 1])
        end
    end
    for i = 1, samples - 2 do
        marks[#marks + 1] = total * i / (samples - 1)
    end
    table.sort(marks)
    local sorted = {}
    for _, s in ipairs(marks) do
        if #sorted == 0 or s - sorted[#sorted] > total * 1e-12 then
            sorted[#sorted + 1] = s
        end
    end
    local function point(s)
        local k = 2
        while k < #points and lengths[k] < s do
            k = k + 1
        end
        local l = lengths[k] - lengths[k - 1]
        return at(points[k - 1], points[k], l > 0 and (s - lengths[k - 1]) / l or 0)
    end
    local sides = {}
    for i = 2, #sorted do
        local p = point((sorted[i - 1] + sorted[i]) / 2)
        sides[i] = F.locate(m, p[1], p[2]) or false
    end
    local data = {}
    local function sample(s, element, gap)
        local p = point(s)
        local f, why
        if element then
            f, why = F.sample(m, p[1], p[2], element)
        else
            why = "outside electric FE domain"
        end
        if m.hybrid and E.clearance(m, p[1], p[2]) <= 0 then
            f, why = nil, "source exclusion"
        end
        check(f or outside == "nan", why)
        f = f or { V = 0 / 0, Ex = 0 / 0, Ey = 0 / 0, E = 0 / 0, Dx = 0 / 0, Dy = 0 / 0, D = 0 / 0 }
        if gap then
            for _, k in ipairs({ "Ex", "Ey", "E", "Dx", "Dy", "D" }) do
                f[k] = 0 / 0
            end
        end
        if m.fe.axisymmetric then
            f.Er, f.Ez, f.Dr, f.Dz = f.Ex, f.Ey, f.Dx, f.Dy
            f.r, f.z = p[1], p[2]
        end
        f.x, f.y, f.s, f.t = p[1], p[2], s, s / total
        data[#data + 1] = f
    end
    for i, s in ipairs(sorted) do
        local left = i > 1 and sides[i] or sides[2]
        local right = i < #sorted and sides[i + 1] or sides[i]
        sample(s, left)
        if left ~= right then
            local a, b = left and m.fe.triangles[left], right and m.fe.triangles[right]
            local different = m.hybrid
                or not a
                or not b
                or abs(a.Ex - b.Ex)
                        + abs(a.Ey - b.Ey)
                        + abs(a.Dx - b.Dx) / E.epsilon0
                        + abs(a.Dy - b.Dy) / E.epsilon0
                    > 1e-9 * max(1e-30, a.E, b.E)
            if different then
                sample(s, left, true)
                sample(s, right)
            end
        end
    end
    return { model = E.snapshot(m), points = P.copy(points), data = data, length = total }
end
-- A zero-thickness voltage path can have fields pointing in the same direction
-- on both sides. A gradient-only walk would then cross it; physical electrode
-- membership must stop the line even when no solid conductor triangle exists.
local function on_electrode(m, point)
    local tolerance = max(8 * m.fe.stats.mesh_tolerance, 1e-10 * m.fe.options.h)
    for _, path in ipairs(m.fe.paths) do
        if path.kind == "electrode" then
            local count = path.closed and #path.points or #path.points - 1
            for i = 1, count do
                local a, b = path.points[i], path.points[i % #path.points + 1]
                local dx, dy = b[1] - a[1], b[2] - a[2]
                local length = sqrt(dx * dx + dy * dy)
                local along = ((point[1] - a[1]) * dx + (point[2] - a[2]) * dy) / length
                if
                    along >= -tolerance
                    and along <= length + tolerance
                    and abs(dx * (point[2] - a[2]) - dy * (point[1] - a[1]))
                        <= tolerance * length
                then
                    return true
                end
            end
        end
    end
    return false
end
V.on_electrode = on_electrode
local function trace(m, o, seed, sign)
    local f = m.fe
    local points = { seed }
    local p = seed
    local length = 0
    if on_electrode(m, seed) then
        return points, "electrode"
    end
    local current = F.locate(m, p[1], p[2])
    local previous
    if not current then
        return points, "domain"
    end
    for _ = 1, o.max_steps do
        local t = f.triangles[current]
        if t.electrode then
            return points, "electrode"
        end
        if t.E <= o.null_field then
            return points, "null field"
        end
        local dx, dy = sign * t.Ex / t.E, sign * t.Ey / t.E
        local w = F.barycentric(m, t, p[1], p[2])
        local step = math.huge
        for j = 1, 3 do
            local d = t.gx[j] * dx + t.gy[j] * dy
            if d < -math.max(1e-10, 100 * f.tolerance) / f.options.h then
                step = min(step, max(0, w[j]) / -d)
            end
        end
        local allowed = o.max_length - length
        local stop = "length"
        for _, v in ipairs({
            { p[1], dx, o.bounds[1], o.bounds[3] },
            {
                p[2],
                dy,
                o.bounds[2],
                o.bounds[4],
            },
        }) do
            if v[2] ~= 0 then
                local a = ((v[2] > 0 and v[4] or v[3]) - v[1]) / v[2]
                if a < allowed then
                    allowed = max(0, a)
                    stop = "window"
                end
            end
        end
        if allowed < step then
            if allowed > 0 then
                points[#points + 1] = { p[1] + allowed * dx, p[2] + allowed * dy }
            end
            return points, stop
        end
        check(step < math.huge, "unbounded triangle walk")
        if step > f.options.h * 1e-13 then
            p = { p[1] + step * dx, p[2] + step * dy }
            points[#points + 1] = p
            length = length + step
        end
        if on_electrode(m, p) then
            return points, "electrode"
        end
        local nextid
        for _, id in
            ipairs(
                F.index(m).query(
                    p[1] - f.options.h * 1e-8,
                    p[2] - f.options.h * 1e-8,
                    p[1] + f.options.h * 1e-8,
                    p[2] + f.options.h * 1e-8
                )
            )
        do
            if id ~= current and id ~= previous then
                local q = f.triangles[id]
                local weights = F.barycentric(m, q, p[1], p[2])
                if min(weights[1], weights[2], weights[3]) >= -1e-9 then
                    if q.electrode then
                        return points, "electrode"
                    end
                    if q.E > o.null_field then
                        local enter = true
                        for j = 1, 3 do
                            if
                                weights[j] <= 1e-9
                                and sign * (q.gx[j] * q.Ex + q.gy[j] * q.Ey) / q.E
                                    < -math.max(1e-10, 100 * f.tolerance) / f.options.h
                            then
                                enter = false
                            end
                        end
                        if enter then
                            nextid = id
                            break
                        end
                    end
                end
            end
        end
        if not nextid then
            return points, "interface stagnation or domain"
        end
        previous, current = current, nextid
    end
    return points, "steps"
end
function V.lines(m, o)
    E.freeze(m)
    local f = m.fe
    local seeds = o.seeds
    if not seeds then
        local positive, negative = {}, {}
        for _, e in ipairs(f.edges) do
            if e.electrode or e.outer and e.condition == "voltage" then
                for _, id in ipairs(e.triangles) do
                    local tri = f.triangles[id]
                    if not tri.electrode and tri.E > o.null_field then
                        local a, b = f.nodes[e[1]], f.nodes[e[2]]
                        local middle = { (a[1] + b[1]) / 2, (a[2] + b[2]) / 2 }
                        local center = { 0, 0 }
                        for _, n in ipairs(tri.ids) do
                            center[1] = center[1] + f.nodes[n][1] / 3
                            center[2] = center[2] + f.nodes[n][2] / 3
                        end
                        local dot = tri.Ex * (center[1] - middle[1])
                            + tri.Ey * (center[2] - middle[2])
                        local bucket = dot > 0 and positive or negative
                        bucket[#bucket + 1] =
                            { a = a, b = b, center = center, length = distance(a, b) }
                    end
                end
            end
        end
        local candidates = {}
        if o.seed_sign == "all" then
            for _, p in ipairs(positive) do
                candidates[#candidates + 1] = p
            end
            for _, p in ipairs(negative) do
                candidates[#candidates + 1] = p
            end
        elseif o.seed_sign == "negative" or o.seed_sign == "auto" and #positive == 0 then
            candidates = negative
        else
            candidates = positive
        end
        seeds = {}
        local sum = 0
        for _, p in ipairs(candidates) do
            sum = sum + p.length
        end
        if sum > 0 then
            for i = 1, o.lines do
                local target = (i - 0.5) * sum / o.lines
                local accumulated = 0
                for _, p in ipairs(candidates) do
                    if accumulated + p.length >= target then
                        seeds[#seeds + 1] =
                            at(at(p.a, p.b, (target - accumulated) / p.length), p.center, 1e-6)
                        break
                    end
                    accumulated = accumulated + p.length
                end
            end
        end
    end
    local result = {}
    for _, seed in ipairs(seeds) do
        if
            seed[1] >= o.bounds[1]
            and seed[1] <= o.bounds[3]
            and seed[2] >= o.bounds[2]
            and seed[2] <= o.bounds[4]
        then
            local backward, forward, bs, fs = { seed }, { seed }, "disabled", "disabled"
            if o.direction ~= "forward" then
                backward, bs = trace(m, o, seed, -1)
            end
            if o.direction ~= "backward" then
                forward, fs = trace(m, o, seed, 1)
            end
            local points = {}
            for j = #backward, 2, -1 do
                points[#points + 1] = backward[j]
            end
            for _, p in ipairs(forward) do
                points[#points + 1] = p
            end
            if #points > 1 then
                result[#result + 1] = { points = points, start = bs, finish = fs }
            end
        end
    end
    return result
end
--- Omit vector arrows crossing physical dielectric/charge/electrode boundaries.
-- Ordinary P1 gradient changes remain raw local vector samples.
function V.crosses_interface(m, a, b)
    -- A segment can cross a charged layer exactly at a mesh vertex, where its
    -- neighboring intervals need not share an edge. Audit the physical layer
    -- directly as well as the element-pair topology below.
    local dx, dy = b[1] - a[1], b[2] - a[2]
    for _, path in ipairs(m.fe.paths) do
        if path.kind == "electrode" or path.kind == "charge" and path.density ~= 0 then
            local count = path.closed and #path.points or #path.points - 1
            for i = 1, count do
                local p, q = path.points[i], path.points[i % #path.points + 1]
                local ux, uy = q[1] - p[1], q[2] - p[2]
                local det = dx * uy - dy * ux
                if det ~= 0 then
                    local t = ((p[1] - a[1]) * uy - (p[2] - a[2]) * ux) / det
                    local u = ((p[1] - a[1]) * dy - (p[2] - a[2]) * dx) / det
                    if t > 1e-12 and t < 1 - 1e-12 and u >= -1e-12 and u <= 1 + 1e-12 then
                        return true
                    end
                end
            end
        end
    end
    local cuts = V.cuts(m, a, b)
    local previous
    for i = 2, #cuts do
        local p = at(a, b, (cuts[i - 1] + cuts[i]) / 2)
        local id = F.locate(m, p[1], p[2])
        local t = id and m.fe.triangles[id]
        if not t or t.electrode then
            return true
        end
        if previous then
            local k = math.min(previous, id) .. ":" .. math.max(previous, id)
            local before, after = m.fe.triangles[previous].tensor, t.tensor
            -- Crossing exactly at a mesh vertex need not visit adjacent faces.
            -- Compare the complete law, including rotated off-diagonal entries.
            if
                m.fe.interface_pairs[k]
                or before.xx ~= after.xx
                or before.xy ~= after.xy
                or before.yy ~= after.yy
            then
                return true
            end
        end
        previous = id
    end
    return false
end
--- Equipotentials of the solved P1 potential, without resampling its gradients.
-- Shared edge/vertex identities join triangles before ordinary TikZ drawing.
function V.contours(m, o)
    E.freeze(m)
    local f = m.fe
    local low, high = math.huge, -math.huge
    for _, v in ipairs(f.solution) do
        low = min(low, v)
        high = max(high, v)
    end
    local levels = o.levels or {}
    if not o.levels and low < high then
        for i = 1, o.count do
            levels[i] = low + (high - low) * i / (o.count + 1)
        end
    end
    local result = {}
    for _, level in ipairs(levels) do
        local vertices, segments, adj = {}, {}, {}
        local function vertex(a, b)
            local va, vb = f.solution[a], f.solution[b]
            local t = (level - va) / (vb - va)
            local k
            if t == 0 then
                k = "n" .. a
            elseif t == 1 then
                k = "n" .. b
            else
                k = "e" .. min(a, b) .. ":" .. max(a, b)
            end
            vertices[k] = vertices[k] or at(f.nodes[a], f.nodes[b], t)
            return k
        end
        for ti, tri in ipairs(f.triangles) do
            local hits = {}
            for j = 1, 3 do
                local a, b = tri.ids[j], tri.ids[j % 3 + 1]
                if (f.solution[a] >= level) ~= (f.solution[b] >= level) then
                    hits[#hits + 1] = vertex(a, b)
                end
            end
            if #hits == 2 and hits[1] ~= hits[2] then
                local a, b = vertices[hits[1]], vertices[hits[2]]
                local lo, hi = 0, 1
                for j = 1, 2 do
                    local d = b[j] - a[j]
                    if d == 0 then
                        if a[j] < o.bounds[j] or a[j] > o.bounds[j + 2] then
                            lo, hi = 1, 0
                        end
                    else
                        local u, w = (o.bounds[j] - a[j]) / d, (o.bounds[j + 2] - a[j]) / d
                        lo = max(lo, min(u, w))
                        hi = min(hi, max(u, w))
                    end
                end
                if lo < hi then
                    if lo > 0 then
                        hits[1] = "c" .. ti .. ":1"
                        vertices[hits[1]] = at(a, b, lo)
                    end
                    if hi < 1 then
                        hits[2] = "c" .. ti .. ":2"
                        vertices[hits[2]] = at(a, b, hi)
                    end
                    segments[#segments + 1] = hits
                    for _, k in ipairs(hits) do
                        adj[k] = adj[k] or {}
                        adj[k][#adj[k] + 1] = #segments
                    end
                end
            end
        end
        local used = {}
        local function walk(start)
            local points, k = {}, start
            while true do
                points[#points + 1] = vertices[k]
                local nextedge
                for _, id in ipairs(adj[k]) do
                    if not used[id] then
                        nextedge = id
                        break
                    end
                end
                if not nextedge then
                    break
                end
                used[nextedge] = true
                local e = segments[nextedge]
                k = e[1] == k and e[2] or e[1]
                if k == start then
                    points[#points + 1] = vertices[k]
                    break
                end
            end
            if #points > 1 then
                result[#result + 1] = { points = points, level = level }
            end
        end
        for _, e in ipairs(segments) do
            for _, k in ipairs(e) do
                if #adj[k] == 1 and not used[adj[k][1]] then
                    walk(k)
                end
            end
        end
        for id, e in ipairs(segments) do
            if not used[id] then
                walk(e[1])
            end
        end
    end
    return result
end
return V
