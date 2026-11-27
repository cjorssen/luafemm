-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Evaluated-field views for singular sources with a raw P1 correction.
-- RK stages remain on one analytic-plus-affine element extension. Accepted
-- chords stop on the first mesh face; the next cell supplies its own one-sided
-- direction. Contours subdivide actual dielectric triangles, evaluate the
-- nonlinear potential, and solve edge roots. Neither view interpolates a
-- singular nodal potential or substitutes one constant total triangle field.
local V = {}
local H = require("luafemm-electric-hybrid")
local F = require("luafemm-electric-fe")
local E = require("luafemm-electrostatic")
local min, max, sqrt, abs = math.min, math.max, math.sqrt, math.abs
local function check(ok, message)
    assert(ok, "luafemm electric hybrid view: " .. message)
end
local function distance(a, b)
    return sqrt((a[1] - b[1]) ^ 2 + (a[2] - b[2]) ^ 2)
end
local function at(a, b, t)
    return { a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]) }
end
local function limits(m, o)
    local b = m.fe.bounds
    return {
        max(b[1], o.bounds[1]),
        max(b[2], o.bounds[2]),
        min(b[3], o.bounds[3]),
        min(b[4], o.bounds[4]),
    }
end
local function inside(b, p)
    return p[1] >= b[1] and p[1] <= b[3] and p[2] >= b[2] and p[2] <= b[4]
end
local function normal(m, o, p, id, sign)
    if E.clearance(m, p[1], p[2]) <= 0 then
        return nil, "source"
    end
    local s, why = H.sample(m, p[1], p[2], id)
    if not s then
        return nil, why
    end
    if s.electrode then
        return nil, "electrode"
    end
    if s.E <= o.null_field then
        return nil, "null field"
    end
    return { sign * s.Ex / s.E, sign * s.Ey / s.E }
end
local function rk(m, o, p, id, sign, h)
    local a = normal(m, o, p, id, sign)
    if not a then
        return
    end
    local b = normal(m, o, { p[1] + h * a[1] / 2, p[2] + h * a[2] / 2 }, id, sign)
    if not b then
        return
    end
    local c = normal(m, o, { p[1] + h * b[1] / 2, p[2] + h * b[2] / 2 }, id, sign)
    if not c then
        return
    end
    local d = normal(m, o, { p[1] + h * c[1], p[2] + h * c[2] }, id, sign)
    if not d then
        return
    end
    return {
        p[1] + h * (a[1] + 2 * b[1] + 2 * c[1] + d[1]) / 6,
        p[2] + h * (a[2] + 2 * b[2] + 2 * c[2] + d[2]) / 6,
    }
end
local function entering(m, o, p, id, sign)
    local t = m.fe.triangles[id]
    local w = F.barycentric(m, t, p[1], p[2])
    if min(w[1], w[2], w[3]) < -1e-9 then
        return
    end
    local n = normal(m, o, p, id, sign)
    if not n then
        return
    end
    for j = 1, 3 do
        if w[j] <= 1e-9 and t.gx[j] * n[1] + t.gy[j] * n[2] < -1e-9 / m.fe.options.h then
            return
        end
    end
    return n, w
end
local function trace(m, o, seed, sign)
    local bounds, points = limits(m, o), { seed }
    if not inside(bounds, seed) then
        return points, "outside"
    end
    local p, current, length, h = seed, F.locate(m, seed[1], seed[2]), 0, o.step
    if not current then
        return points, "domain"
    end
    local face_tol = max(1e-12, m.fe.options.h * 1e-10)
    for _ = 1, o.max_steps do
        if require("luafemm-electric-fe-view").on_electrode(m, p) then
            return points, "electrode"
        end
        if E.clearance(m, p[1], p[2]) < o.tolerance then
            return points, "source"
        end
        local n, w = entering(m, o, p, current, sign)
        if not n then
            for _, id in
                ipairs(
                    F.index(m)
                        .query(p[1] - face_tol, p[2] - face_tol, p[1] + face_tol, p[2] + face_tol)
                )
            do
                n, w = entering(m, o, p, id, sign)
                if n then
                    current = id
                    break
                end
            end
        end
        if not n then
            return points, "interface stagnation or domain"
        end
        local t = m.fe.triangles[current]
        local exit = math.huge
        for j = 1, 3 do
            local d = t.gx[j] * n[1] + t.gy[j] * n[2]
            if d < -1e-10 / m.fe.options.h then
                exit = min(exit, max(0, w[j]) / -d)
            end
        end
        h = min(
            h,
            o.step,
            o.max_length - length,
            max(o.tolerance, 1.2 * exit),
            max(o.tolerance, 0.45 * E.clearance(m, p[1], p[2]))
        )
        if h < o.tolerance * 1e-5 then
            return points, "length"
        end
        local z, err, accepted
        for _ = 1, 40 do
            local full = rk(m, o, p, current, sign, h)
            local mid = rk(m, o, p, current, sign, h / 2)
            z = mid and rk(m, o, mid, current, sign, h / 2)
            if full and z then
                err = max(distance(full, z) / 15, distance(mid, at(p, z, 0.5)))
                local wm, wz = F.barycentric(m, t, mid[1], mid[2]), F.barycentric(m, t, z[1], z[2])
                local hidden_exit = min(wm[1], wm[2], wm[3]) < -1e-10
                    and min(wz[1], wz[2], wz[3]) >= 0
                if err <= o.tolerance and not hidden_exit then
                    accepted = true
                    break
                end
            end
            h = h / 2
            if h < o.tolerance * 1e-5 then
                return points, "critical point"
            end
        end
        check(accepted, "streamline tolerance could not be attained")
        local wz = F.barycentric(m, t, z[1], z[2])
        local fraction, stop = 1, nil
        for j = 1, 3 do
            if wz[j] < 0 and wz[j] < w[j] then
                fraction = min(fraction, max(0, w[j]) / (w[j] - wz[j]))
            end
        end
        for j = 1, 2 do
            local d = z[j] - p[j]
            if d ~= 0 then
                local a = ((d > 0 and bounds[j + 2] or bounds[j]) - p[j]) / d
                if a >= 0 and a <= fraction then
                    fraction, stop = a, "window"
                end
            end
        end
        z = at(p, z, fraction)
        local ds = distance(p, z)
        if ds < face_tol then
            return points, stop or "interface stagnation or domain"
        end
        points[#points + 1], length, p = z, length + ds, z
        if stop then
            return points, stop
        end
        if length >= o.max_length - o.tolerance * 1e-5 then
            return points, "length"
        end
        h = min(o.step, h * min(2, 0.9 * sqrt(o.tolerance / max(err, 1e-30))))
    end
    return points, "steps"
end
function V.lines(m, o)
    H.solve(m)
    local seeds = o.seeds or {}
    if not o.seeds then
        local sign = o.seed_sign
        if sign == "auto" then
            sign = "negative"
            for _, s in ipairs(m.sources) do
                if s.q > 0 then
                    sign = "positive"
                    break
                end
            end
        end
        for _, s in ipairs(m.sources) do
            if
                sign == "all"
                or sign == "positive" and s.q > 0
                or sign == "negative" and s.q < 0
            then
                for j = 1, o.lines do
                    local a = m.fe.axisymmetric and (-0.5 * math.pi + math.pi * (j - 0.5) / o.lines)
                        or 2 * math.pi * (j - 0.5) / o.lines
                    seeds[#seeds + 1] =
                        { s.x + 1.05 * s.radius * math.cos(a), s.y + 1.05 * s.radius * math.sin(a) }
                end
            end
        end
    end
    local result = {}
    for _, seed in ipairs(seeds) do
        local back, forward, bs, fs = { seed }, { seed }, "disabled", "disabled"
        if o.direction ~= "forward" then
            back, bs = trace(m, o, seed, -1)
        end
        if o.direction ~= "backward" then
            forward, fs = trace(m, o, seed, 1)
        end
        local points = {}
        for j = #back, 2, -1 do
            points[#points + 1] = back[j]
        end
        for _, p in ipairs(forward) do
            points[#points + 1] = p
        end
        if #points > 1 then
            result[#result + 1] = { points = points, start = bs, finish = fs }
        end
    end
    return result
end
local function segment_distance(s, a, b)
    local dx, dy = b[1] - a[1], b[2] - a[2]
    local t = max(0, min(1, ((s.x - a[1]) * dx + (s.y - a[2]) * dy) / (dx * dx + dy * dy)))
    return sqrt((s.x - a[1] - t * dx) ^ 2 + (s.y - a[2] - t * dy) ^ 2)
end
-- A contour cell must not straddle an exclusion disk, even if all vertices
-- happen to lie outside it. Full interior cells can be discarded immediately.
local function mask(m, points)
    for _, s in ipairs(m.sources) do
        local all = true
        for _, p in ipairs(points) do
            all = all and distance(p, { s.x, s.y }) <= s.radius
        end
        if all then
            return true, true
        end
        local positive, negative = false, false
        for j = 1, 3 do
            local a, b = points[j], points[j % 3 + 1]
            if segment_distance(s, a, b) <= s.radius then
                return true, false
            end
            local side = (b[1] - a[1]) * (s.y - a[2]) - (b[2] - a[2]) * (s.x - a[1])
            positive, negative = positive or side > 0, negative or side < 0
        end
        if not (positive and negative) then
            return true, false
        end
    end
    return false, false
end
local function clip(a, b, box)
    local lo, hi = 0, 1
    for j = 1, 2 do
        local d = b[j] - a[j]
        if d == 0 then
            if a[j] < box[j] or a[j] > box[j + 2] then
                return
            end
        else
            local u, v = (box[j] - a[j]) / d, (box[j + 2] - a[j]) / d
            lo, hi = max(lo, min(u, v)), min(hi, max(u, v))
        end
    end
    if lo < hi then
        return at(a, b, lo), at(a, b, hi)
    end
end
-- Join roots geometrically across adaptive hanging edges. Nearby buckets avoid
-- rounding-bin seams; tolerance is many orders below the drawing tolerance.
local function join(segments, level, tolerance)
    local points, adj, bins = {}, {}, {}
    local function vertex(p)
        local x, y = math.floor(p[1] / tolerance), math.floor(p[2] / tolerance)
        for i = x - 1, x + 1 do
            for j = y - 1, y + 1 do
                for _, id in ipairs(bins[i .. ":" .. j] or {}) do
                    if distance(p, points[id]) < tolerance then
                        return id
                    end
                end
            end
        end
        local id = #points + 1
        points[id], adj[id] = p, {}
        local k = x .. ":" .. y
        bins[k] = bins[k] or {}
        bins[k][#bins[k] + 1] = id
        return id
    end
    local edges = {}
    for _, s in ipairs(segments) do
        local a, b = vertex(s[1]), vertex(s[2])
        if a ~= b then
            local id = #edges + 1
            edges[id] = { a, b }
            adj[a][#adj[a] + 1] = id
            adj[b][#adj[b] + 1] = id
        end
    end
    local used, out = {}, {}
    local function walk(start)
        local p, id = {}, start
        while true do
            p[#p + 1] = points[id]
            local nextedge
            for _, e in ipairs(adj[id]) do
                if not used[e] then
                    nextedge = e
                    break
                end
            end
            if not nextedge then
                break
            end
            used[nextedge] = true
            local e = edges[nextedge]
            id = e[1] == id and e[2] or e[1]
            if id == start then
                p[#p + 1] = points[id]
                break
            end
        end
        if #p > 1 then
            out[#out + 1] = { points = p, level = level }
        end
    end
    for id, list in ipairs(adj) do
        if #list ~= 2 then
            walk(id)
        end
    end
    for id, e in ipairs(edges) do
        if not used[id] then
            walk(e[1])
        end
    end
    return out
end
function V.contours(m, o)
    H.solve(m)
    local f, bounds = m.fe, limits(m, o)
    if bounds[1] >= bounds[3] or bounds[2] >= bounds[4] then
        return {}
    end
    local span = max(bounds[3] - bounds[1], bounds[4] - bounds[2])
    local target = (bounds[3] - bounds[1]) / o.grid
    local low, high = math.huge, -math.huge
    for ti, t in ipairs(f.triangles) do
        if not t.electrode then
            local a, b, c = f.nodes[t.ids[1]], f.nodes[t.ids[2]], f.nodes[t.ids[3]]
            -- Interior probes also see curvature when boundary vertices have
            -- equal prescribed potentials. This scale controls drawing only.
            for _, p in ipairs({
                a,
                b,
                c,
                at(a, b, 0.5),
                at(b, c, 0.5),
                at(c, a, 0.5),
                { (a[1] + b[1] + c[1]) / 3, (a[2] + b[2] + c[2]) / 3 },
            }) do
                if inside(bounds, p) and E.clearance(m, p[1], p[2]) > 0 then
                    local s = H.sample(m, p[1], p[2], ti)
                    if s then
                        low, high = min(low, s.V), max(high, s.V)
                    end
                end
            end
        end
    end
    if not (low < high) then
        return {}
    end
    local levels = o.levels or {}
    if not o.levels then
        for j = 1, o.count do
            levels[j] = low + (high - low) * j / (o.count + 1)
        end
    end
    local scale = max(high - low, 1e-300)
    local cells, visited = {}, 0
    local function subdivide(ti, a, b, c, depth)
        visited = visited + 1
        check(visited <= o.contour_cells, "hybrid contour cell budget exceeded")
        if
            max(a[1], b[1], c[1]) < bounds[1]
            or min(a[1], b[1], c[1]) > bounds[3]
            or max(a[2], b[2], c[2]) < bounds[2]
            or min(a[2], b[2], c[2]) > bounds[4]
        then
            return
        end
        local masked, full = mask(m, { a, b, c })
        if full then
            return
        end
        local ab, bc, ca = at(a, b, 0.5), at(b, c, 0.5), at(c, a, 0.5)
        local size = max(distance(a, b), distance(b, c), distance(c, a))
        local values, nonlinear = {}, 0
        if not masked then
            for j, p in ipairs({
                a,
                b,
                c,
                ab,
                bc,
                ca,
                { (a[1] + b[1] + c[1]) / 3, (a[2] + b[2] + c[2]) / 3 },
            }) do
                values[j] = assert(H.sample(m, p[1], p[2], ti)).V
            end
            nonlinear = max(
                abs(values[4] - (values[1] + values[2]) / 2),
                abs(values[5] - (values[2] + values[3]) / 2),
                abs(values[6] - (values[3] + values[1]) / 2),
                abs(values[7] - (values[1] + values[2] + values[3]) / 3)
            ) / scale
        end
        if size > target or nonlinear > o.contour_tolerance then
            check(
                depth < o.contour_depth,
                "hybrid contour tolerance not reached within maximum depth"
            )
            subdivide(ti, a, ab, ca, depth + 1)
            subdivide(ti, ab, b, bc, depth + 1)
            subdivide(ti, ca, bc, c, depth + 1)
            subdivide(ti, ab, bc, ca, depth + 1)
        elseif not masked then
            cells[#cells + 1] = { ti = ti, points = { a, b, c }, values = values }
        end
    end
    for ti, t in ipairs(f.triangles) do
        if not t.electrode then
            subdivide(ti, f.nodes[t.ids[1]], f.nodes[t.ids[2]], f.nodes[t.ids[3]], 0)
        end
    end
    local result = {}
    for _, level in ipairs(levels) do
        local segments = {}
        for _, cell in ipairs(cells) do
            local hits = {}
            for j = 1, 3 do
                local k = j % 3 + 1
                local va, vb = cell.values[j], cell.values[k]
                if (va >= level) ~= (vb >= level) then
                    local a, b = cell.points[j], cell.points[k]
                    local lo, hi = 0, 1
                    for _ = 1, 36 do
                        local mid = (lo + hi) / 2
                        local p = at(a, b, mid)
                        local v = assert(H.sample(m, p[1], p[2], cell.ti)).V
                        if (v >= level) == (va >= level) then
                            lo = mid
                        else
                            hi = mid
                        end
                    end
                    hits[#hits + 1] = at(a, b, (lo + hi) / 2)
                end
            end
            if #hits == 2 then
                local a, b = clip(hits[1], hits[2], bounds)
                if a then
                    segments[#segments + 1] = { a, b }
                end
            end
        end
        for _, path in ipairs(join(segments, level, span * 1e-9)) do
            result[#result + 1] = path
        end
    end
    return result
end
return V
