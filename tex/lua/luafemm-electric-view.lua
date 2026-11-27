-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Electric view geometry, independent of TeX and of display transforms.
-- Streamlines use adaptive RK4 step doubling on the unit E vector. Scalar
-- contours use a conforming triangular sampling grid, never a magnetic A map.
-- All spatial tolerances refer to model coordinates; sampled fields remain SI.
-- @module luafemm-electric-view
local R = {}
local E = require("luafemm-electrostatic")
local P = require("luafemm-physical")
local sqrt, abs, min, max = math.sqrt, math.abs, math.min, math.max
local function check(v, s)
    assert(v, "luafemm electric view: " .. s)
end
function R.options(m, o)
    o = o or {}
    local r = {}
    r.bounds = {
        E.number(o.xmin or m.bounds[1], "xmin"),
        E.number(o.ymin or m.bounds[2], "ymin"),
        E.number(o.xmax or m.bounds[3], "xmax"),
        E.number(o.ymax or m.bounds[4], "ymax"),
    }
    local b = r.bounds
    check(b[1] < b[3] and b[2] < b[4], "empty view window")
    r.arrow_position = E.number(o.arrow_position or 0.45, "arrow position")
    check(
        r.arrow_position >= 0 and r.arrow_position <= 0.95,
        "arrow position must lie between 0 and .95"
    )
    r.step = E.positive(o.step or 0.08, "step")
    r.tolerance = E.positive(o.tolerance or 0.0001, "trace tolerance")
    r.max_length = E.positive(o.max_length or 40, "maximum length")
    r.max_steps = E.integer(o.max_steps or 10000, "maximum steps", 1, 1000000)
    r.null_field = E.positive(o.null_field or 1e-9, "null field")
    r.lines = E.integer(o.lines or 16, "lines per charge", 1, 1000)
    r.grid = E.integer(o.grid or 80, "contour resolution", 4, 1000)
    r.contour_tolerance = E.positive(o.contour_tolerance or 0.002, "contour tolerance")
    r.contour_depth = E.integer(o.contour_depth or 8, "maximum contour depth", 1, 14)
    r.contour_cells = E.integer(o.contour_cells or 250000, "maximum contour cells", 1, 2000000)
    r.count = E.integer(o.count or 11, "level count", 1, 100)
    r.spacing = E.positive(o.spacing or 0.5, "vector spacing")
    r.vector_length = E.positive(o.vector_length or 0.28, "vector length")
    r.vector_scale = E.positive(o.vector_scale or 1e-5, "vector scale")
    r.vector_mode = o.vector_mode or "normalized"
    check(r.vector_mode == "normalized" or r.vector_mode == "scaled", "unknown vector mode")
    r.direction = o.direction or "both"
    check(
        r.direction == "both" or r.direction == "forward" or r.direction == "backward",
        "unknown direction"
    )
    r.seed_sign = o.seed_sign or "auto"
    check(
        r.seed_sign == "auto"
            or r.seed_sign == "positive"
            or r.seed_sign == "negative"
            or r.seed_sign == "all",
        "unknown seed sign"
    )
    r.seeds = o.seeds
    r.levels = o.levels
    if r.seeds then
        check(#r.seeds > 0 and #r.seeds <= 10000, "invalid seed count")
        for _, s in ipairs(r.seeds) do
            E.number(s[1], "seed x")
            E.number(s[2], "seed y")
        end
    end
    if r.levels then
        check(#r.levels > 0 and #r.levels <= 100, "invalid level count")
        local seen = {}
        for _, v in ipairs(r.levels) do
            E.number(v, "potential level")
            check(not seen[v], "duplicate level")
            seen[v] = true
        end
    end
    return r
end
local function inside(b, x, y)
    return x >= b[1] and y >= b[2] and x <= b[3] and y <= b[4]
end
local function normfield(m, o, x, y)
    if E.clearance(m, x, y) <= 0 then
        return nil, "source"
    end
    local f, why = P.sample(m, x, y)
    if not f then
        return nil, why
    end
    if f.E <= o.null_field then
        return nil, "null field"
    end
    return { f.Ex / f.E, f.Ey / f.E }
end
-- First intersection of a segment with a view edge or exclusion circle.
-- This is used for both streamlines and arrows, so neither jumps over a mask.
local function clip(m, b, a, z)
    local dx, dy = z[1] - a[1], z[2] - a[2]
    local t, reason = 1, nil
    for _, v in ipairs({ { dx, b[1] - a[1], b[3] - a[1] }, { dy, b[2] - a[2], b[4] - a[2] } }) do
        if v[1] ~= 0 then
            local q = (v[1] > 0 and v[3] or v[2]) / v[1]
            if q >= 0 and q < t then
                t, reason = q, "window"
            end
        end
    end
    local aa = dx * dx + dy * dy
    if aa > 0 then
        for _, s in ipairs(m.sources) do
            local x, y = a[1] - s.x, a[2] - s.y
            local bb = x * dx + y * dy
            local cc = x * x + y * y - s.radius * s.radius
            local disc = bb * bb - aa * cc
            if disc >= 0 then
                local q = (-bb - sqrt(disc)) / aa
                if q >= 0 and q <= t then
                    t, reason = q, "source"
                end
            end
        end
    end
    return { a[1] + t * dx, a[2] + t * dy }, reason
end
local function rk(m, o, p, h)
    local a, why = normfield(m, o, p[1], p[2])
    if not a then
        return nil, why
    end
    local b
    b, why = normfield(m, o, p[1] + h * a[1] / 2, p[2] + h * a[2] / 2)
    if not b then
        return nil, why
    end
    local c
    c, why = normfield(m, o, p[1] + h * b[1] / 2, p[2] + h * b[2] / 2)
    if not c then
        return nil, why
    end
    local d
    d, why = normfield(m, o, p[1] + h * c[1], p[2] + h * c[2])
    if not d then
        return nil, why
    end
    return {
        p[1] + h * (a[1] + 2 * b[1] + 2 * c[1] + d[1]) / 6,
        p[2] + h * (a[2] + 2 * b[2] + 2 * c[2] + d[2]) / 6,
    }
end
local function halftrace(m, o, seed, sign)
    local points = { { seed[1], seed[2] } }
    if not inside(o.bounds, seed[1], seed[2]) then
        return points, "outside"
    end
    local h, length = o.step, 0
    for _ = 1, o.max_steps do
        local a = points[#points]
        local n, why = normfield(m, o, a[1], a[2])
        if not n then
            return points, why
        end
        local clearance = E.clearance(m, a[1], a[2])
        if clearance < o.tolerance then
            return points, "source"
        end
        h = min(h, o.max_length - length, max(o.tolerance, clearance * 0.45))
        if h < o.tolerance * 0.01 then
            return points, "length"
        end
        local z, error, accepted
        for _ = 1, 40 do
            local full = rk(m, o, a, sign * h)
            local mid = rk(m, o, a, sign * h / 2)
            z = mid and rk(m, o, mid, sign * h / 2)
            if full and z then
                error = sqrt((full[1] - z[1]) ^ 2 + (full[2] - z[2]) ^ 2) / 15
                if error <= o.tolerance then
                    accepted = true
                    break
                end
            end
            h = h / 2
            if h < o.tolerance * 0.001 then
                return points, "critical point"
            end
        end
        check(accepted, "streamline tolerance could not be attained")
        local stop
        z, stop = clip(m, o.bounds, a, z)
        local ds = sqrt((z[1] - a[1]) ^ 2 + (z[2] - a[2]) ^ 2)
        if ds < o.tolerance * 0.001 then
            return points, stop or "critical point"
        end
        points[#points + 1] = z
        length = length + ds
        if stop then
            return points, stop
        end
        if length >= o.max_length - o.tolerance * 0.01 then
            return points, "length"
        end
        h = min(o.step, h * min(2, 0.9 * (o.tolerance / max(error, 1e-30)) ^ 0.2))
    end
    return points, "steps"
end
--- Return polylines ordered in the physical E direction, with stop diagnostics.
function R.lines(m, options)
    if m.hybrid then
        return require("luafemm-electric-hybrid-view").lines(m, R.options(m, options))
    end
    if m.fe then
        return require("luafemm-electric-fe-view").lines(m, R.options(m, options))
    end
    E.freeze(m)
    local o = R.options(m, options)
    if m.classic then
        return require("luafemm-electric-classics").lines(m, o)
    end
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
                or (sign == "positive" and s.q > 0)
                or (sign == "negative" and s.q < 0)
            then
                for i = 1, o.lines do
                    local a = 2 * math.pi * (i - 0.5) / o.lines
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
            back, bs = halftrace(m, o, seed, -1)
        end
        if o.direction ~= "backward" then
            forward, fs = halftrace(m, o, seed, 1)
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
--- Conservative masking: discard cells touching any source disk.
local function masked(m, x0, y0, x1, y1)
    for _, s in ipairs(m.sources) do
        local x, y = max(x0, min(x1, s.x)), max(y0, min(y1, s.y))
        if (s.x - x) ^ 2 + (s.y - y) ^ 2 <= s.radius * s.radius then
            return true
        end
    end
    return false
end
function R.contours(m, options, sampler)
    if not sampler and m.hybrid then
        return require("luafemm-electric-hybrid-view").contours(m, R.options(m, options))
    end
    if not sampler and m.fe then
        return require("luafemm-electric-fe-view").contours(m, R.options(m, options))
    end
    E.freeze(m)
    local o = R.options(m, options)
    local b = o.bounds
    local nx = o.grid
    local ny = max(4, math.ceil(nx * (b[4] - b[2]) / (b[3] - b[1])))
    check((nx + 1) * (ny + 1) <= 1000000, "contour sampling budget exceeded")
    local grid, low, high = {}, math.huge, -math.huge
    local function id(i, j)
        return j * (nx + 1) + i + 1
    end
    for j = 0, ny do
        for i = 0, nx do
            local x, y = b[1] + (b[3] - b[1]) * i / nx, b[2] + (b[4] - b[2]) * j / ny
            local f = sampler and sampler.sample(x, y) or (not sampler and P.sample(m, x, y))
            local v = f and f.V
            grid[id(i, j)] = { x, y, v }
            if v and (sampler and sampler.clearance(x, y) or E.clearance(m, x, y)) > 0 then
                low = min(low, v)
                high = max(high, v)
            end
        end
    end
    local levels = o.levels or {}
    if not o.levels and low < high then
        -- Symmetric levels give an immediately recognizable signed dipole.
        if low < 0 and high > 0 then
            local a = min(-low, high)
            low, high = -a, a
        end
        for i = 1, o.count do
            levels[i] = low + (high - low) * i / (o.count + 1)
        end
    end
    local result = {}
    for _, level in ipairs(levels) do
        local edges, vertices, adj = {}, {}, {}
        local function crossing(a, c)
            local key = min(a, c) .. ":" .. max(a, c)
            if not vertices[key] then
                local p, q = grid[a], grid[c]
                local t = (level - p[3]) / (q[3] - p[3])
                vertices[key] = { p[1] + t * (q[1] - p[1]), p[2] + t * (q[2] - p[2]) }
                adj[key] = {}
            end
            return key
        end
        local function triangle(a, c, d)
            local ids = { a, c, d }
            local hits = {}
            for k = 1, 3 do
                local u, v = ids[k], ids[k % 3 + 1]
                if (grid[u][3] >= level) ~= (grid[v][3] >= level) then
                    hits[#hits + 1] = crossing(u, v)
                end
            end
            if #hits == 2 then
                local n = #edges + 1
                edges[n] = { hits[1], hits[2] }
                for _, key in ipairs(hits) do
                    adj[key][#adj[key] + 1] = n
                end
            end
        end
        for j = 0, ny - 1 do
            for i = 0, nx - 1 do
                local a, c, d, e = id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)
                if
                    grid[a][3]
                    and grid[c][3]
                    and grid[d][3]
                    and grid[e][3]
                    and not (
                        sampler and sampler.mask(grid[a][1], grid[a][2], grid[d][1], grid[d][2])
                        or (
                            not sampler
                            and masked(m, grid[a][1], grid[a][2], grid[d][1], grid[d][2])
                        )
                    )
                then
                    triangle(a, c, d)
                    triangle(a, d, e)
                end
            end
        end
        local used = {}
        local function walk(start)
            local points, key = {}, start
            while true do
                points[#points + 1] = vertices[key]
                local nextedge
                for _, n in ipairs(adj[key]) do
                    if not used[n] then
                        nextedge = n
                        break
                    end
                end
                if not nextedge then
                    break
                end
                used[nextedge] = true
                local edge = edges[nextedge]
                key = edge[1] == key and edge[2] or edge[1]
                if key == start then
                    points[#points + 1] = vertices[key]
                    break
                end
            end
            if #points > 1 then
                result[#result + 1] = { points = points, level = level }
            end
        end
        -- Sorted edge order preserves deterministic output (pairs does not).
        for _, edge in ipairs(edges) do
            for _, key in ipairs(edge) do
                if #adj[key] == 1 and not used[adj[key][1]] then
                    walk(key)
                end
            end
        end
        for n, edge in ipairs(edges) do
            if not used[n] then
                walk(edge[1])
            end
        end
    end
    return result
end
function R.vectors(m, options)
    E.freeze(m)
    local o = R.options(m, options)
    local b = o.bounds
    local result = {}
    local nx, ny = math.floor((b[3] - b[1]) / o.spacing), math.floor((b[4] - b[2]) / o.spacing)
    check((nx + 1) * (ny + 1) <= 100000, "vector sampling budget exceeded")
    for j = 0, ny do
        for i = 0, nx do
            local x, y = b[1] + (i + 0.5) * o.spacing, b[2] + (j + 0.5) * o.spacing
            if inside(b, x, y) and E.clearance(m, x, y) > 0 then
                local f = P.sample(m, x, y)
                if f and f.E > o.null_field then
                    local l = o.vector_mode == "normalized" and o.vector_length
                        or f.E * o.vector_scale
                    local z, why = clip(m, b, { x, y }, { x + l * f.Ex / f.E, y + l * f.Ey / f.E })
                    -- Suppress cropped arrows: their length would misrepresent E.
                    if
                        not why
                        and not (m.fe and require("luafemm-electric-fe-view").crosses_interface(
                            m,
                            { x, y },
                            z
                        ))
                        and not (
                            m.classic
                            and require("luafemm-electric-classics").crosses_jump(m, { x, y }, z)
                        )
                    then
                        result[#result + 1] = { points = { { x, y }, z }, magnitude = f.E }
                    end
                end
            end
        end
    end
    return result
end
--- Uniform arclength samples of one flattened path, with owned source data.
-- Invalid samples become explicit gaps only under the requested nan policy.
function R.profile(m, points, options)
    if m.fe then
        return require("luafemm-electric-fe-view").profile(m, points, options)
    end
    local o = options or {}
    local n = E.integer(o.samples or 201, "profile samples", 2, 100000)
    local outside = o.outside or "error"
    check(outside == "error" or outside == "nan", "unknown outside policy")
    check(#points >= 2, "profile requires a nonempty contour")
    local lengths = { 0 }
    local total = 0
    for i = 2, #points do
        total = total
            + sqrt((points[i][1] - points[i - 1][1]) ^ 2 + (points[i][2] - points[i - 1][2]) ^ 2)
        lengths[i] = total
    end
    check(total > 0, "profile requires positive length")
    E.freeze(m)
    local targets = {}
    for i = 0, n - 1 do
        targets[#targets + 1] = total * i / (n - 1)
    end
    -- Insert a gap even when a narrow excluded disk falls between two samples.
    -- Intersections use the actual flattened path, not chords across bends.
    for i = 2, #points do
        local a, b = points[i - 1], points[i]
        local dx, dy = b[1] - a[1], b[2] - a[2]
        local aa = dx * dx + dy * dy
        if aa > 0 then
            for _, source in ipairs(m.sources) do
                local x, y = a[1] - source.x, a[2] - source.y
                local bb = x * dx + y * dy
                local cc = x * x + y * y - source.radius ^ 2
                local disc = bb * bb - aa * cc
                if disc >= 0 then
                    local lo, hi = max(0, (-bb - sqrt(disc)) / aa), min(1, (-bb + sqrt(disc)) / aa)
                    if lo <= hi then
                        check(outside == "nan", "invalid profile sample: source exclusion")
                        targets[#targets + 1] = lengths[i - 1] + (lo + hi) * sqrt(aa) / 2
                    end
                end
            end
        end
    end
    local events = m.classic
            and require("luafemm-electric-classics").profile_events(m, points, lengths)
        or {}
    local scheduled = {}
    for _, s in ipairs(targets) do
        local at_interface = false
        for _, event in ipairs(events) do
            if abs(s - event.s) <= 1e-12 * total then
                at_interface = true
                break
            end
        end
        if not at_interface then
            scheduled[#scheduled + 1] = { s = s, order = 0 }
        end
    end
    for _, event in ipairs(events) do
        scheduled[#scheduled + 1] = { s = event.s, side = event.left, order = 1 }
        if event.jump then
            scheduled[#scheduled + 1] = { s = event.s, side = event.left, gap = true, order = 2 }
            scheduled[#scheduled + 1] = { s = event.s, side = event.right, order = 3 }
        end
    end
    table.sort(scheduled, function(a, b)
        return a.s == b.s and a.order < b.order or a.s < b.s
    end)
    local data = {}
    local k = 2
    for _, target in ipairs(scheduled) do
        local s = target.s
        while k < #points and lengths[k] < s do
            k = k + 1
        end
        local l = lengths[k] - lengths[k - 1]
        local t = l > 0 and (s - lengths[k - 1]) / l or 0
        local a, b = points[k - 1], points[k]
        local x, y = a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])
        local f, why = P.sample(m, x, y, target.side)
        if E.clearance(m, x, y) <= 0 then
            f, why = nil, "source exclusion"
        end
        check(f or outside == "nan", "invalid profile sample: " .. tostring(why))
        f = f or { V = 0 / 0, Ex = 0 / 0, Ey = 0 / 0, E = 0 / 0, Dx = 0 / 0, Dy = 0 / 0, D = 0 / 0 }
        if target.gap then
            for _, key in ipairs({ "Ex", "Ey", "E", "Dx", "Dy", "D" }) do
                f[key] = 0 / 0
            end
        end
        f.x, f.y, f.s, f.t = x, y, s, s / total
        data[#data + 1] = f
    end
    return { model = E.snapshot(m), points = P.copy(points), data = data, length = total }
end
return R
