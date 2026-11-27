-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Planar residual indicators and conforming local red/green refinement.
-- Element current residuals and tangential H jumps estimate field error.
-- The estimator is not a certified bound, especially at material corners.
-- Marking is deterministic bulk marking; every split edge is shared by both
-- neighbours. No material boundary, hole or narrow feature is moved or erased.
local A = {}
local M = require("luafemm")
local copy = require("luafemm-time-state").copy
local function check(v, s)
    assert(v, "luafemm adaptive: " .. s)
end
local function edge(a, b)
    if a > b then
        a, b = b, a
    end
    return a .. ":" .. b
end
local function eligible(m)
    check(
        m.options.cache.mode == "off",
        "explicit adaptation requires cache mode off; save .ansx answers"
    )
    check(
        m.options.problem_type == "planar"
            and m.options.field_model == "open"
            and not m.import_document
            and not m.import_answer
            and not m.operating
            and not m.sliding
            and not m.transient
            and not m.network
            and not require("luafemm-flux-drive").active(m)
            and (not m.pairs or #m.pairs == 0),
        "requires an ordinary planar model without periodicity, motion, imported state or controllers"
    )
end
function A.estimate(m)
    eligible(m)
    M.solve(m)
    local indicators, edges, fields = {}, {}, {}
    local total, energy = 0, 0
    for i, e in ipairs(m.triangles) do
        local x, y = 0, 0
        for _, id in ipairs(e.ids) do
            x = x + m.nodes[id][1] / 3
            y = y + m.nodes[id][2] / 3
        end
        local f = m.options.frequency > 0 and require("luafemm-harmonic").field(m, e, x, y)
            or { Hx = { e.hx, 0 }, Hy = { e.hy, 0 }, J = { e.current, 0 } }
        fields[i] = f
        local h2 = 0
        for j = 1, 3 do
            local u, v = e.ids[j], e.ids[j % 3 + 1]
            local p, q = m.nodes[u], m.nodes[v]
            h2 = math.max(h2, (p[1] - q[1]) ^ 2 + (p[2] - q[2]) ^ 2)
            local key = edge(u, v)
            local list = edges[key] or { u = u, v = v }
            edges[key] = list
            list[#list + 1] = i
        end
        local nu = (e.p.tensor[1] + e.p.tensor[3]) / 2
        if e.p.bh then
            nu = M.constitutive(e.p, e.b)
        end
        indicators[i] = h2 * e.area * M.mu0 * (f.J[1] ^ 2 + f.J[2] ^ 2) / nu
        local b2 = e.b ^ 2
        if m.options.frequency > 0 then
            b2 = f.Bx[1] ^ 2 + f.Bx[2] ^ 2 + f.By[1] ^ 2 + f.By[2] ^ 2
        end
        energy = energy + e.area * nu * b2 / M.mu0
    end
    for _, ids in pairs(edges) do
        if #ids == 2 then
            local p, q = m.nodes[ids.u], m.nodes[ids.v]
            local dx, dy = q[1] - p[1], q[2] - p[2]
            local l = math.sqrt(dx * dx + dy * dy)
            local a, b = fields[ids[1]], fields[ids[2]]
            local jump = 0
            for part = 1, 2 do
                jump = jump
                    + (((a.Hx[part] - b.Hx[part]) * dx + (a.Hy[part] - b.Hy[part]) * dy) / l)
                        ^ 2
            end
            for _, id in ipairs(ids) do
                local p = m.triangles[id].p
                local nu = (p.tensor[1] + p.tensor[3]) / 2
                indicators[id] = indicators[id] + 0.5 * l * l * M.mu0 * jump / nu
            end
        end
    end
    for _, v in ipairs(indicators) do
        total = total + v
    end
    return indicators, math.sqrt(total / math.max(energy, 1e-300)), total
end
--- Split marked cells conformingly, including unmarked neighbours sharing edges.
function A.split(m, marked, maxnodes)
    local nodes, split, counts = {}, {}, {}
    for i, p in ipairs(m.nodes) do
        nodes[i] = { p[1], p[2] }
    end
    for i, e in ipairs(m.triangles) do
        for j = 1, 3 do
            local a, b = e.ids[j], e.ids[j % 3 + 1]
            local k = edge(a, b)
            counts[k] = (counts[k] or 0) + 1
            if marked[i] and not split[k] then
                local p, q = nodes[a], nodes[b]
                nodes[#nodes + 1] = { (p[1] + q[1]) / 2, (p[2] + q[2]) / 2 }
                split[k] = #nodes
            end
        end
    end
    check(#nodes <= maxnodes, "node budget exceeded before refinement")
    local elements = {}
    local function tri(a, b, c)
        elements[#elements + 1] = { a, b, c }
    end
    for _, e in ipairs(m.triangles) do
        local a, b, c = table.unpack(e.ids)
        local u, v, w = split[edge(a, b)], split[edge(b, c)], split[edge(c, a)]
        local n = (u and 1 or 0) + (v and 1 or 0) + (w and 1 or 0)
        if n == 0 then
            tri(a, b, c)
        elseif n == 3 then
            tri(a, u, w)
            tri(u, b, v)
            tri(w, v, c)
            tri(u, v, w)
        else
            if n == 1 then
                while not u do
                    a, b, c = b, c, a
                    u, v, w = v, w, u
                end
                tri(a, u, c)
                tri(u, b, c)
            else
                while w do
                    a, b, c = b, c, a
                    u, v, w = v, w, u
                end
                tri(b, v, u)
                tri(a, u, c)
                tri(u, v, c)
            end
        end
    end
    local boundary = copy(m.exterior_nodes)
    for k, id in pairs(split) do
        if counts[k] == 1 then
            boundary[id] = true
        end
    end
    return nodes, elements, boundary, split
end
--- Refine transactionally: failure leaves the previous solved model available.
-- Options: passes integer 1..12, fraction in (0,1], tolerance>=0, max_nodes>=3.
function A.run(m, o)
    eligible(m)
    o = o or {}
    local passes, fraction, tolerance, maxnodes =
        tonumber(o.passes or 3),
        tonumber(o.fraction or 0.5),
        tonumber(o.tolerance or 0.05),
        tonumber(o.max_nodes or 100000)
    check(passes and passes % 1 == 0 and passes >= 1 and passes <= 12, "passes must be 1..12")
    check(fraction and fraction > 0 and fraction <= 1, "fraction must be in (0,1]")
    check(tolerance and tolerance >= 0 and tolerance < math.huge, "invalid tolerance")
    check(
        maxnodes and maxnodes % 1 == 0 and maxnodes >= 3 and maxnodes <= 1000000,
        "invalid node budget"
    )
    local start, mem = os.clock(), collectgarbage("count")
    local history = {}
    for pass = 0, passes do
        local indicator, relative, total = A.estimate(m)
        history[#history + 1] = { nodes = #m.nodes, triangles = #m.triangles, indicator = relative }
        m.stats.adaptive_error, m.stats.adaptive_passes = relative, pass
        if relative <= tolerance or pass == passes then
            break
        end
        local order = {}
        for i in ipairs(indicator) do
            order[i] = i
        end
        table.sort(order, function(a, b)
            return indicator[a] == indicator[b] and a < b or indicator[a] > indicator[b]
        end)
        local marked, sum = {}, 0
        for _, i in ipairs(order) do
            marked[i] = true
            sum = sum + indicator[i]
            if sum >= fraction * total then
                break
            end
        end
        local nodes, elements, boundary, split = A.split(m, marked, maxnodes)
        local options = copy(m.options)
        options.cache = { mode = "off" }
        local nextmodel = M.new(options)
        for _, key in ipairs({
            "materials",
            "regions",
            "boundaries",
            "frontiers",
            "domain",
            "topology",
            "native_topology",
            "circuits",
            "point_constraints",
            "frame",
            "component_basis",
        }) do
            nextmodel[key] = copy(m[key])
        end
        nextmodel.mesh_stats = copy(m.mesh_stats or {})
        local stats = nextmodel.mesh_stats
        stats.boundary_nodes = copy(boundary)
        local segments = {}
        for _, s in ipairs(stats.segments or {}) do
            local middle = split[edge(s[1], s[2])]
            if middle then
                local t = copy(s)
                t[2] = middle
                segments[#segments + 1] = t
                t = copy(s)
                t[1] = middle
                segments[#segments + 1] = t
            else
                segments[#segments + 1] = copy(s)
            end
        end
        stats.segments, stats.nodes, stats.triangles = segments, #nodes, #elements
        stats.min_angle, stats.max_area = 180, 0
        for _, ids in ipairs(elements) do
            local a, b, c = nodes[ids[1]], nodes[ids[2]], nodes[ids[3]]
            local twice = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
            stats.max_area = math.max(stats.max_area, twice / (2 * options.unit ^ 2))
            for j = 1, 3 do
                local p, q, r = nodes[ids[j]], nodes[ids[j % 3 + 1]], nodes[ids[(j + 1) % 3 + 1]]
                local ux, uy, vx, vy = q[1] - p[1], q[2] - p[2], r[1] - p[1], r[2] - p[2]
                local cosine = (ux * vx + uy * vy)
                    / math.sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy))
                stats.min_angle = math.min(
                    stats.min_angle,
                    math.deg(math.acos(math.max(-1, math.min(1, cosine))))
                )
            end
        end
        check(
            stats.min_angle + 1e-8 >= options.min_angle,
            "refinement violates minimum angle; reduce the target or refine geometry"
        )
        M.prepare_mesh(nextmodel, nodes, elements, boundary)
        nextmodel.adaptive_mesh = { nodes = copy(nodes), elements = copy(elements) }
        M.solve(nextmodel)
        -- Preserve the active model's identity held by named TikZ consumers.
        for key in pairs(m) do
            m[key] = nil
        end
        for key, value in pairs(nextmodel) do
            m[key] = value
        end
    end
    m.adaptive_history = history
    m.stats.adaptive_seconds = os.clock() - start
    m.stats.adaptive_memory_kib = math.max(0, collectgarbage("count") - mem)
    return m
end
return A
