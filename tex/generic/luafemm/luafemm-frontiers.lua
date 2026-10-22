-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Native computational contours and tagged, rigidly paired boundaries.
-- Geometry is passive model-coordinate data. Pairing is strong nodal equality,
-- never a penalty approximation; paired edges are synchronised before assembly.
-- @module luafemm-frontiers
local F = {}
local top = require("luafemm-topology")
local function check(v, s)
    assert(v, "luafemm frontiers: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
local function copy(t)
    if type(t) ~= "table" then
        return t
    end
    local r = {}
    for k, v in pairs(t) do
        r[k] = copy(v)
    end
    return r
end
local function editable(m)
    check(m and not m.nodes, "declare frontiers before meshing")
    check(not m.import_document, "imported topology is immutable")
    check(m.options.field_model == "open", "frontiers require field model=open")
    check(m.options.mesher == "delaunay", "path frontiers require mesher=delaunay")
    check(not m.options.boundary, "path frontiers cannot use a boundary callback")
end
local function points(m, p)
    check(type(p) == "table" and #p >= 2, "path needs at least two points")
    for _, v in ipairs(p) do
        check(finite(v[1]) and finite(v[2]), "invalid path coordinate")
        check(
            v[1] >= m.options.xmin
                and v[1] <= m.options.xmax
                and v[2] >= m.options.ymin
                and v[2] <= m.options.ymax,
            "path outside problem bounds"
        )
    end
end

--- Set one compound computational domain, including even-odd holes.
-- The problem rectangle remains a staging/meshing box, not an extra domain.
-- @tparam table m Native unmeshed open-field model.
-- @tparam table contours Closed polygon loops, without repeated endpoints.
-- @tparam[opt="even odd"] string fill_rule PGF-compatible fill rule.
function F.domain(m, contours, fill_rule)
    editable(m)
    check(not m.domain, "only one computational domain is allowed")
    check(#contours > 0, "empty domain")
    for _, p in ipairs(contours) do
        points(m, p)
        require("luafemm-mesh").validate(p)
    end
    fill_rule = fill_rule or "even odd"
    check(fill_rule == "even odd" or fill_rule == "nonzero", "unknown domain fill rule")
    m.domain = { contours = copy(contours), fill_rule = fill_rule }
end

--- Attach a named physical condition to an open or closed polygonal path.
-- Internal paths must bound faces; only Dirichlet internal conditions are valid.
-- open: first exterior multipole, H.t = order/(mu0*radius)*A, radius in model units.
-- @tparam table m Native unmeshed model.
-- @tparam string name Unique physical path name.
-- @tparam table p Model-coordinate vertices, without a repeated closure.
-- @tparam boolean closed Whether the last vertex connects to the first.
-- @tparam[opt] table o Boundary type/data; values SI, geometric lengths model units.
function F.path(m, name, p, closed, o)
    editable(m)
    points(m, p)
    check(type(name) == "string" and name:match("%S"), "frontier needs a name")
    m.frontiers = m.frontiers or {}
    check(not m.frontiers[name], "duplicate frontier " .. name)
    o = copy(o or {})
    local kind = o.type or "dirichlet"
    check(
        kind == "dirichlet"
            or kind == "neumann"
            or kind == "robin"
            or kind == "periodic"
            or kind == "open",
        "unknown frontier type"
    )
    local c = {
        type = kind,
        value = o.value or 0,
        dx = o.dx or 0,
        dy = o.dy or 0,
        coefficient = o.coefficient or 0,
    }
    for _, k in ipairs({ "value", "dx", "dy", "coefficient" }) do
        check(finite(c[k]), "invalid " .. k)
    end
    check(c.coefficient >= 0, "negative Robin coefficient")
    check(kind == "robin" or c.coefficient == 0, "coefficient requires robin")
    if kind == "open" then
        check(
            c.value == 0 and c.dx == 0 and c.dy == 0,
            "open boundaries do not accept imposed loads"
        )
        local radius, order = o.radius, o.order or 1
        check(
            finite(radius) and radius > 0 and finite(order) and order > 0,
            "open boundary needs positive radius and order"
        )
        c.type, c.coefficient = "robin", order / (4 * math.pi * 1e-7 * radius * m.options.unit)
        c.value, c.dx, c.dy = 0, 0, 0
        local cx, cy = o.center_x or 0, o.center_y or 0
        check(finite(cx) and finite(cy), "invalid open boundary centre")
        for _, v in ipairs(p) do
            local r = math.sqrt((v[1] - cx) ^ 2 + (v[2] - cy) ^ 2)
            check(
                math.abs(r - radius) <= math.max(0.001, o.tolerance or 0, radius * 3e-4),
                "open boundary vertices must lie on the specified circle"
            )
        end
        check(closed, "open boundary approximation requires a complete closed circle")
        c.open_radius = radius
    end
    if kind == "periodic" then
        check(c.value == 0 and c.dx == 0 and c.dy == 0, "periodic paths cannot prescribe loads")
        c.path = name
    end
    m.frontiers[name] = { points = copy(p), closed = closed or false, condition = c }
end

--- Pair named periodic paths using a rotation about (cx,cy), then translation.
-- The transform maps first to second; sign is +1 or -1.
-- @tparam table m Native unmeshed model.
-- @tparam string first Source frontier identifier.
-- @tparam string second Target frontier identifier.
-- @tparam[opt] table o angle (degrees), x/y/cx/cy/tolerance (model units), sign.
function F.pair(m, first, second, o)
    editable(m)
    o = o or {}
    check(first ~= second, "pair needs two distinct paths")
    local p = {
        first = first,
        second = second,
        angle = o.angle or 0,
        x = o.x or 0,
        y = o.y or 0,
        cx = o.cx or 0,
        cy = o.cy or 0,
        sign = o.sign or 1,
        tolerance = o.tolerance or 0.001,
    }
    for _, k in ipairs({ "angle", "x", "y", "cx", "cy", "sign", "tolerance" }) do
        check(finite(p[k]), "invalid pair " .. k)
    end
    check(p.tolerance > 0, "pair tolerance must be positive")
    check(p.sign == 1 or p.sign == -1, "pair sign must be 1 or -1")
    m.pairs = m.pairs or {}
    m.pairs[#m.pairs + 1] = p
end

--- Compile once per physical preparation, retaining input regions unmodified.
-- @tparam table m Model whose physical declarations are complete.
function F.compile(m)
    if m.import_document then
        return
    end
    if not (m.domain or m.frontiers or m.circuits) then
        return
    end
    check(
        m.options.field_model == "open" and m.options.mesher == "delaunay",
        "static geometry requires open/delaunay"
    )
    local used = {}
    for _, p in ipairs(m.pairs or {}) do
        for _, name in ipairs({ p.first, p.second }) do
            check(
                m.frontiers and m.frontiers[name] and m.frontiers[name].condition.type == "periodic",
                "pair requires a periodic path: " .. tostring(name)
            )
            check(not used[name], "a frontier cannot belong to two pairs")
            used[name] = true
        end
    end
    for name, f in pairs(m.frontiers or {}) do
        check(f.condition.type ~= "periodic" or used[name], "unpaired frontier " .. name)
    end
    m.topology = top.native(m)
    m.native_topology = true
    require("luafemm-static").excite(m, m.topology)
end
local function transform(p, v, inverse)
    local c, s = math.cos(math.rad(p.angle)), math.sin(math.rad(p.angle))
    if inverse then
        local x, y = v[1] - p.cx - p.x, v[2] - p.cy - p.y
        return { p.cx + c * x + s * y, p.cy - s * x + c * y }
    end
    local x, y = v[1] - p.cx, v[2] - p.cy
    return { p.cx + c * x - s * y + p.x, p.cy + s * x + c * y + p.y }
end
local function on(m, name, v, epsilon)
    for _, s in ipairs(m.topology.segments) do
        if
            s.condition
            and s.condition.path == name
            and top.distance(v, m.topology.points[s[1]], m.topology.points[s[2]])
                <= (epsilon or m.topology.tolerance * 30)
        then
            return true
        end
    end
    return false
end
local function boundary_edges(elements)
    local edges = {}
    for ti, e in ipairs(elements) do
        for j = 1, 3 do
            local a, b = e[j], e[j % 3 + 1]
            local key = math.min(a, b) .. ":" .. math.max(a, b)
            local r = edges[key]
                or { a, b, triangle = ti, opposite = e[(j + 1) % 3 + 1], count = 0 }
            r.count = r.count + 1
            edges[key] = r
        end
    end
    return edges
end

--- Insert the union of transformed boundary vertices on both paired traces.
-- Splitting a boundary triangle preserves conformity and positive area. The
-- carved mesh audit recomputes area, quality and Euler checks afterwards.
-- @tparam table m Model with its compiled arrangement.
-- @tparam table nodes Mutable vertices in model units.
-- @tparam table elements Mutable counterclockwise index triples.
-- @tparam table stats Mesh diagnostics and protected segment indices.
-- @return Vertices, elements and recomputed carved-mesh diagnostics.
function F.synchronize(m, nodes, elements, stats)
    if not m.pairs or #m.pairs == 0 then
        return nodes, elements, stats
    end
    local epsilon = m.topology.tolerance * 30
    local targets = {}
    local edges = boundary_edges(elements)
    local exposed = {}
    for _, e in pairs(edges) do
        if e.count == 1 then
            exposed[e[1]] = true
            exposed[e[2]] = true
        end
    end
    for _, p in ipairs(m.pairs) do
        local counts = { 0, 0 }
        for id, v in ipairs(nodes) do
            for k, name in ipairs({ p.first, p.second }) do
                if on(m, name, v) then
                    check(exposed[id], "periodic conditions require exterior paths")
                    counts[k] = counts[k] + 1
                    local q = transform(p, v, k == 2)
                    local target = k == 1 and p.second or p.first
                    check(on(m, target, q, p.tolerance), "paired paths are not rigid images")
                    -- Project within the declared coordinate tolerance onto
                    -- the actual target trace (PGF coordinates are rounded).
                    local best, projected = math.huge
                    for _, e in pairs(edges) do
                        if
                            e.count == 1
                            and on(m, target, nodes[e[1]])
                            and on(m, target, nodes[e[2]])
                        then
                            local a, b = nodes[e[1]], nodes[e[2]]
                            local d, t = top.distance(q, a, b)
                            if d < best then
                                t = math.max(0, math.min(1, t))
                                best = d
                                projected = { a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]) }
                                for _, v in ipairs({ a, b }) do
                                    if (v[1] - q[1]) ^ 2 + (v[2] - q[2]) ^ 2 < p.tolerance ^ 2 then
                                        projected = { v[1], v[2] }
                                        break
                                    end
                                end
                            end
                        end
                    end
                    check(projected and best <= p.tolerance, "pair target is not exterior")
                    targets[#targets + 1] = projected
                end
            end
        end
        check(counts[1] >= 2 and counts[2] >= 2, "periodic path absent from mesh")
    end
    local added = 0
    for _, q in ipairs(targets) do
        local found = false
        for _, v in ipairs(nodes) do
            if (v[1] - q[1]) ^ 2 + (v[2] - q[2]) ^ 2 <= epsilon ^ 2 then
                found = true
                break
            end
        end
        if not found then
            edges = boundary_edges(elements)
            local edge
            for _, e in pairs(edges) do
                if e.count == 1 and top.distance(q, nodes[e[1]], nodes[e[2]]) <= epsilon then
                    edge = e
                    break
                end
            end
            check(edge, "cannot synchronise paired boundary")
            local id = #nodes + 1
            check(id <= 40000, "paired mesh exceeds the 40000-vertex limit")
            nodes[id] = q
            local a, b, c = edge[1], edge[2], edge.opposite
            elements[edge.triangle] = { a, id, c }
            elements[#elements + 1] = { id, b, c }
            for si, s in ipairs(stats.segments) do
                if (s[1] == a and s[2] == b) or (s[1] == b and s[2] == a) then
                    stats.segments[si] = { a, id }
                    stats.segments[#stats.segments + 1] = { id, b }
                    break
                end
            end
            added = added + 1
        end
    end
    nodes, elements, stats = top.carve(m.topology, nodes, elements, stats)
    stats.paired_vertices = added
    return nodes, elements, stats
end

-- Signed union-find: value(i)=weight(i)*value(root). Odd signed cycles fix zero.
local function union_find(n)
    local parent, weight, zero = {}, {}, {}
    for i = 1, n do
        parent[i], weight[i] = i, 1
    end
    local function root(i)
        if parent[i] ~= i then
            local r, s = root(parent[i])
            weight[i] = weight[i] * s
            parent[i] = r
        end
        return parent[i], weight[i]
    end
    local function join(a, b, sign)
        local ra, sa = root(a)
        local rb, sb = root(b)
        if ra == rb then
            if sa ~= sign * sb then
                zero[ra] = true
            end
            return
        end
        if ra > rb then
            ra, rb = rb, ra
            sa, sb = sb, sa
        end
        parent[rb], weight[rb] = ra, sign * sa * sb
        zero[ra] = zero[ra] or zero[rb]
    end
    return root, join, zero
end

--- Eliminate signed periodic constraints and choose only genuine null gauges.
-- All residual/tangent contributions must subsequently use nodal dof_sign.
-- @tparam table m Model with rigid pairing declarations.
-- @tparam table nodes SI vertices; signed constraint metadata are attached.
-- @tparam table elements Counterclockwise vertex index triples.
-- @tparam table fixed Prescribed potentials; propagated in place.
-- @tparam[opt] table gauges Preliminary additive gauges, removed before pairing.
-- @tparam[opt] table components Preliminary geometric components.
-- @tparam[opt] table edges Weak boundary edges, including Robin anchors.
-- @return Representatives, surviving null modes and signed nodal pairs.
function F.reduce(m, nodes, elements, fixed, gauges, components, edges)
    if not m.pairs or #m.pairs == 0 then
        return nil, nil, nil
    end
    for _, id in pairs(gauges or {}) do
        fixed[id] = nil
    end
    local root, join, zero = union_find(#nodes)
    local nr, nj, nz = union_find(#nodes)
    for _, e in ipairs(elements) do
        nj(e[1], e[2], 1)
        nj(e[1], e[3], 1)
    end
    local epsilon = m.topology.tolerance * 30
    local pairs_list = {}
    for _, p in ipairs(m.pairs) do
        local matched = {}
        for i, v in ipairs(nodes) do
            local a = { v[1] / m.options.unit, v[2] / m.options.unit }
            if on(m, p.first, a) then
                local q = transform(p, a, false)
                local match, best = nil, p.tolerance ^ 2
                for j, w in ipairs(nodes) do
                    local distance = (w[1] / m.options.unit - q[1]) ^ 2
                        + (w[2] / m.options.unit - q[2]) ^ 2
                    if
                        distance <= best
                        and on(m, p.second, { w[1] / m.options.unit, w[2] / m.options.unit })
                    then
                        match, best = j, distance
                    end
                end
                check(match, "paired mesh has an unmatched node")
                check(
                    not matched[match],
                    "pair tolerance merged distinct vertices; reduce tolerance"
                )
                matched[match] = true
                join(i, match, p.sign)
                nj(i, match, p.sign)
                pairs_list[#pairs_list + 1] = { i, match, p.sign }
            end
        end
        for j, w in ipairs(nodes) do
            if on(m, p.second, { w[1] / m.options.unit, w[2] / m.options.unit }) then
                check(matched[j], "paired target has an unmatched node")
            end
        end
    end
    local prescribed, anchored = {}, {}
    for id in ipairs(nodes) do
        local value = fixed[id]
        if value ~= nil then
            local r, s = root(id)
            local v = value * s
            check(
                prescribed[r] == nil or math.abs(prescribed[r] - v) < 1e-12,
                "conflicting periodic Dirichlet values"
            )
            prescribed[r] = v
            anchored[nr(id)] = true
        end
    end
    for _, e in ipairs(edges or {}) do
        if e.coefficient > 0 then
            anchored[nr(e.ids[1])] = true
        end
    end
    local nulls = {}
    for i in ipairs(nodes) do
        local r, s = root(i)
        local nc, ns = nr(i)
        if zero[r] then
            check(
                not prescribed[r] or math.abs(prescribed[r]) < 1e-12,
                "antiperiodic corner must have zero potential"
            )
            prescribed[r] = 0
        end
        if not anchored[nc] and not nz[nc] and not nulls[nc] then
            nulls[nc] = { node = i }
            prescribed[r] = 0
        end
        nodes[i].null_component = nc
        nodes[i].null_sign = ns
    end
    local reps = {}
    for i, node in ipairs(nodes) do
        local r, s = root(i)
        reps[i] = r
        node.dof_sign = s
        fixed[i] = prescribed[r] and s * prescribed[r] or nil
    end
    m.periodic_nodes = pairs_list
    m.periodic_nulls = nulls
    return reps, nulls, pairs_list
end

--- Verify the air/multipole assumptions on native circular open boundaries.
-- @tparam table triangles Prepared elements, with material/source data.
-- @tparam table edges Prepared weak edges with optional open-radius tags.
function F.open_check(triangles, edges)
    local open = {}
    for _, edge in ipairs(edges) do
        if edge.open_radius then
            open[math.min(edge.ids[1], edge.ids[2]) .. ":" .. math.max(edge.ids[1], edge.ids[2])] =
                true
        end
    end
    if not next(open) then
        return
    end
    local current, scale = 0, 0
    for _, e in ipairs(triangles) do
        current = current + e.current * e.area
        scale = scale + math.abs(e.current * e.area)
        for j = 1, 3 do
            local a, b = e.ids[j], e.ids[j % 3 + 1]
            if open[math.min(a, b) .. ":" .. math.max(a, b)] then
                check(
                    not e.p.bh
                        and e.p.mur_x == 1
                        and e.p.mur_y == 1
                        and e.p.hc == 0
                        and e.current == 0,
                    "open multipole boundary must be surrounded by source-free air"
                )
            end
        end
    end
    check(
        math.abs(current) <= 1e-10 * math.max(scale, 1e-30),
        "open multipole boundary requires zero net longitudinal current"
    )
end

--- Check the Neumann compatibility against each surviving signed null mode.
-- @tparam table m Model carrying its surviving signed null modes.
-- @tparam table triangles Prepared elements.
-- @tparam table edges Prepared weak boundary edges.
-- @tparam table nodes Vertices carrying component/null-sign metadata.
function F.compatibility(m, triangles, edges, nodes)
    for nc in pairs(m.periodic_nulls or {}) do
        local sum, scale = 0, 0
        for _, e in ipairs(triangles) do
            if nodes[e.ids[1]].null_component == nc then
                local v = nodes[e.ids[1]].null_sign * e.current * e.area
                sum, scale = sum + v, scale + math.abs(v)
            end
        end
        for _, e in ipairs(edges) do
            if nodes[e.ids[1]].null_component == nc then
                local v = nodes[e.ids[1]].null_sign * e.length * (e.g[1] + e.g[2]) / 2
                sum, scale = sum - v, scale + math.abs(v)
            end
        end
        check(
            math.abs(sum) <= 1e-10 * math.max(scale, 1e-30),
            "incompatible periodic Neumann current balance"
        )
    end
end
return F
