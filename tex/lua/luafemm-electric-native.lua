-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Native PSLG cells and point constraints, independent of document parsing.
-- Flood triangle adjacency across untagged edges. Material/hole seeds classify
-- connected cells; unlabelled exterior is removed, unlabelled enclosed cells fail.
-- @module luafemm-electric-native
local N = {}
local E = require("luafemm-electrostatic")
local function check(v, s)
    assert(v, "luafemm electric native: " .. s)
end
function N.on(p, a, b, tol)
    local dx, dy = b[1] - a[1], b[2] - a[2]
    local l = dx * dx + dy * dy
    if l == 0 then
        return false
    end
    local t = ((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / l
    return t >= -tol / math.sqrt(l)
        and t <= 1 + tol / math.sqrt(l)
        and math.abs(dx * (p[2] - a[2]) - dy * (p[1] - a[1])) <= tol * math.sqrt(l)
end
function N.key(a, b)
    if a > b then
        a, b = b, a
    end
    return a .. ":" .. b
end
function N.generate(m, constraints)
    local f = m.fe
    local d = f.native
    if d.solution then
        local nodes, ids, labels = {}, {}, {}
        for i, p in ipairs(d.solution.nodes) do
            nodes[i] = { p[1], p[2] }
        end
        for i, t in ipairs(d.solution.triangles) do
            ids[i] = {}
            for j = 1, 3 do
                ids[i][j] = E.integer(t[j], "solution vertex", 0, #nodes - 1) + 1
            end
            labels[i] = E.integer(t[4], "solution label", 0, #d.labels - 1) + 1
            local a, b, c = nodes[ids[i][1]], nodes[ids[i][2]], nodes[ids[i][3]]
            local area = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
            check(area ~= 0, "degenerate solution element")
            if area < 0 then
                ids[i][2], ids[i][3] = ids[i][3], ids[i][2]
            end
        end
        local used, seen = {}, {}
        for _, t in ipairs(ids) do
            local sorted = { t[1], t[2], t[3] }
            table.sort(sorted)
            local key = table.concat(sorted, ":")
            check(not seen[key], "duplicate solution element")
            seen[key] = true
            for _, id in ipairs(t) do
                used[id] = true
            end
        end
        for i = 1, #nodes do
            check(used[i], "unused solution node")
        end
        check(#ids > 0, "empty solution mesh")
        f.native_labels = labels
        f.periodic_pairs, f.periodic_nodes = {}, {}
        return nodes, ids, { mesh_tolerance = 1e-10 }
    end
    for _, s in ipairs(d.lines) do
        constraints[#constraints + 1] = { contours = { { s.a, s.b } }, mesh_size = s.h or 0 }
    end
    f.options.seed_points = {}
    for _, p in ipairs(d.points) do
        f.options.seed_points[#f.options.seed_points + 1] = { p[1], p[2] }
    end
    local nodes, ids, stats = require("luafemm-electric-periodic").generate(f, constraints)
    local edges, neighbors = {}, {}
    for i, t in ipairs(ids) do
        neighbors[i] = {}
        for j = 1, 3 do
            local a, b = t[j], t[j % 3 + 1]
            local k = N.key(a, b)
            edges[k] = edges[k] or { a, b, triangles = {} }
            table.insert(edges[k].triangles, i)
        end
    end
    local tol = math.max(
        stats.mesh_tolerance * 8,
        1e-9 * math.max(f.bounds[3] - f.bounds[1], f.bounds[4] - f.bounds[2])
    )
    local exterior = {}
    for _, e in pairs(edges) do
        local blocked = false
        for _, s in ipairs(d.lines) do
            if N.on(nodes[e[1]], s.a, s.b, tol) and N.on(nodes[e[2]], s.a, s.b, tol) then
                blocked = true
                break
            end
        end
        if not blocked then
            if #e.triangles == 2 then
                local a, b = e.triangles[1], e.triangles[2]
                table.insert(neighbors[a], b)
                table.insert(neighbors[b], a)
            else
                exterior[e.triangles[1]] = true
            end
        end
    end
    local cells, of = {}, {}
    for i = 1, #ids do
        if not of[i] then
            local ci = #cells + 1
            local c = { i }
            cells[ci] = c
            of[i] = ci
            local j = 1
            while c[j] do
                local t = c[j]
                c.exterior = c.exterior or exterior[t]
                for _, v in ipairs(neighbors[t]) do
                    if not of[v] then
                        of[v] = ci
                        c[#c + 1] = v
                    end
                end
                j = j + 1
            end
        end
    end
    local function seed(p, tag)
        local found
        for i, t in ipairs(ids) do
            local a, b, c = nodes[t[1]], nodes[t[2]], nodes[t[3]]
            local function cross(u, v)
                return (v[1] - u[1]) * (p[2] - u[2]) - (v[2] - u[2]) * (p[1] - u[1])
            end
            if
                cross(a, b) >= -tol * tol
                and cross(b, c) >= -tol * tol
                and cross(c, a) >= -tol * tol
            then
                check(not found or found == of[i], "label lies on a physical interface")
                found = of[i]
            end
        end
        check(found, "label outside native geometry")
        check(not cells[found].tag, "multiple labels or holes in one connected cell")
        cells[found].tag = tag
    end
    local default
    for i, p in ipairs(d.labels) do
        if p[6] == 2 then
            check(not default, "multiple default labels")
            default = i
        end
        seed(p, i)
    end
    for _, p in ipairs(d.holes) do
        seed(p, 0)
    end
    local kept, labels, used = {}, {}, {}
    for i, t in ipairs(ids) do
        local c = cells[of[i]]
        local tag = c.tag
        if tag == nil and not c.exterior then
            tag = default
            check(tag, "unlabelled enclosed native cell")
        end
        if tag and tag > 0 then
            check(not c.exterior, "native material cell is not enclosed by declared segments")
            kept[#kept + 1] = t
            labels[#labels + 1] = tag
            for _, id in ipairs(t) do
                used[id] = true
            end
        end
    end
    check(#kept > 0, "no active dielectric elements")
    local compact, map = {}, {}
    for i, p in ipairs(nodes) do
        if used[i] then
            map[i] = #compact + 1
            compact[#compact + 1] = p
        end
    end
    for _, t in ipairs(kept) do
        for j = 1, 3 do
            t[j] = map[t[j]]
        end
    end
    for _, p in ipairs(f.periodic_pairs) do
        check(map[p[1]] and map[p[2]], "periodic boundary outside active cells")
        p[1], p[2] = map[p[1]], map[p[2]]
    end
    f.periodic_nodes = {}
    for _, p in ipairs(f.periodic_pairs) do
        f.periodic_nodes[p[1]], f.periodic_nodes[p[2]] = true, true
    end
    f.native_labels = labels
    return compact, kept, stats
end
--- Native qp is C/m: a transverse line or a ring with Q=2*pi*r*qp.
-- Explicit constraints win over loads. Conflicting prescribed values fail;
-- silent native precedence is intentionally not used to hide contradictory data.
function N.points(m, nodes, fixed, owners, floating, fix)
    local f = m.fe
    local d = f.native
    if not d then
        return
    end
    f.point_loads = {}
    local tol = 1e-8 * math.max(f.bounds[3] - f.bounds[1], f.bounds[4] - f.bounds[2])
    if f.axisymmetric then
        for _, s in ipairs(d.lines) do
            if s.a[1] == 0 and s.b[1] == 0 then
                for id, p in ipairs(nodes) do
                    if N.on(p, s.a, s.b, tol) then
                        if s.conductor > 0 then
                            local c = d.conductors[s.conductor]
                            fix(id, c.vc, c.conductorname, true)
                        elseif s.boundary > 0 then
                            fix(id, d.boundaries[s.boundary].vs, "outer boundary", true)
                        end
                    end
                end
            end
        end
    end
    for _, p in ipairs(d.points) do
        if p[3] > 0 or p[5] > 0 then
            local id
            for i, q in ipairs(nodes) do
                if math.abs(p[1] - q[1]) < tol and math.abs(p[2] - q[2]) < tol then
                    check(not id, "ambiguous native point")
                    id = i
                end
            end
            check(id, "native point is outside active mesh")
            if p[5] > 0 then
                local c = d.conductors[p[5]]
                fix(id, c.vc, c.conductorname, true)
            end
            if p[3] > 0 then
                local prop = d.pointprops[p[3]]
                if prop.qp == 0 then
                    fix(id, prop.vp, "outer boundary", true)
                elseif floating[id] then
                    check(false, "point charge on a floating conductor is unsupported")
                elseif fixed[id] == nil then
                    f.point_loads[id] = (f.axisymmetric and 2 * math.pi * m.unit * p[1] or 1)
                        * prop.qp
                end
            end
        end
    end
end
return N
