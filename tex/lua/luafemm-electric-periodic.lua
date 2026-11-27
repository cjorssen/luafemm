-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Matching translated straight exterior paths, including mesher subdivisions.
-- Pairs have identical oriented vectors and opposite exterior normals. In an
-- axisymmetric model only axial translations preserve the radial measure.
-- @module luafemm-electric-periodic
local P = {}
local mesh = require("luafemm-mesh")
local copy = require("luafemm-physical").copy
local function check(v, s)
    assert(v, "luafemm electric periodic: " .. s)
end
function P.generate(f, constraints)
    local groups, order = {}, {}
    for _, p in ipairs(f.paths) do
        if p.condition == "periodic" or p.condition == "antiperiodic" then
            check(
                not p.closed and #p.points == 2,
                "paired boundary must be one straight open segment"
            )
            check(
                type(p.pair) == "string" and p.pair:match("%S"),
                "paired boundary needs a pair name"
            )
            if not groups[p.pair] then
                groups[p.pair] = {}
                order[#order + 1] = p.pair
            end
            local g = groups[p.pair]
            g[#g + 1] = p
        end
    end
    local span = math.max(f.bounds[3] - f.bounds[1], f.bounds[4] - f.bounds[2])
    local tol = math.max(1e-9 * span, f.options.mesh_tolerance * 8)
    local pairspec = {}
    local function side(p)
        local a, b = p.points[1], p.points[2]
        for _, spec in ipairs({ { 1, 1, -1 }, { 1, 3, 1 }, { 2, 2, -1 }, { 2, 4, 1 } }) do
            if
                math.abs(a[spec[1]] - f.bounds[spec[2]]) < tol
                and math.abs(b[spec[1]] - f.bounds[spec[2]]) < tol
            then
                return spec[1], spec[3]
            end
        end
    end
    for _, name in ipairs(order) do
        local g = groups[name]
        check(#g == 2, "a periodic pair needs exactly two paths")
        local a, b = g[1], g[2]
        check(a.condition == b.condition, "pair conditions disagree")
        local sa, na = side(a)
        local sb, nb = side(b)
        check(sa and sa == sb and na == -nb, "paired paths must lie on opposite exterior sides")
        for j = 1, 2 do
            check(
                math.abs((a.points[2][j] - a.points[1][j]) - (b.points[2][j] - b.points[1][j]))
                    < tol,
                "paired paths need identical oriented vectors"
            )
        end
        check(not f.axisymmetric or sa == 2, "axisymmetric periodicity requires axial translation")
        pairspec[#pairspec + 1] = { a, b, a.condition == "periodic" and 1 or -1 }
    end
    local function coordinate(p, path)
        local a, b = path.points[1], path.points[2]
        local dx, dy = b[1] - a[1], b[2] - a[2]
        local l = dx * dx + dy * dy
        check(l > 0, "zero-length periodic path")
        local t = ((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / l
        if
            t >= -tol / math.sqrt(l)
            and t <= 1 + tol / math.sqrt(l)
            and math.abs(dx * (p[2] - a[2]) - dy * (p[1] - a[1])) <= tol * math.sqrt(l)
        then
            return t
        end
    end
    local work = copy(constraints)
    for pass = 1, 12 do
        local nodes, ids, stats = mesh.generate(f.options, work)
        local pairs, marked, missing = {}, {}, {}
        for _, g in ipairs(pairspec) do
            local sides = { {}, {} }
            for i, p in ipairs(nodes) do
                for side = 1, 2 do
                    local t = coordinate(p, g[side])
                    if t then
                        sides[side][#sides[side] + 1] = { t, i }
                    end
                end
            end
            for side = 1, 2 do
                local other = 3 - side
                for _, u in ipairs(sides[side]) do
                    local found
                    for _, v in ipairs(sides[other]) do
                        if math.abs(u[1] - v[1]) < tol / span then
                            found = v[2]
                            break
                        end
                    end
                    if found then
                        if side == 1 then
                            pairs[#pairs + 1] = { u[2], found, g[3] }
                            marked[u[2]], marked[found] = true, true
                        end
                    else
                        local a, b = g[other].points[1], g[other].points[2]
                        local q = { a[1] + u[1] * (b[1] - a[1]), a[2] + u[1] * (b[2] - a[2]) }
                        missing[#missing + 1] = {
                            contours = { { a, q }, { q, b } },
                            mesh_size = span * 2,
                            fill_rule = "nonzero",
                        }
                    end
                end
            end
        end
        if #missing == 0 then
            f.periodic_pairs, f.periodic_nodes = pairs, marked
            return nodes, ids, stats
        end
        for _, r in ipairs(missing) do
            work[#work + 1] = r
        end
    end
    error(
        "luafemm electric periodic: paired mesh synchronization did not converge; simplify boundary refinements"
    )
end
return P
