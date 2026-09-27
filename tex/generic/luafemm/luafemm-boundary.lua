-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Exterior boundary data and P1 edge assembly, independent of TeX.
-- With outward normal n, the positive tangent is t=(-ny,nx).
-- Neumann prescribes H.t=g; Robin prescribes H.t=c*Az+g, c>=0.
-- Boundary integrals use SI coordinates, also for affine data g.
-- @module luafemm-boundary
local B = {}
local sides = { "left", "right", "bottom", "top" }
local valid_sides = { left = true, right = true, bottom = true, top = true }
local kinds = { dirichlet = true, neumann = true, robin = true }
local function check(value, message)
    assert(value, "luafemm boundary: " .. message)
end
local function finite(value)
    return type(value) == "number" and value == value and math.abs(value) < math.huge
end

--- Initialise the four sides with homogeneous Dirichlet data.
-- @tparam table m Newly constructed model with validated domain options.
function B.configure(m)
    check(not m.options.boundary or type(m.options.boundary) == "function", "invalid callback")
    m.boundaries = {}
    for _, side in ipairs(sides) do
        m.boundaries[side] = { type = "dirichlet", value = 0, dx = 0, dy = 0, coefficient = 0 }
    end
end

--- Replace one complete side condition before meshing (or all four sides).
-- @tparam table m Unmeshed model.
-- @tparam string side left, right, bottom, top or all.
-- @tparam[opt] table p type, value, dx, dy, coefficient; omitted numbers are zero.
-- Dirichlet: Az=value+dx*x+dy*y. Neumann/Robin: g=value+dx*x+dy*y.
function B.set(m, side, p)
    check(m and m.boundaries, "create a problem before declaring boundary conditions")
    check(not m.nodes, "declare boundary conditions before meshing")
    check(not m.options.boundary, "cannot combine side conditions with the boundary callback")
    check(side == "all" or valid_sides[side], "unknown side " .. tostring(side))
    p = p or {}
    check(type(p) == "table", "expected a condition table")
    local condition = { type = p.type == nil and "dirichlet" or p.type }
    check(kinds[condition.type], "unknown type " .. tostring(condition.type))
    for _, key in ipairs({ "value", "dx", "dy", "coefficient" }) do
        condition[key] = p[key] == nil and 0 or p[key]
        check(finite(condition[key]), "non-finite " .. key)
    end
    for key in pairs(p) do
        check(condition[key] ~= nil, "unknown option " .. tostring(key))
    end
    check(condition.coefficient >= 0, "Robin coefficient must be nonnegative")
    check(
        condition.type == "robin" or condition.coefficient == 0,
        "coefficient requires type=robin"
    )
    for _, name in ipairs(sides) do
        if side == "all" or name == side then
            local own = {}
            for key, value in pairs(condition) do
                own[key] = value
            end
            m.boundaries[name] = own
        end
    end
end

local function data(condition, point)
    local value = condition.value + condition.dx * point[1] + condition.dy * point[2]
    check(finite(value), "non-finite evaluated boundary value")
    return value
end

--- Reconstruct exterior edges and prescribed nodal values without mutating m.
-- The exterior tag is geometric, not a list of eliminated degrees of freedom.
-- The deterministic edge order is shared by fresh and cached meshes.
-- @tparam table m Model with current side declarations.
-- @tparam table nodes SI vertex coordinates.
-- @tparam table elements Counterclockwise triples of vertex indices.
-- @tparam table exterior Set of geometric exterior vertex indices.
-- @return Fixed-value map, weak boundary edges, optional gauge-node index.
function B.prepare(m, nodes, elements, exterior)
    local o = m.options
    local x0, x1, y0, y1 = o.xmin * o.unit, o.xmax * o.unit, o.ymin * o.unit, o.ymax * o.unit
    local epsilon = math.max(x1 - x0, y1 - y0) * 1e-10
    local candidates = {}
    for _, ids in ipairs(elements) do
        for j = 1, 3 do
            local a, b = ids[j], ids[j % 3 + 1]
            if exterior[a] and exterior[b] then
                if a > b then
                    a, b = b, a
                end
                local key = a .. ":" .. b
                local edge = candidates[key] or { ids = { a, b }, count = 0 }
                edge.count = edge.count + 1
                candidates[key] = edge
            end
        end
    end
    local edges = {}
    for _, edge in pairs(candidates) do
        if edge.count == 1 then
            edges[#edges + 1] = edge
        end
    end
    table.sort(edges, function(a, b)
        return a.ids[1] < b.ids[1] or (a.ids[1] == b.ids[1] and a.ids[2] < b.ids[2])
    end)
    local fixed, weak, anchored = {}, {}, false
    for _, edge in ipairs(edges) do
        local p, q = nodes[edge.ids[1]], nodes[edge.ids[2]]
        local function distance(axis, value)
            return math.max(math.abs(p[axis] - value), math.abs(q[axis] - value))
        end
        -- Select the closest whole side, not the first side within tolerance:
        -- on a very thin rectangle both opposite sides may pass that test.
        local distances = { distance(1, x0), distance(1, x1), distance(2, y0), distance(2, y1) }
        local best = 1
        for i = 2, 4 do
            if distances[i] < distances[best] then
                best = i
            end
        end
        check(distances[best] <= epsilon, "exposed edge is not on a domain side")
        local side = sides[best]
        local condition = m.boundaries[side]
        if o.boundary or condition.type == "dirichlet" then
            for _, id in ipairs(edge.ids) do
                local value
                if o.boundary then
                    value = fixed[id]
                    if value == nil then
                        value = o.boundary(nodes[id][1], nodes[id][2])
                    end
                else
                    value = data(condition, nodes[id])
                end
                check(finite(value), "non-finite prescribed potential")
                local old = fixed[id]
                check(
                    old == nil
                        or math.abs(old - value)
                            <= 1e-10 * math.max(1e-12, math.abs(old), math.abs(value)),
                    "conflicting Dirichlet values at a corner"
                )
                fixed[id] = value
            end
            anchored = true
        else
            edge.side = side
            edge.length = math.sqrt((q[1] - p[1]) ^ 2 + (q[2] - p[2]) ^ 2)
            edge.g = { data(condition, p), data(condition, q) }
            edge.coefficient = condition.coefficient
            edge.count = nil
            weak[#weak + 1] = edge
            anchored = anchored or condition.coefficient > 0
        end
    end
    local gauge
    if not anchored then
        -- A single reference fixes only the additive constant. The source
        -- compatibility check below must succeed before this system is solved.
        gauge = 1
        fixed[gauge] = 0
    end
    return fixed, weak, gauge
end

--- Reject incompatible all-Neumann loads before eliminating the gauge row.
-- Ampere's law requires integral(Jz dS)=integral(H.t ds).
-- Magnetisation does not enter this balance: curls of a constant test vanish.
-- @tparam table m Prepared triangles, boundary_edges and optional gauge_node.
function B.compatibility(m)
    if not m.gauge_node then
        return
    end
    local balance, scale = 0, 0
    for _, e in ipairs(m.triangles) do
        local load = e.current * e.area
        balance, scale = balance + load, scale + math.abs(load)
    end
    for _, edge in ipairs(m.boundary_edges) do
        balance = balance - edge.length * (edge.g[1] + edge.g[2]) / 2
        scale = scale + edge.length * (math.abs(edge.g[1]) + math.abs(edge.g[2])) / 2
    end
    check(
        math.abs(balance) <= 1e-10 * math.max(scale, 1e-30),
        "incompatible Neumann data: net current must equal the counterclockwise integral of H.t"
    )
end

--- Add exact P1 edge integrals to a mu0-scaled residual and tangent.
-- An affine g uses the consistent edge mass matrix, not endpoint lumping.
-- @tparam table m Prepared model with nodal degrees of freedom and edge data.
-- @tparam table A Current nodal potential in T m.
-- @tparam table R Residual vector, updated in place.
-- @tparam table K Sparse tangent rows, updated in place; nil for residual only.
-- @tparam number mu0 Vacuum permeability used to scale the entire system.
function B.assemble(m, A, R, K, mu0)
    for _, edge in ipairs(m.boundary_edges) do
        local weight = mu0 * edge.length / 6
        for i = 1, 2 do
            local row = m.nodes[edge.ids[i]].dof
            if row then
                for j = 1, 2 do
                    local mass = weight * (i == j and 2 or 1)
                    R[row] = R[row] + mass * (edge.coefficient * A[edge.ids[j]] + edge.g[j])
                    local col = m.nodes[edge.ids[j]].dof
                    if K and col then
                        K[row][col] = (K[row][col] or 0) + mass * edge.coefficient
                    end
                end
            end
        end
    end
end

return B
