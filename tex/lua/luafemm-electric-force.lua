-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Maxwell traction on a closed contour in homogeneous source-free isotropic medium.
-- Partition at every crossed triangle face. Ordinary P1 gradients integrate
-- exactly; a source/correction field uses bounded adaptive total-stress rules.
-- No smoothing or path-independence is implied on a finite mesh. Return total
-- force (N), planar torque (N m), and planar per-length quantities when defined.
-- @module luafemm-electric-force
local M = {}
local E = require("luafemm-electrostatic")
local F = require("luafemm-electric-fe")
local V = require("luafemm-electric-fe-view")
local copy = require("luafemm-physical").copy
local function check(v, s)
    assert(v, "luafemm electric force: " .. s)
end
-- Reject contact as well as proper crossings: touching an electrode at a
-- contour vertex or following a charged sheet has no unambiguous dielectric
-- trace. The tolerance is geometric roundoff, not a user-visible smoothing.
local function touches(a, b, c, d, tol)
    local function distance(p, u, v)
        local dx, dy = v[1] - u[1], v[2] - u[2]
        local l = dx * dx + dy * dy
        local t = l > 0 and math.max(0, math.min(1, ((p[1] - u[1]) * dx + (p[2] - u[2]) * dy) / l))
            or 0
        return (p[1] - u[1] - t * dx) ^ 2 + (p[2] - u[2] - t * dy) ^ 2
    end
    local dx, dy = b[1] - a[1], b[2] - a[2]
    local ux, uy = d[1] - c[1], d[2] - c[2]
    local det = dx * uy - dy * ux
    if det ~= 0 then
        local t = ((c[1] - a[1]) * uy - (c[2] - a[2]) * ux) / det
        local u = ((c[1] - a[1]) * dy - (c[2] - a[2]) * dx) / det
        if t >= 0 and t <= 1 and u >= 0 and u <= 1 then
            return true
        end
    end
    return math.min(distance(a, c, d), distance(b, c, d), distance(c, a, b), distance(d, a, b))
        <= tol * tol
end
function M.calculate(m, points, o)
    o = o or {}
    check(m.fe, "finite-element model required")
    E.freeze(m)
    points = copy(points)
    if
        #points > 1
        and points[1][1] == points[#points][1]
        and points[1][2] == points[#points][2]
    then
        table.remove(points)
    end
    require("luafemm-mesh").validate(points)
    local area = 0
    for i, p in ipairs(points) do
        local q = points[i % #points + 1]
        area = area + p[1] * q[2] - q[1] * p[2]
    end
    local orientation = area > 0 and 1 or -1
    local ox, oy =
        E.number(o.origin_x or 0, "force origin x"), E.number(o.origin_y or 0, "force origin y")
    local fx, fy, torque, epsilon = 0, 0, 0, nil
    local integrate
    if m.hybrid then
        local tolerance = E.positive(o.quadrature_tolerance or 1e-6, "force quadrature tolerance")
        check(tolerance < 1, "force quadrature tolerance must be less than one")
        integrate = require("luafemm-electric-quadrature").new(
            tolerance,
            E.integer(o.quadrature_depth or 10, "maximum force quadrature depth", 1, 16),
            E.integer(o.quadrature_cells or 250000, "maximum force quadrature cells", 1, 2000000),
            "source force"
        )
    end
    for i, a in ipairs(points) do
        local b = points[i % #points + 1]
        -- Closing along the regular axis has zero revolved area, even where
        -- the enclosed body meets the axis. It carries no Maxwell traction.
        if not (m.fe.axisymmetric and a[1] == 0 and b[1] == 0) then
            check(
                not V.crosses_interface(m, a, b),
                "contour crosses an electrode, material interface or charged layer"
            )
            local bounds = m.fe.bounds
            local tol = 1e-10 * math.max(bounds[3] - bounds[1], bounds[4] - bounds[2])
            if m.hybrid then
                for _, s in ipairs(m.sources) do
                    local p = { s.x, s.y }
                    check(not touches(a, b, p, p, tol), "force contour touches a singular source")
                end
            end
            for _, edge in ipairs(m.fe.edges) do
                local physical = edge.electrode or edge.loads ~= 0
                if #edge.triangles == 2 then
                    local u, v = edge.triangles[1], edge.triangles[2]
                    physical = physical
                        or m.fe.interface_pairs[math.min(u, v) .. ":" .. math.max(u, v)]
                end
                if physical and not edge.axis then
                    check(
                        not touches(a, b, m.fe.nodes[edge[1]], m.fe.nodes[edge[2]], tol),
                        "force contour touches an electrode, material interface or charged layer"
                    )
                end
            end
            local cuts = V.cuts(m, a, b)
            for j = 2, #cuts do
                local lo, hi = cuts[j - 1], cuts[j]
                local mid = (lo + hi) / 2
                local x, y = a[1] + mid * (b[1] - a[1]), a[2] + mid * (b[2] - a[2])
                local id = F.locate(m, x, y)
                local t = id and m.fe.triangles[id]
                check(
                    t and not t.electrode and t.density == 0,
                    "contour must remain in source-free dielectric"
                )
                local tensor = t.tensor
                check(
                    math.abs(tensor.xy) < 1e-12 * math.max(tensor.xx, tensor.yy)
                        and math.abs(tensor.xx - tensor.yy)
                            < 1e-12 * math.max(tensor.xx, tensor.yy),
                    "force contour medium must be isotropic"
                )
                local ep = E.epsilon0 * tensor.xx
                check(
                    not epsilon or math.abs(ep / epsilon - 1) < 1e-12,
                    "force contour medium must be homogeneous"
                )
                epsilon = ep
                if integrate then
                    -- Integrate the actual nonlinear total stress. Drawing masks
                    -- are irrelevant; only contact with the physical singularity
                    -- is excluded. Partitioning preserves one-sided FE gradients.
                    local nx = orientation * (b[2] - a[2]) * m.unit
                    local ny = -orientation * (b[1] - a[1]) * m.unit
                    local value = integrate(function(tau)
                        local px, py = a[1] + tau * (b[1] - a[1]), a[2] + tau * (b[2] - a[2])
                        local s = assert(F.sample(m, px, py, id))
                        local xx = ep * (s.Ex * s.Ex - 0.5 * s.E * s.E)
                        local yy = ep * (s.Ey * s.Ey - 0.5 * s.E * s.E)
                        local xy = ep * s.Ex * s.Ey
                        local tx, ty = xx * nx + xy * ny, xy * nx + yy * ny
                        if m.fe.axisymmetric then
                            return { 0, 2 * math.pi * px * m.unit * ty, 0 }
                        end
                        return { tx, ty, (px - ox) * m.unit * ty - (py - oy) * m.unit * tx }
                    end, lo, hi)
                    fx, fy, torque = fx + value[1], fy + value[2], torque + value[3]
                else
                    -- n ds for a counterclockwise contour; coordinates and moment arm in SI.
                    local nx = orientation * (b[2] - a[2]) * (hi - lo) * m.unit
                    local ny = -orientation * (b[1] - a[1]) * (hi - lo) * m.unit
                    local xx = ep * (t.Ex * t.Ex - 0.5 * t.E * t.E)
                    local yy = ep * (t.Ey * t.Ey - 0.5 * t.E * t.E)
                    local xy = ep * t.Ex * t.Ey
                    local tx, ty = xx * nx + xy * ny, xy * nx + yy * ny
                    if m.fe.axisymmetric then
                        fy = fy + 2 * math.pi * x * m.unit * ty
                    else
                        fx, fy = fx + tx, fy + ty
                        torque = torque + (x - ox) * m.unit * ty - (y - oy) * m.unit * tx
                    end
                end
            end
        end
    end
    if m.fe.axisymmetric then
        return { Fz = fy, model = E.snapshot(m), points = points }
    end
    return {
        Fx = fx * m.fe.depth,
        Fy = fy * m.fe.depth,
        torque = torque * m.fe.depth,
        ["Fx per length"] = fx,
        ["Fy per length"] = fy,
        ["torque per length"] = torque,
        model = E.snapshot(m),
        points = points,
    }
end
return M
