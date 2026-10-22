-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Release sensitivity study: mesh error, exterior truncation and nonlinear force.
-- Run separately from the fast suite: make static-convergence.
dofile("tests/bootstrap.lua")
local M, F = require("luafemm"), require("luafemm-frontiers")
local post = require("luafemm-post")
require("lfs").mkdir("build/convergence")
local out = assert(io.open("build/convergence/static.csv", "w"))
out:write("study,spacing,radius,nodes,triangles,value,error\n")
local function record(study, h, radius, m, value, err)
    local row = string.format(
        "%s,%g,%g,%d,%d,%.12g,%.12g\n",
        study,
        h,
        radius,
        #m.nodes,
        #m.triangles,
        value,
        err
    )
    out:write(row)
    out:flush()
    io.write(row)
end

-- A=sin(pi*x), B=(0,-pi*cos(pi*x)). Three symmetric quadrature points
-- integrate the squared field error, including variation inside each triangle.
local previous
for _, h in ipairs({ 0.16, 0.08, 0.04 }) do
    local m = M.new({
        xmin = 0,
        xmax = 1,
        ymin = 0,
        ymax = 1,
        unit = 1,
        h = h,
        mesher = "delaunay",
        tolerance = 1e-10,
        source = function(x)
            return math.pi ^ 2 / M.mu0 * math.sin(math.pi * x)
        end,
    })
    M.boundary(m, "all", { type = "neumann" })
    F.path(m, "left", { { 0, 0 }, { 0, 1 } }, false, { type = "periodic" })
    F.path(m, "right", { { 1, 0 }, { 1, 1 } }, false, { type = "periodic" })
    F.pair(m, "left", "right", { x = 1, sign = -1 })
    M.solve(m)
    local err = 0
    for _, e in ipairs(m.triangles) do
        local sum = 0
        for _, id in ipairs(e.ids) do
            sum = sum + m.nodes[id][1]
        end
        for _, id in ipairs(e.ids) do
            local x = sum / 6 + m.nodes[id][1] / 2
            err = err + e.area / 3 * (e.bx ^ 2 + (e.by + math.pi * math.cos(math.pi * x)) ^ 2)
        end
    end
    err = math.sqrt(err) / (math.pi / math.sqrt(2))
    record("antiperiodic", h, 0, m, err, err)
    if previous then
        assert(previous / err > 1.6, "Field refinement lost first-order convergence")
    end
    previous = err
end

local function circle(r, n, cx)
    local p = {}
    for i = 1, n do
        local a = 2 * math.pi * (i - 1) / n
        p[i] = { cx + r * math.cos(a), r * math.sin(a) }
    end
    return p
end
-- An off-centre uniformly magnetised disk has Bx=1/2 T internally in
-- infinite air. Offsetting it excites exterior modes beyond the Robin n=1
-- approximation. Keep the magnet polygon/local spacing fixed while moving R.
for _, case in ipairs({ { 0.14, 0.9 }, { 0.14, 1.5 }, { 0.14, 2.5 }, { 0.07, 1.5 }, { 0.035, 1.5 } }) do
    local h, radius = case[1], case[2]
    local m = M.new({
        xmin = -radius,
        xmax = radius,
        ymin = -radius,
        ymax = radius,
        unit = 1,
        h = h,
        mesher = "delaunay",
        tolerance = 1e-10,
    })
    local outer = circle(radius, math.ceil(2 * math.pi * radius / (h * 6 / 7)), 0)
    F.domain(m, { outer })
    F.path(m, "open", outer, true, { type = "open", radius = radius })
    M.material(m, "magnet", { mur = 1, coercivity = 1 / M.mu0 })
    M.region(m, "magnet", circle(0.25, 64, 0.4), 0, 0, h / 4)
    M.solve(m)
    local bx, area = 0, 0
    for _, e in ipairs(m.triangles) do
        if e.material == "magnet" then
            bx = bx + e.bx * e.area
            area = area + e.area
        end
    end
    bx = bx / area
    record(h == 0.14 and "open-radius" or "open-mesh", h, radius, m, bx, math.abs(bx - 0.5) / 0.5)
end

-- Same nonlinear armature and weighting convention as the xfemm oracle.
-- Successive changes are sensitivity evidence, not an exact-force error bound.
local make = dofile("tests/u_model.lua")
previous = nil
for _, h in ipairs({ 4, 2, 1.5 }) do
    local m = make(h, 75, 3000)
    m.options.depth = 20
    m.options.interpolation = "femm"
    m.options.mesher = "delaunay"
    for id, p in pairs(m.materials) do
        p.interpolation = "femm"
        M.material(m, id, p)
    end
    post.select(
        m,
        "armature",
        { material = "iron", contours = { { { -41, 32 }, { 41, 32 }, { 41, 49 }, { -41, 49 } } } }
    )
    local r = post.measure(m, { name = "force", selection = "armature", quantity = "force" })
    local change = previous and math.abs(r.Fy - previous) / math.abs(r.Fy) or 0
    record("armature-force", h, 75, m, r.Fy, change)
    previous = r.Fy
end
out:close()
print("PASS: static convergence study; see build/convergence/static.csv")
