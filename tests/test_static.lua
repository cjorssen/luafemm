-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Analytic tensor, periodic, open-boundary and circuit regression problems.
dofile("tests/bootstrap.lua")
local M = require("luafemm")
local F = require("luafemm-frontiers")
local S = require("luafemm-static")
local function near(a, b, t)
    assert(math.abs(a - b) < t, string.format("%.16g ~= %.16g", a, b))
end
local function fails(fn, word)
    local ok, err = pcall(fn)
    assert(not ok and tostring(err):find(word, 1, true), tostring(err))
end
local function model(o)
    o = o or {}
    o.xmin = o.xmin or 0
    o.xmax = o.xmax or 1
    o.ymin = o.ymin or 0
    o.ymax = o.ymax or 1
    o.unit = o.unit or 1
    o.h = o.h or 0.18
    o.mesher = "delaunay"
    o.depth = o.depth or 0.2
    o.tolerance = 1e-10
    return M.new(o)
end
-- An oblique polygon and a full rotated tensor reproduce any affine Az.
do
    local m = model()
    local poly = { { 0.1, 0.1 }, { 0.9, 0.2 }, { 0.8, 0.9 }, { 0.2, 0.8 } }
    F.domain(m, { poly })
    F.path(m, "outer", poly, true, { value = 0.3, dx = 0.4, dy = -0.2 })
    M.material(m, "tensor", { mur_x = 3, mur_y = 11, material_angle = 37 })
    M.region(m, "tensor", poly)
    M.solve(m)
    for i, p in ipairs(m.nodes) do
        near(m.A[i], 0.3 + 0.4 * p[1] - 0.2 * p[2], 1e-10)
    end
    local hx, hy = S.law(m.materials.tensor, -0.2, -0.4)
    for _, e in ipairs(m.triangles) do
        near(e.hx, hx / M.mu0, 1e-4)
        near(e.hy, hy / M.mu0, 1e-4)
    end
end
print("Oblique domain / rotated tensor affine patch OK")
-- Linear sheet homogenisation has the exact arithmetic/harmonic limits.
do
    local m = model()
    M.material(m, "sheets", { mur = 100, lamination = "parallel x", lamination_fill = 0.8 })
    near(m.materials.sheets.mur_x, 80.2, 1e-12)
    near(m.materials.sheets.mur_y, 1 / (0.8 / 100 + 0.2), 1e-12)
    fails(function()
        M.material(m, "bad", { mur = 100, mur_y = 5, bh = { { 0, 0 }, { 1, 100 } } })
    end, "nonlinear anisotropy")
end
-- Different boundary subdivisions are synchronised, both periodic signs.
for _, sign in ipairs({ 1, -1 }) do
    local m = model()
    M.boundary(m, "all", { type = "neumann" })
    F.path(m, "left", { { 0, 0 }, { 0, 0.27 }, { 0, 1 } }, false, { type = "periodic" })
    F.path(m, "right", { { 1, 0 }, { 1, 0.63 }, { 1, 1 } }, false, { type = "periodic" })
    F.pair(m, "left", "right", { x = 1, sign = sign })
    if sign == 1 then
        M.boundary(m, "bottom", { value = 0 })
        M.boundary(m, "top", { value = 0.2 })
    else
        -- A=x-1/2, with homogeneous Neumann data at horizontal edges.
        -- A constant imposed J is incompatible with this analytic field;
        -- set a closed internal Dirichlet contour to provide excitation.
        F.path(
            m,
            "drive",
            { { 0.2, 0.3 }, { 0.4, 0.3 }, { 0.4, 0.7 }, { 0.2, 0.7 } },
            true,
            { value = -0.5, dx = 1 }
        )
    end
    M.solve(m)
    assert(m.mesh_stats.paired_vertices > 0)
    for _, pair in ipairs(m.periodic_nodes) do
        near(m.A[pair[1]], sign * m.A[pair[2]], 1e-12)
    end
    local saved_tex = tex
    local output = {}
    tex = {
        sprint = function(s)
            output[#output + 1] = s
        end,
    }
    m.frame = { 1, 0, 0, 1, 0, 0 }
    require("luafemm-tex").tikz_repeated_field(m, 6, 2, 1)
    tex = saved_tex
    assert(#output > 0 and #output % 2 == 0)
    if sign == 1 then
        for i, p in ipairs(m.nodes) do
            near(m.A[i], 0.2 * p[2], 1e-9)
        end
    end
end
print("Periodic and antiperiodic strong nodal constraints / mesh synchronisation OK")
-- Signed-periodic manufactured Poisson field: A=sin(pi*x), J=pi^2/mu0*sin(pi*x).
do
    local m = model({
        h = 0.08,
        source = function(x)
            return math.pi ^ 2 / M.mu0 * math.sin(math.pi * x)
        end,
    })
    M.boundary(m, "all", { type = "neumann" })
    F.path(m, "left", { { 0, 0 }, { 0, 1 } }, false, { type = "periodic" })
    F.path(m, "right", { { 1, 0 }, { 1, 1 } }, false, { type = "periodic" })
    F.pair(m, "left", "right", { x = 1, sign = -1 })
    M.solve(m)
    for i, p in ipairs(m.nodes) do
        near(m.A[i], math.sin(math.pi * p[1]), 0.015)
    end
end
print("Antiperiodic manufactured field OK")
-- Circular first-multipole exterior: uniform magnet disk, analytic central B.
local function disk(radius, n)
    local p = {}
    for i = 1, n do
        local a = 2 * math.pi * (i - 1) / n
        p[i] = { radius * math.cos(a), radius * math.sin(a) }
    end
    return p
end
do
    local m = model({ xmin = -2, xmax = 2, ymin = -2, ymax = 2, h = 0.18 })
    local outer = disk(2, 80)
    F.domain(m, { outer })
    F.path(m, "infinity", outer, true, { type = "open", radius = 2 })
    M.material(m, "magnet", { mur = 1, coercivity = 1 / M.mu0 })
    M.region(m, "magnet", disk(0.5, 48))
    M.solve(m)
    local bx, by = M.sample(m, 0, 0)
    near(bx, 0.5, 0.015)
    near(by, 0, 0.01)
end
print("Circular open boundary / magnet dipole analytic field OK")
-- Three-point source averaging converges; radial and angular presets agree.
do
    local m = model()
    M.material(m, "pm", { coercivity = 1000 })
    local r = M.region(m, "pm", { { 0.2, 0.2 }, { 0.9, 0.2 }, { 0.2, 0.9 } })
    S.region(r, { magnetization = "radial" })
    local a, b = S.direction(r, { { 0.2, 0.2 }, { 0.9, 0.2 }, { 0.2, 0.9 } })
    near(a, b, 1e-14)
    assert(a * a + b * b < 1)
    S.region(r, { magnetization = "angular", order = 1 })
    local c, d = S.direction(r, { { 0.2, 0.2 }, { 0.9, 0.2 }, { 0.2, 0.9 } })
    near(a, c, 1e-14)
    near(b, d, 1e-14)
end
-- A closed two-leg coil: ampere turns, ohmic law, energy/linkage identity.
local function coil(kind, drive)
    local m = model({ h = 0.12 })
    M.material(m, "copper", { conductivity = 5.8e7 })
    S.circuit(m, "coil", {
        type = kind,
        drive = drive,
        current = drive == "current" and 2 or 0,
        voltage = drive == "voltage" and 0.01 or 0,
    })
    for k, x in ipairs({ 0.15, 0.65 }) do
        local r =
            M.region(m, "copper", { { x, 0.2 }, { x + 0.2, 0.2 }, { x + 0.2, 0.8 }, { x, 0.8 } })
        S.region(r, { circuit = "coil", turns = k == 1 and 20 or -20 })
    end
    M.solve(m)
    local c = S.result(m, "coil")
    local R = 2 * 400 * 0.2 / (5.8e7 * 0.12)
    near(c.resistance, R, 1e-16)
    near(c.voltage, c.current * R, 1e-15)
    local energy = 0
    for _, e in ipairs(m.triangles) do
        energy = energy + (e.bx * e.hx + e.by * e.hy) / 2 * e.area * 0.2
    end
    near(energy, c.linkage * c.current / 2, math.max(1e-12, energy * 1e-8))
    assert(c.inductance > 0)
    return m
end
coil("series", "current")
coil("series", "voltage")
coil("parallel", "current")
print("Named DC circuits: signed turns, voltage/current drive, linkage/energy identity OK")
print("Static extension tests passed")

-- Round trips exercise the documented representable subset of FEMM records.
local E = require("luafemm-fem")
local X = require("luafemm-static-exchange")
local C = require("luafemm-cache")
require("lfs").mkdir("build/static")
do
    local m = coil("series", "current")
    local doc = E.document(m)
    assert(#doc.circuits == 1 and doc.circuits[1].circuitname == "coil")
    local other = E.model(doc, { h = 0.12, unit = 1 })
    M.solve(other)
    near(S.result(other, "coil").linkage, S.result(m, "coil").linkage, 1e-11)
    M.export_ans(m, "build/static/coil.ans")
    for _, sign in ipairs({ 1, -1 }) do
        local q = model({ unit = 0.001 })
        M.boundary(q, "all", { type = "neumann" })
        M.boundary(q, "bottom", { value = 0 })
        if sign == 1 then
            M.boundary(q, "top", { value = 0.0002 })
        end
        F.path(q, "a", { { 0, 0 }, { 0, 0.27 }, { 0, 1 } }, false, { type = "periodic" })
        F.path(q, "b", { { 1, 0 }, { 1, 0.63 }, { 1, 1 } }, false, { type = "periodic" })
        F.pair(q, "a", "b", { x = 1, sign = sign })
        local d = E.document(q)
        local imported = E.model(d, { h = 0.18 })
        M.solve(q)
        M.solve(imported)
        for _, pt in ipairs({ { 0.3, 0.3 }, { 0.6, 0.6 } }) do
            near(
                select(3, M.sample(q, pt[1], pt[2])),
                select(3, M.sample(imported, pt[1], pt[2])),
                1e-10
            )
        end
        M.export_fem(q, "build/static/periodic-" .. sign .. ".fem")
        M.export_ans(q, "build/static/periodic-" .. sign .. ".ans")
    end
    local r = { angle = 17, magnetization = { type = "angular", order = -2, x = 0.3, y = -0.7 } }
    local law, a = X.parse_expression(X.expression(r), 1)
    near(law.center_x, 0.3, 1e-14)
    near(law.center_y, -0.7, 1e-14)
    near(law.order, -2, 1e-14)
    near(a, 17, 1e-14)
    fails(function()
        X.parse_expression('os.execute("echo no")', 1)
    end, "expressions")
    local q = model()
    M.material(q, "rotated", { mur_x = 2, mur_y = 8, material_angle = 30 })
    M.region(q, "rotated", { { 0.1, 0.1 }, { 0.9, 0.1 }, { 0.9, 0.9 }, { 0.1, 0.9 } })
    fails(function()
        E.document(q)
    end, "rotated anisotropic")
end
print("FEM circuits/periodic round trips and passive angular expressions OK")
-- Cache restoration must reconstruct paired degrees, sources and circuit data.
do
    local function cached(mode, current)
        local m = model({ cache = { mode = mode, file = "build/static/cache" } })
        S.circuit(m, "coil", { current = current })
        M.material(m, "wire", { conductivity = 5.8e7 })
        for k, x in ipairs({ 0.1, 0.7 }) do
            local r =
                M.region(m, "wire", { { x, 0.2 }, { x + 0.1, 0.2 }, { x + 0.1, 0.8 }, { x, 0.8 } })
            S.region(r, { circuit = "coil", turns = k == 1 and 10 or -10 })
        end
        M.solve(m)
        return m
    end
    local a = cached("refresh", 2)
    local b = cached("frozen", 2)
    assert(b.cache_info.status == "solution")
    near(S.result(a, "coil").linkage, S.result(b, "coil").linkage, 1e-14)
    fails(function()
        cached("frozen", 3)
    end, "frozen")
end
print("Static source cache restoration and invalidation OK")
-- A pure antiperiodic boundary removes the additive nullspace itself: cached
-- validation must not reintroduce the all-Neumann gauge used before pairing.
do
    local function problem(mode)
        local m = model({ cache = { mode = mode, file = "build/static/antiperiodic-cache" } })
        M.boundary(m, "all", { type = "neumann" })
        F.path(m, "a", { { 0, 0 }, { 0, 1 } }, false, { type = "periodic" })
        F.path(m, "b", { { 1, 0 }, { 1, 1 } }, false, { type = "periodic" })
        F.pair(m, "a", "b", { x = 1, sign = -1 })
        M.material(m, "source", { current = 2e6 })
        M.region(m, "source", { { 0.1, 0.2 }, { 0.3, 0.2 }, { 0.3, 0.8 }, { 0.1, 0.8 } })
        M.solve(m)
        return m
    end
    local a, b = problem("refresh"), problem("frozen")
    assert(b.cache_info.status == "solution")
    for i in ipairs(a.A) do
        near(a.A[i], b.A[i], 1e-14)
    end
end
-- Separate named parallel branches distribute current by resistance; every
-- branch contains its own signed outgoing/returning leg pair.
do
    local m = model({ h = 0.1 })
    M.material(m, "wire", { conductivity = 5.8e7 })
    S.circuit(m, "parallel", { type = "parallel", current = 3 })
    for branch = 1, 2 do
        for leg = 1, 2 do
            local x = 0.05 + (branch - 1) * 0.5 + (leg - 1) * 0.2
            local r =
                M.region(m, "wire", { { x, 0.2 }, { x + 0.1, 0.2 }, { x + 0.1, 0.8 }, { x, 0.8 } })
            S.region(r, {
                circuit = "parallel",
                winding = "branch" .. branch,
                turns = leg == 1 and 10 or -10,
                length = 0.2 * branch,
            })
        end
    end
    M.solve(m)
    local c = m.static_circuits.parallel
    near(c.windings.branch1.current, 2, 1e-12)
    near(c.windings.branch2.current, 1, 1e-12)
    local result = S.result(m, "parallel")
    local energy = 0
    for _, e in ipairs(m.triangles) do
        energy = energy + (e.bx * e.hx + e.by * e.hy) * e.area * 0.2 / 2
    end
    near(energy, result.linkage * 3 / 2, energy * 1e-8)
end
-- Native open boundaries refuse magnetised edges and an omitted return wire.
do
    local m = model({ xmin = -2, xmax = 2, ymin = -2, ymax = 2 })
    local outer = disk(2, 40)
    F.domain(m, { outer })
    F.path(m, "outer", outer, true, { type = "open", radius = 2 })
    M.material(m, "wire", { current = 1 })
    M.region(m, "wire", { { -0.2, -0.2 }, { 0.2, -0.2 }, { 0.2, 0.2 }, { -0.2, 0.2 } })
    fails(function()
        M.solve(m)
    end, "zero net longitudinal current")
end
print("Antiperiodic cache, multi-branch DC split and open-boundary assumptions OK")

-- FEMM sheet directions use their along-sheet scalar as the base law.
do
    local d = E.read("tests/fixtures/uniform.fem")
    d.materials[1].mu_x, d.materials[1].mu_y = 100, 2
    d.materials[1].lamtype, d.materials[1].lamfill = 1, 0.8
    local m = E.model(d)
    near(m.materials["fem:1"].mur_x, 80.2, 1e-12)
    near(m.materials["fem:1"].mur_y, 1 / (0.8 / 100 + 0.2), 1e-12)
    local q = model()
    M.material(q, "recoil", { mur_x = 3, mur_y = 3, remanence = 0.9 })
    near(q.materials.recoil.hc, 0.9 / (M.mu0 * 3), 1e-8)
end
