-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
require("lfs").mkdir("build")
require("lfs").mkdir("build/interop")
local M, F, T = require("luafemm"), require("luafemm-fem"), require("luafemm-topology")
local function near(a, b, t)
    assert(math.abs(a - b) < (t or 1e-8), tostring(a) .. " != " .. tostring(b))
end
local function reject(fn, part)
    local ok, e = pcall(fn)
    assert(not ok and tostring(e):find(part, 1, true), tostring(e))
end
local d = F.read("tests/fixtures/uniform.fem")
assert(F.serialize(F.parse(F.serialize(d))) == F.serialize(d))
local m = F.model(d, { h = 3 })
assert(not m.nodes)
M.export_fem(m, "build/interop/uniform.fem")
assert(not m.nodes)
reject(function()
    M.export_ans(m, "build/interop/unsolved.ans")
end, "converged")
M.solve(m)
for _, e in ipairs(m.triangles) do
    near(e.bx, 0)
    near(e.by, -0.25)
end
M.export_ans(m, "build/interop/uniform-lua.ans")
local r = M.import_fem("build/interop/uniform.fem", { h = 2.5 })
M.solve(r)
near(select(2, M.sample(r, 8, 9)), -0.25)
m.options.depth = 13
reject(function()
    M.export_ans(m, "build/interop/stale.ans")
end, "stale")
-- Negative format/physics cases must fail before meshing or solving.
for field, v in pairs({ frequency = 50, problemtype = "axisymmetric", prevsoln = "previous.ans" }) do
    local x = F.copy(d)
    x.header[field] = v
    reject(function()
        F.model(x)
    end, "unsupported")
end
local x = F.copy(d)
x.labels[2] = F.copy(x.labels[1])
reject(function()
    F.model(x)
end, "duplicate")
x = F.copy(d)
x.materials[1].mu_y = 2
reject(function()
    F.model(x)
end, "anisotropic")
x = F.copy(d)
x.labels[1][10] = 'os.execute("false")'
reject(function()
    F.model(x)
end, "expressions")
reject(function()
    F.parse("[Format] = 4\n[NumPoints] = 2\n0 0 0 0", "broken.fem")
end, "broken.fem:")
-- All six length units describe the same physical square and affine potential.
for name, u in pairs(F.units) do
    x = F.copy(d)
    x.header.lengthunits = name
    local scale = 0.001 / u
    for _, p in ipairs(x.points) do
        p[1], p[2] = p[1] * scale, p[2] * scale
    end
    x.labels[1][1], x.labels[1][2] = 10 * scale, 10 * scale
    x.boundaries[1].a_1 = 0.25 * u
    x.header.depth = 12 * scale
    local mm = F.model(x, { h = 4 })
    M.solve(mm)
    near(select(2, M.sample(mm, 8, 9)), -0.25)
    near(mm.options.depth, 12)
end
-- Two circles made from arcs, with a real unmeshed hole.
x = F.copy(d)
x.points = { { 10, 0, 0, 0 }, { -10, 0, 0, 0 }, { 2, 0, 0, 0 }, { -2, 0, 0, 0 } }
x.segments = {}
x.arcs = {
    { 0, 1, 180, 15, 1, 0, 0 },
    { 1, 0, 180, 15, 1, 0, 0 },
    { 2, 3, 180, 15, 1, 0, 0 },
    {
        3,
        2,
        180,
        15,
        1,
        0,
        0,
    },
}
x.holes = { { 0, 0, 0 } }
x.labels = { { 5, 0, 1, 0, 0, 0, 0, 1, 0 } }
F.write_text("build/interop/annulus.fem", F.serialize(x))
local ann = F.model(
    x,
    { h = 2, curve_tolerance = 0.05, cache = { mode = "refresh", file = "build/interop/annulus" } }
)
M.solve(ann)
assert(T.characteristic(ann.topology) == 0)
for _, e in ipairs(ann.triangles) do
    near(e.bx, 0)
    near(e.by, -0.25)
end
reject(function()
    M.sample(ann, 0, 0)
end, "outside")
M.export_ans(ann, "build/interop/annulus-lua.ans")
local warm = F.model(
    x,
    { h = 2, curve_tolerance = 0.05, cache = { mode = "frozen", file = "build/interop/annulus" } }
)
M.solve(warm)
assert(warm.cache_info.status == "solution")
M.export_ans(warm, "build/interop/annulus-warm.ans")
-- Paint order, disjoint faces and nested cavities are compiled without meshing.
local native =
    M.new({ xmin = -10, xmax = 10, ymin = -10, ymax = 10, h = 2, depth = 5, mesher = "delaunay" })
M.material(native, "iron", { mur = 100 })
M.region(native, "iron", { { -7, -7 }, { 7, -7 }, { 7, 7 }, { -7, 7 } })
M.region(native, "air", { { -2, -2 }, { 2, -2 }, { 2, 2 }, { -2, 2 } })
local doc, top = F.document(native)
assert(#doc.labels == 3 and T.characteristic(top) == 1)
M.solve(native)
M.export_ans(native, "build/interop/nested-lua.ans")
-- Series excitation is the specified total NI, not NI plus the material source.
x = F.copy(d)
x.materials[1].j_re = 3
x.circuits = { { circuitname = "winding", circuittype = 1, totalamps_re = 2, totalamps_im = 0 } }
x.labels[1][5], x.labels[1][8] = 1, -30
local coil = F.model(x, { h = 3 })
M.solve(coil)
near(coil.stats.net_current, -60)
M.export_ans(coil, "build/interop/series-lua.ans")
-- FEMM cubic law reduces to the exact line and has a consistent derivative.
local B = require("luafemm-bh")
local curve = B.prepare({ { 0, 0 }, { 1, 100 }, { 2, 200 } })
for i = 0, 30 do
    local h, der = B.evaluate(curve, i / 10)
    near(h, i * 10)
    near(der, 100)
end
curve = B.prepare({ { 0, 0 }, { 0.5, 80 }, { 1, 220 }, { 1.5, 600 }, { 2, 1600 } })
for i = 1, 190 do
    local b = i / 100
    local _, der = B.evaluate(curve, b)
    local hp = B.evaluate(curve, b + 1e-6)
    local hm = B.evaluate(curve, b - 1e-6)
    near(der, (hp - hm) / 2e-6, 1e-5)
end
print("PASS: FEM/ANS, six units, arcs, exclusions, caches, circuits, topology and cubic laws")

-- Constant mixed data anchor the otherwise free additive potential.
x = F.copy(d)
x.boundaries[1] = { bdryname = "Mixed", bdrytype = 2, c0 = 100000, c1 = -200, c0i = 0, c1i = 0 }
local mixed = F.model(x, { h = 4 })
M.solve(mixed)
for _, v in ipairs(mixed.A) do
    near(v, 0.002, 1e-10)
end
M.export_fem(mixed, "build/interop/mixed.fem")
M.export_ans(mixed, "build/interop/mixed-lua.ans")
-- Two disconnected naturally bounded faces need two independent gauges.
x = F.copy(d)
x.boundaries = {}
x.points = {}
x.segments = {}
x.labels = {}
for part = 0, 1 do
    local first = #x.points
    for _, p in ipairs({ { 0, 0 }, { 4, 0 }, { 4, 4 }, { 0, 4 } }) do
        x.points[#x.points + 1] = { p[1] + part * 10, p[2], 0, 0 }
    end
    for i = 0, 3 do
        x.segments[#x.segments + 1] = { first + i, first + (i + 1) % 4, -1, 0, 0, 0 }
    end
    x.labels[#x.labels + 1] = { 2 + part * 10, 2, 1, 0, 0, 0, 0, 1, 0 }
end
local disconnected = F.model(x, { h = 1 })
M.solve(disconnected)
local gauges = 0
for _ in pairs(disconnected.gauge_nodes) do
    gauges = gauges + 1
end
assert(gauges == 2 and T.characteristic(disconnected.topology) == 2)
x.materials[1].j_re = 1
reject(function()
    M.solve(F.model(x, { h = 1 }))
end, "incompatible Neumann")
-- A union of adjacent child faces is one hole in its parent, not several.
local a = M.new({ xmin = -5, xmax = 5, ymin = -5, ymax = 5, depth = 1, h = 1, mesher = "delaunay" })
M.region(a, "air", { { -3, -2 }, { 0, -2 }, { 0, 2 }, { -3, 2 } })
M.region(a, "air", { { 0, -2 }, { 3, -2 }, { 3, 2 }, { 0, 2 } })
local adoc, atop = F.document(a)
assert(T.characteristic(atop) == 1)
local aimport = F.model(adoc, { h = 1 })
M.solve(aimport)
-- An unknown physical scalar is preserved by the codec and refused by the adapter.
x = F.copy(d)
x.materials[1].futurephysics = 1
assert(F.parse(F.serialize(x)).materials[1].futurephysics == 1)
reject(function()
    F.model(x)
end, "unsupported material field")
-- Native callbacks and nonconstant mixed data have no exact FEM representation.
local callback = M.new({
    depth = 1,
    source = function()
        return 0
    end,
})
reject(function()
    F.document(callback)
end, "callbacks")
local varying = M.new({ depth = 1 })
M.boundary(varying, "left", { type = "neumann", dy = 1 })
reject(function()
    F.document(varying)
end, "varying mixed")
print("PASS: mixed boundaries, disconnected gauges, shared child cycles and strict export")

local duplicate_bh = F.serialize(d):gsub("<BHPoints> = 0", "<BHPoints> = 0\n<BHPoints> = 0")
reject(function()
    F.parse(duplicate_bh)
end, "duplicate property bhpoints")
