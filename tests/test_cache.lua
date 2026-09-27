-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Exercise disk round trips, layered invalidation and fail-closed frozen reads.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local profiles = require("luafemm-profiles")
local lfs = require("lfs")
lfs.mkdir("build")
lfs.mkdir("build/tests")
local stem = "build/tests/cache-test"
local file = stem .. ".lfc"
local function bytes()
    local stream = assert(io.open(file, "rb"))
    local result = stream:read("*a")
    stream:close()
    return result
end
local function write(text)
    local stream = assert(io.open(file, "wb"))
    assert(stream:write(text))
    assert(stream:close())
end
local function fails(fn, pattern)
    local ok, message = pcall(fn)
    assert(not ok and tostring(message):find(pattern, 1, true), tostring(message))
end
local function model(mode, changes)
    local c = changes or {}
    local options = {
        xmin = -6,
        xmax = 6,
        ymin = -6,
        ymax = 6,
        h = c.h or 2,
        mesher = c.mesher or "delaunay",
        tolerance = c.tolerance or 1e-7,
        max_newton = c.max_newton or 60,
        unit = c.unit or 0.001,
        mesh_tolerance = c.mesh_tolerance,
        cache = { file = stem, mode = mode },
    }
    local m = f.new(options)
    m.options.log = function(message)
        if message:match("^Newton") then
            m.iteration_logs = (m.iteration_logs or 0) + 1
        end
    end
    f.material(m, "iron", { bh = { { 0, 0 }, { 1, c.bh or 100 }, { 2, 20000 } } })
    f.material(m, "magnet", { mur = 4, coercivity = c.hc or 45000 })
    f.region(m, "iron", { { -3, -4 }, { 3, -4 }, { 3, -2 }, { -3, -2 } })
    f.region(
        m,
        "magnet",
        { { -2, -1 }, { c.x or 2, -1 }, { c.x or 2, 2 }, { -2, 2 } },
        c.current or 0,
        c.angle or 90,
        c.local_h or 0
    )
    return m
end
local function same(a, b)
    assert(#a.A == #b.A and #a.triangles == #b.triangles)
    for i, value in ipairs(a.A) do
        assert(value == b.A[i], "potential did not round-trip exactly")
    end
    local aa, bb = { f.sample(a, 0.31, 3.17) }, { f.sample(b, 0.31, 3.17) }
    for i = 1, #aa do
        assert(aa[i] == bb[i], "derived field did not round-trip exactly")
    end
    profiles.register(a, "cold", { { -4, 3 }, { 0, 4 }, { 4, 3 } }, { samples = 19 })
    profiles.register(b, "warm", { { -4, 3 }, { 0, 4 }, { 4, 3 } }, { samples = 19 })
    for _, component in ipairs({ "B", "Bx", "By", "H", "Hx", "Hy", "Az", "Bt", "Bn" }) do
        assert(profiles.coordinates("cold", component) == profiles.coordinates("warm", component))
    end
    local ac, bc = f.contours(a, 9), f.contours(b, 9)
    assert(#ac == #bc)
    for i, path in ipairs(ac) do
        assert(#path == #bc[i] and path.level == bc[i].level)
        for j, point in ipairs(path) do
            assert(point.x == bc[i][j].x and point.y == bc[i][j].y)
        end
    end
end

for _, mesher in ipairs({ "grid", "delaunay" }) do
    os.remove(file)
    local cold = model("auto", { mesher = mesher })
    f.solve(cold)
    assert(cold.cache_info.status == "miss" and cold.iteration_logs > 0)
    local saved = bytes()
    local warm = model("auto", { mesher = mesher })
    warm.frame = { 2, 0, 0, 3, 17, -23 }
    f.material(warm, "unused", { mur = 789 })
    local mesh = f.mesh
    f.mesh = function()
        error("unexpected mesh generation")
    end
    f.solve(warm)
    f.mesh = mesh
    assert(warm.cache_info.status == "solution" and not warm.iteration_logs)
    assert(warm.frame[5] == 17 and bytes() == saved)
    same(cold, warm)
    -- Repeated solve is idempotent; profiles can outlive their source picture.
    f.solve(warm)
    assert(not warm.iteration_logs and bytes() == saved)
    local frozen = model("frozen", { mesher = mesher })
    f.mesh(frozen)
    assert(not frozen.stats and frozen.cache_info.status == "mesh")
    f.solve(frozen)
    assert(frozen.cache_info.status == "solution" and bytes() == saved)
    same(cold, frozen)
    if mesher == "grid" then
        assert(not frozen.mesh_stats and #frozen.xs == #cold.xs)
    else
        assert(frozen.mesh_stats.min_angle == cold.mesh_stats.min_angle)
    end
    local refresh = model("refresh", { mesher = mesher })
    f.solve(refresh)
    assert(refresh.cache_info.status == "refresh" and refresh.iteration_logs > 0)
    same(cold, refresh)
    print(mesher .. ": cold, warm, frozen, mesh-first and refresh agree exactly")
end

local baseline = bytes()
for name, change in pairs({
    current = { current = 2e6 },
    bh = { bh = 200 },
    coercivity = { hc = 30000 },
    angle = { angle = 30 },
    tolerance = { tolerance = 1e-9 },
    iterations = { max_newton = 80 },
}) do
    write(baseline)
    local stale = model("frozen", change)
    fails(function()
        f.solve(stale)
    end, "frozen cache unavailable")
    assert(not stale.nodes and bytes() == baseline)
    local changed = model("auto", change)
    local mesh = f.mesh
    f.mesh = function()
        error("physical change regenerated the mesh")
    end
    f.solve(changed)
    f.mesh = mesh
    assert(changed.cache_info.status == "mesh" and changed.iteration_logs > 0)
    local reference = model("off", change)
    f.solve(reference)
    same(changed, reference)
    local hot = model("frozen", change)
    f.solve(hot)
    same(changed, hot)
    print(name .. ": reused mesh and recomputed solution")
end
for name, change in pairs({
    geometry = { x = 2.5 },
    spacing = { h = 1.5 },
    local_spacing = { local_h = 1 },
    unit = { unit = 0.002 },
    mesh_tolerance = { mesh_tolerance = 1e-8 },
    mesher = { mesher = "grid" },
}) do
    write(baseline)
    fails(function()
        f.solve(model("frozen", change))
    end, "frozen cache unavailable")
    local changed = model("auto", change)
    f.solve(changed)
    assert(changed.cache_info.status == "miss" and changed.iteration_logs > 0)
    print(name .. ": invalidated mesh and solution")
end

os.remove(file)
fails(function()
    f.solve(model("frozen"))
end, "frozen cache unavailable")
assert(not io.open(file, "rb"))
local mesh_only = model("auto")
f.mesh(mesh_only)
local mesh_bytes = bytes()
local frozen_mesh = model("frozen")
f.mesh(frozen_mesh)
fails(function()
    f.solve(frozen_mesh)
end, "solution missing")
assert(not frozen_mesh.stats and bytes() == mesh_bytes)
local finish = model("auto")
f.solve(finish)
assert(finish.cache_info.status == "mesh")

-- A disabled cache neither reads nor replaces an existing record.
write("This is deliberately not a cache")
local disabled = model("off")
f.solve(disabled)
assert(disabled.cache_info.status == "off" and bytes() == "This is deliberately not a cache")
for _, bad in ipairs({
    "",
    baseline:sub(1, -2),
    baseline .. "extra",
    "return os.execute('echo never execute cache contents')",
}) do
    write(bad)
    fails(function()
        f.solve(model("frozen"))
    end, "frozen cache unavailable")
    assert(bytes() == bad)
    local recover = model("auto")
    f.solve(recover)
    assert(recover.cache_info.status == "miss")
end

-- Recompute the transport checksum after deliberately changing metadata;
-- compatibility must be checked independently of file integrity.
local function with_payload(payload)
    local a, b = 1, 0
    for i = 1, #payload do
        a = (a + payload:byte(i)) % 65521
        b = (b + a) % 65521
    end
    return string.format("LUAFEMM-CACHE %d %08x\n", #payload, b * 65536 + a) .. payload
end
local payload = baseline:match("^[^\n]+\n(.*)$")
write((baseline:gsub("0%.6%.0%-dev", "0.7.0-dev")))
fails(function()
    f.solve(model("frozen"))
end, "checksum mismatch")
write(with_payload((payload:gsub("0%.6%.0%-dev", "0.7.0-dev"))))
fails(function()
    f.solve(model("frozen"))
end, "incompatible package version")
local recover = model("auto")
f.solve(recover)
assert(recover.cache_info.status == "miss")
local invalid_node, replacements =
    payload:gsub("s5:nodesm(%d+):n1;m2:n1;n[^;]+;", "s5:nodesm%1:n1;m2:n1;n1e99;", 1)
assert(replacements == 1)
write(with_payload(invalid_node))
fails(function()
    f.solve(model("frozen"))
end, "node outside domain")
local repaired = model("auto")
f.solve(repaired)
assert(repaired.cache_info.status == "miss")

-- A failed nonlinear solve must never be advertised as a saved solution.
local difficult = model("refresh", { current = 2e10, max_newton = 1 })
fails(function()
    f.solve(difficult)
end, "Newton failed to converge")
local partial = model("frozen", { current = 2e10, max_newton = 1 })
f.mesh(partial)
assert(partial.cache_info.status == "mesh" and not partial.stats)
fails(function()
    f.solve(partial)
end, "solution missing")

-- An interrupted replacement must preserve the previous readable result.
write(baseline)
local rename = os.rename
os.rename = function()
    return nil, "simulated rename failure"
end
fails(function()
    f.solve(model("refresh"))
end, "cannot replace")
os.rename = rename
assert(bytes() == baseline)
for name in lfs.dir("build/tests") do
    assert(not name:match("^cache%-test%.lfc%.tmp%-"), "temporary cache leaked")
end

for _, key in ipairs({ "boundary", "source" }) do
    fails(function()
        f.new({
            [key] = function()
                return 0
            end,
            cache = { file = stem },
        })
    end, "callbacks require cache mode off")
end
fails(function()
    f.new({ cache = { mode = "unknown" } })
end, "unknown mode")
fails(function()
    f.new({ cache = {} })
end, "filename stem")
local unwritable = model("auto")
unwritable.options.cache.file = "build/tests/absent-cache-directory/cache"
unwritable.cache_info.file = unwritable.options.cache.file .. ".lfc"
fails(function()
    f.solve(unwritable)
end, "cannot write")
os.remove(file)
print("All persistent cache tests passed.")
