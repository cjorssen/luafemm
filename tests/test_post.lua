-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f, post = require("luafemm"), require("luafemm-post")
local profiles = require("luafemm-profiles")
local function near(a, b, t)
    assert(math.abs(a - b) < t, string.format("%.12g != %.12g (tol %.3g)", a, b, t))
end
local function fails(fn, pattern)
    local ok, err = pcall(fn)
    assert(not ok and tostring(err):find(pattern, 1, true), tostring(err))
end
local function rectangle(x0, y0, x1, y1)
    return { { x0, y0 }, { x1, y0 }, { x1, y1 }, { x0, y1 } }
end
-- Exact affine patch: independently known volume, average fields and energy.
local m = f.new({
    unit = 1,
    depth = 0.3,
    h = 0.2,
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1,
    boundary = function(x, y)
        return x + 2 * y
    end,
})
post.select(m, "all")
local r = post.measure(m, { name = "patch", selection = "all" })
near(r.area, 1, 1e-12)
near(r.volume, 0.3, 1e-12)
near(r.Bx, 2, 1e-10)
near(r.By, -1, 1e-10)
near(r.energy, 5 * 0.3 / (2 * f.mu0), 1e-5)
near(r.energy, r.coenergy, 1e-6)
post.select(
    m,
    "ring",
    { contours = { rectangle(0, 0, 1, 1), rectangle(0.2, 0.2, 0.8, 0.8) }, fill_rule = "even odd" }
)
near(post.measure(m, { name = "ring-result", selection = "ring" }).area, 0.64, 1e-12)
post.select(m, "empty", { contours = { rectangle(2, 2, 3, 3) } })
fails(function()
    post.measure(m, { name = "bad", selection = "empty" })
end, "no elements")
profiles.register(
    m,
    "box",
    { { 0.1, 0.1 }, { 0.9, 0.1 }, { 0.9, 0.9 }, { 0.1, 0.9 }, { 0.1, 0.1 } }
)
r = post.measure(m, { name = "zero", quantity = "force", method = "maxwell", profile = "box" })
near(r.force, 0, 1e-3)
near(r.torque, 0, 1e-3)
fails(function()
    post.value("patch", "Fx")
end, "unavailable")
fails(function()
    post.measure(m, { name = "bad", quantity = "force", selection = "all" })
end, "touches")

-- Exact primitive checks at knots, between knots, and above the table.
local p = { bh = { { 0, 0 }, { 1, 100 }, { 2, 400 } }, hc = 0 }
near(post.primitive(p, 1.5), 137.5, 1e-10)
near(post.primitive(p, 3), 300 + 400 + 1 / (2 * f.mu0), 1e-8)
f.material(f.new(), "unused", {})
local cm = f.new({ interpolation = "femm" })
f.material(cm, "curve", p)
for _, b in ipairs({ 0.3, 1.3, 2.5 }) do
    local eps = 1e-5
    local derivative = (
        post.primitive(cm.materials.curve, b + eps) - post.primitive(cm.materials.curve, b - eps)
    ) / (2 * eps)
    local nu = f.constitutive(cm.materials.curve, b)
    near(derivative, nu * b / f.mu0, 1e-5)
end

local function body(kind, h)
    local a = f.new({
        unit = 1,
        depth = 0.3,
        h = h,
        xmin = -1,
        xmax = 1,
        ymin = -1,
        ymax = 1,
        boundary = function(x)
            return -0.2 * x
        end,
    })
    f.material(a, "body", kind == "magnet" and { mur = 1, coercivity = 1000 } or { mur = 1 })
    f.region(a, "body", rectangle(-0.2, -0.1, 0.2, 0.1), kind == "wire" and 2000 or 0)
    post.select(a, "body", { material = "body" })
    profiles.register(
        a,
        "around",
        { { -0.45, -0.35 }, { 0.45, -0.35 }, { 0.45, 0.35 }, { -0.45, 0.35 }, { -0.45, -0.35 } }
    )
    return a
end
-- J x B and m x B provide nonzero independent force and torque benchmarks.
for _, kind in ipairs({ "wire", "magnet" }) do
    local a = body(kind, 0.05)
    local expected = kind == "wire" and -2000 * 0.08 * 0.3 * 0.2 or 1000 * 0.08 * 0.3 * 0.2
    for _, method in ipairs({ "weighted", "maxwell", "lorentz" }) do
        if kind == "wire" or method ~= "lorentz" then
            local q = post.measure(a, {
                name = "load",
                selection = "body",
                quantity = "force",
                method = method,
                profile = "around",
            })
            if kind == "wire" then
                near(q.Fx, expected, 0.04 * math.abs(expected))
                near(q.Fy, 0, 0.01)
            else
                near(q.torque, expected, 0.04 * expected)
                near(q.force, 0, 0.3)
            end
        end
    end
    if kind == "magnet" then
        fails(function()
            post.measure(
                a,
                { name = "bad", selection = "body", quantity = "force", method = "lorentz" }
            )
        end, "nonmagnetic")
    end
end
-- Reversing a stress contour must not reverse the outward mechanical force.
local a = body("wire", 0.05)
local forward = post.measure(
    a,
    { name = "forward", quantity = "force", method = "maxwell", profile = "around" }
)
local reversed = {}
for i = #profiles.get("around").points, 1, -1 do
    reversed[#reversed + 1] = profiles.get("around").points[i]
end
profiles.register(a, "reverse", reversed)
local backward = post.measure(a, {
    name = "backward",
    quantity = "force",
    method = "maxwell",
    profile = "reverse",
    origin_y = 0.7,
})
near(forward.Fx, backward.Fx, 1e-10)
near(backward.torque, forward.torque + 0.7 * forward.Fx, 1e-10)
profiles.register(
    a,
    "inside",
    { { -0.1, -0.05 }, { 0.1, -0.05 }, { 0.1, 0.05 }, { -0.1, 0.05 }, { -0.1, -0.05 } }
)
fails(function()
    post.measure(a, { name = "bad", quantity = "force", method = "maxwell", profile = "inside" })
end, "source-free air")
profiles.register(a, "open", { { -0.4, -0.4 }, { 0.4, 0.4 } })
fails(function()
    post.measure(a, { name = "bad", quantity = "force", method = "maxwell", profile = "open" })
end, "closed")
fails(function()
    post.measure(f.new(), { name = "bad" })
end, "depth")
fails(function()
    post.measure(m, { name = "bad", selection = "body" })
end, "another problem")
assert(post.coordinates("forward,backward", "Fx", "parameter"):find("(", 1, true))
-- Virtual work at fixed currents with homogeneous Dirichlet exterior data.
-- This differentiates total coenergy of independently meshed displaced models.
local function coils(cx)
    local model =
        f.new({ unit = 1, depth = 0.3, h = 0.05, xmin = -1, xmax = 1, ymin = -1, ymax = 1 })
    f.material(model, "left", { mur = 1 })
    f.material(model, "right", { mur = 1 })
    f.region(model, "left", rectangle(-0.5, -0.1, -0.3, 0.1), 1000)
    f.region(model, "right", rectangle(cx - 0.1, -0.1, cx + 0.1, 0.1), -1000)
    post.select(model, "moving", { material = "right" })
    local energy = post.measure(model, { name = "system" }).coenergy
    local force =
        post.measure(model, { name = "moving-force", selection = "moving", quantity = "force" }).Fx
    return energy, force
end
local em = coils(0.35)
local ep = coils(0.45)
local _, force = coils(0.4)
local derivative = (ep - em) / 0.1
near(force, derivative, math.abs(force) * 0.08)
-- The declared magnet energy reference must include -Hc m dot B.
local magnet = f.new({
    unit = 1,
    depth = 0.3,
    h = 0.2,
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1,
    boundary = function(_, y)
        return 0.1 * y
    end,
})
f.material(magnet, "magnet", { mur = 1, coercivity = 100000 })
f.region(magnet, "magnet", rectangle(0, 0, 1, 1))
local energy = post.measure(magnet, { name = "magnet-energy" })
near(energy.energy, 0.3 * (0.1 ^ 2 / (2 * f.mu0) - 100000 * 0.1), 1e-5)
near(energy.coenergy, 0.3 * 0.1 ^ 2 / (2 * f.mu0), 1e-5)
assert(energy.energy < 0)
magnet.options.depth = 0.6
fails(function()
    post.measure(magnet, { name = "stale" })
end, "stale solution")
print(
    "PASS: integrals, nonlinear energy, dipole torque, conductor force, contour orientation and diagnostics"
)
