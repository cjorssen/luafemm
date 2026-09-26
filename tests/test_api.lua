-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Ownership and the escaped TeX boundary are part of the public RC contract.
dofile("tests/bootstrap.lua")
local femm = require("luafemm")
local bridge = require("luafemm-tex")

local options = { h = 1, xmin = -3, xmax = 3, ymin = -3, ymax = 3 }
local model = femm.new(options)
assert(options.unit == nil and options.mesher == nil, "constructor mutated caller options")
options.h = 100
assert(model.options.h == 1)
local material = { bh = { { 0, 0 }, { 1, 100 } } }
femm.material(model, "iron", material)
material.bh[2][2] = 999
assert(model.materials.iron.bh[2][2] == 100 and material.hc == nil)
local points = { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } }
femm.region(model, "iron", points)
points[1][1] = -50
assert(model.regions[1].points[1][1] == -1)
femm.solve(model)
assert(model.stats.peak == 0)
assert(not pcall(femm.sample, model, 0 / 0, 0))
assert(not pcall(femm.sample, model, 0, math.huge))
assert(not pcall(femm.sample, femm.new(), 0, 0))

assert(bridge.number("1e-3", "test") == 0.001)
assert(bridge.number("", "test", 17) == 17)
for _, value in ipairs({ "", "1+2", "nan", "inf", "1e999", "0); error('injection') --" }) do
    local ok, message = pcall(bridge.number, value, "numeric-test")
    assert(not ok and message:find("expected a finite number for numeric%-test"))
end
-- Grid meshing has no Delaunay diagnostics. Never substitute solver seconds
-- when a caller requests an unavailable mesh statistic with the same name.
femm.current = model
local ok, message = pcall(bridge.statistic, "seconds", true)
assert(not ok and message:find("statistic unavailable: seconds", 1, true))
femm.current = nil
assert(not pcall(femm.sample, nil, 0, 0))
assert(rawget(_G, "luafemm") == nil, "requiring Lua modules must not create globals")
print("Input ownership, field validation and escaped numeric boundary passed.")
