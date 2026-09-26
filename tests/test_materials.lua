-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local c = require("luafemm-materials")
local function near(a, b, t)
    assert(math.abs(a - b) < t, string.format("%g != %g", a, b))
end
local records = c.list()
assert(#records == 246)
local currents, laminated = 0, 0
for _, r in ipairs(records) do
    local m = f.new()
    f.material(m, "test", { library = r.id })
    local p = m.materials.test
    assert(p.library.name == r.name and #r.bh == r.fields.BHPoints)
    near(p.hc, r.fields.H_c, 1e-10)
    near(p.current, r.fields.J_re * 1e6, 1e-10)
    if p.current ~= 0 then
        currents = currents + 1
    end
    if r.fields.LamFill < 1 then
        laminated = laminated + 1
        for i, v in ipairs(r.bh) do
            near(p.bh[i][1], r.fields.LamFill * v[1] + (1 - r.fields.LamFill) * f.mu0 * v[2], 1e-14)
        end
    end
    local nu, d = f.constitutive(p, 0.37)
    assert(nu > 0 and d > 0)
    -- Each imported instance must be independent of catalogue data.
    if p.bh then
        p.bh[2][1] = 99
        assert(c.get(r.id).bh[2][1] ~= 99)
    end
    p.library.fields.Mu_x = -1
    assert(c.get(r.id).fields.Mu_x > 0)
end
assert(currents == 69 and laminated == 10)
assert(c.get("Copper").id == "copper" and c.properties("copper").mur == 1)
assert(not pcall(c.get, "Supermalloy"))
assert(#c.get("supermalloy").bh == 12)
assert(#c.get("supermalloy-metals-handbook-dc-magnetization-curves").bh == 34)
assert(not pcall(c.get, "unknown"))
local m = f.new({ h = 5 })
f.material(m, "wire", { library = "10awg-15-percent-cca" })
f.region(m, "wire", { { 0, 0 }, { 2, 0 }, { 2, 2 }, { 0, 2 } })
assert(m.regions[1].current == 1e6)
f.region(m, "wire", { { 4, 0 }, { 6, 0 }, { 6, 2 }, { 4, 2 } }, 0)
assert(m.regions[2].current == 0)
f.material(m, "custom", { library = "copper", mur = 2 })
assert(m.materials.custom.mur == 2)
print(
    "All 246 materials validated; 69 material currents, 10 fill factors, duplicate names and independent copies checked."
)
