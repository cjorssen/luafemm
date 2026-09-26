-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local f = require("luafemm")
local p = require("luafemm-profiles")
local function near(a, b, t)
    assert(math.abs(a - b) < (t or 1e-7), string.format("%g != %g", a, b))
end
local m = f.new({
    unit = 1,
    mesher = "delaunay",
    h = 0.2,
    xmin = -1,
    xmax = 1,
    ymin = -1,
    ymax = 1,
    boundary = function(x, y)
        return x + 2 * y
    end,
})
-- B=(2,-1) everywhere, so components along an arbitrary polyline are known.
p.register(m, "broken", { { -0.5, -0.5 }, { 0.5, -0.5 }, { 0.5, 0.5 } }, { samples = 5 })
local rows = p.sample("broken")
assert(m.stats and #rows == 5)
near(p.get("broken").length, 2)
for _, r in ipairs(rows) do
    near(r.B, math.sqrt(5))
    near(r.Bx, 2)
    near(r.By, -1)
    near(r.Az, r.x + 2 * r.y)
    near(r.Hx, 2 / f.mu0, 0.02)
    near(r.Hy, -1 / f.mu0, 0.02)
end
near(rows[1].Bt, 2)
near(rows[1].Bn, -1)
near(rows[3].Bt, -1)
near(rows[3].Bn, -2) -- outgoing tangent at the corner
near(rows[5].s, 2)
near(rows[5].t, 1)
p.register(m, "reverse", { { 0.5, 0.5 }, { 0.5, -0.5 }, { -0.5, -0.5 } }, { samples = 5 })
local r = p.sample("reverse")
near(r[1].Bt, 1)
near(r[1].Bn, 2)
p.register(m, "offset", { { -0.5, 0 }, { 0.5, 0 } }, { samples = 3, offset = 0.1 })
near(p.sample("offset")[2].y, 0.1)
p.register(m, "outside", { { 0, 0 }, { 2, 0 } }, { samples = 3 })
assert(not pcall(p.sample, "outside"))
p.register(m, "outside", { { 0, 0 }, { 2, 0 } }, { samples = 3, outside = "nan" })
local text = p.coordinates("outside", "B", "s")
assert(text:find("nan", 1, true))
assert(not pcall(p.coordinates, "broken", "bad"))
assert(not pcall(p.coordinates, "broken", "B", "bad"))
assert(not pcall(p.register, m, "bad", { { 0, 0 }, { 0, 0 } }, { samples = 3 }))
assert(not pcall(p.register, m, "bad", { { 0, 0 }, { 1, 0 } }, { samples = 1 }))
local paths = require("luafemm-path")
local frame = { 1, 0, 0, 1, 0, 0 }
local bezier = "\\pgfsyssoftpath@movetotoken{0pt}{0pt}\\pgfsyssoftpath@curvetosupportatoken{0pt}{1pt}"
    .. "\\pgfsyssoftpath@curvetosupportbtoken{1pt}{1pt}\\pgfsyssoftpath@curvetotoken{1pt}{0pt}"
local poly, closed = paths.read(bezier, frame, 0.005, true)
assert(not closed and #poly > 8)
near(poly[1][1], 0)
near(poly[#poly][1], 1)
assert(not pcall(paths.read, bezier, frame, 0.005), "region must still require closure")
assert(
    not pcall(paths.read, bezier .. "\\pgfsyssoftpath@movetotoken{0pt}{0pt}", frame, 0.005, true)
)
local circle = "\\pgfsyssoftpath@movetotoken{0pt}{0pt}\\pgfsyssoftpath@linetotoken{1pt}{0pt}"
    .. "\\pgfsyssoftpath@linetotoken{1pt}{1pt}\\pgfsyssoftpath@closepathtoken{0pt}{0pt}"
poly, closed = paths.read(circle, frame, 0.01, true)
assert(closed and #poly == 4)
near(poly[1][1], poly[4][1])
near(poly[1][2], poly[4][2])
-- Registry lifetime is independent of the mutable current-model pointer.
f.current = f.new()
assert(p.get("broken").model == m)
near(p.sample("broken")[1].Bx, 2)
print(
    "Profiles: analytic components, arc length, reversal, offsets, outside policy, Bezier paths and model lifetime passed."
)
