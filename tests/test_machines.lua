-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
dofile("tests/bootstrap.lua")
local C = require("luafemm-components")
local P = require("luafemm-profiles")
local f = require("luafemm")
local mesh = require("luafemm-mesh")
local function near(a, b, t)
    assert(math.abs(a - b) < (t or 1e-8), string.format("%g != %g", a, b))
end
for _, kind in ipairs({ "rotor", "stator" }) do
    for _, shape in ipairs({ "rectangle", "trapezoid", "radial" }) do
        local c = C.make(kind, {
            shaft_radius = 8,
            housing_thickness = 2,
            slots = { count = 12, shape = shape, corner_radius = 0.1, layers = 2 },
            winding = { distribution = "sine", ampere_turns = 800 },
            conductors = { { slot = 1, layer = 2, material = "aluminum-1100" } },
        })
        local total, positive = 0, 0
        for _, r in ipairs(c.regions) do
            for _, p in ipairs(r.contours or { r.points }) do
                mesh.validate(p)
            end
            if r.turns then
                local j = C.current(r.points, 0.001, r.turns)
                local recovered = j / C.current(r.points, 1, 1) * 1e-6
                near(recovered, r.turns)
                total = total + recovered
                positive = positive + math.max(recovered, 0)
            end
        end
        near(total, 0)
        near(positive, 800)
        assert(c.anchors["slot 12 conductor 2"])
    end
end
local c = C.make("stator", {
    slots = { placement = "explicit", angles = { 20, 170 }, layers = 2 },
    overrides = { { index = 2, depth = 10, width = 6 } },
    coils = { { slots = { 1, 2 }, layer = 2, turns = 60, current = -2 } },
})
near(c.anchors["slot 1 opening"][1], 30 * math.cos(math.rad(20)))
local ni = {}
for _, r in ipairs(c.regions) do
    if r.turns then
        ni[#ni + 1] = r.turns
    end
end
near(ni[1], -120)
near(ni[2], 120)
for _, o in ipairs({
    { slots = { count = 100 } },
    { slots = { count = 2, depth = 40 } },
    { slots = { count = 2, clearance = 10 } },
    { slots = { count = 2, layers = 1.5 } },
    { slots = { count = 2, corner_radius = 10 } },
    { slots = { placement = "explicit", angles = { 0, 360 } } },
    { slots = { count = 2 }, coils = { { slots = { 1, 3 }, ampere_turns = 10 } } },
    {
        slots = { count = 2, current_density = 1 },
        coils = { { slots = { 1, 2 }, ampere_turns = 10 } },
    },
    {
        slots = { count = 2 },
        coils = { { slots = { 1, 2 }, ampere_turns = 10, turns = 2, current = 3 } },
    },
    { slots = { count = 3 }, winding = { distribution = "uniform", axis = 10 } },
    { outer_shape = "triangle" },
    { slots = { count = 2 }, overrides = { { index = 1, typo = 10 } } },
}) do
    assert(not pcall(C.make, "stator", o), "invalid machine accepted")
end
-- Independent affine field: exact radial projection and Fourier amplitude/phase.
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
P.circle(m, "analytic-polar", 0.1, -0.1, 0.5, { samples = 129, pole_pairs = 2, angle_origin = 0 })
local rows = P.sample("analytic-polar")
near(rows[1].Br, 2, 1e-6)
near(rows[1].Btheta, -1, 1e-6)
near(rows[33].Br, -1, 1e-6)
near(rows[#rows].angle, 360)
near(rows[#rows].electrical_angle, 720)
local spec = P.harmonics("analytic-polar", "Br", 15)
near(spec[1].amplitude, math.sqrt(5), 1e-6)
near(spec[1].phase, math.deg(math.atan(1, 2)), 1e-5)
near(spec.mean, 0, 1e-7)
near(spec.thd, 0, 1e-6)
assert(not pcall(P.harmonics, "analytic-polar", "Br", 64))
-- End-to-end winding comparison at equal total positive ampere-turns.
local model = require("machine-model")
local results = {}
for _, distribution in ipairs({ "uniform", "sine" }) do
    local machine = model(
        { count = 12, start_angle = 15 },
        { distribution = distribution, ampere_turns = 200 },
        4
    )
    P.circle(machine, distribution, 0, 0, 29, { samples = 241 })
    local spectrum = P.harmonics(distribution, "Br", 15)
    assert(spectrum[1].amplitude > 0.01)
    results[distribution] = spectrum
    print(
        string.format(
            "%s: B1=%.7g T, THD(2..15)=%.7g",
            distribution,
            spectrum[1].amplitude,
            spectrum.thd
        )
    )
end
assert(
    results.sine.thd < results.uniform.thd * 0.75,
    "distributed turns did not reduce harmonic content"
)
print(
    "Machine geometry, conductor currents, invalid inputs, polar projections and harmonic reduction passed."
)

-- Refine geometry and field together; compare the fundamental and entire scan.
local refined = model(
    { count = 12, start_angle = 15 },
    { distribution = "sine", ampere_turns = 200 },
    2,
    false,
    { curve_tolerance = 0.02 }
)
P.circle(refined, "sine-refined", 0, 0, 29, { samples = 241 })
local refined_spectrum = P.harmonics("sine-refined", "Br", 15)
local change = math.abs(refined_spectrum[1].amplitude / results.sine[1].amplitude - 1)
assert(change < 0.03, "fundamental did not stabilize under refinement")
assert(
    math.abs(refined_spectrum.thd - results.sine.thd) < 0.03,
    "harmonic distortion did not stabilize"
)
local coarse, fine = P.sample("sine"), P.sample("sine-refined")
local energy, error, flux = 0, 0, 0
for i = 1, #fine - 1 do
    energy = energy + fine[i].Br ^ 2
    error = error + (coarse[i].Br - fine[i].Br) ^ 2
    flux = flux + fine[i].Br
end
assert(math.sqrt(error / energy) < 0.1, "angular scan refinement discrepancy")
assert(
    math.abs(flux) / (#fine - 1) < refined_spectrum[1].amplitude * 0.01,
    "closed-circle net flux discrepancy"
)
print(
    string.format(
        "Refinement: fundamental change %.4g, relative scan RMS %.4g",
        change,
        math.sqrt(error / energy)
    )
)
-- Material topology must not depend on which concentric node is declared first.
local a = model(
    { count = 2, start_angle = 90 },
    {},
    5,
    false,
    { coils = { { slots = { 1, 2 }, ampere_turns = 100 } } }
)
local b = model(
    { count = 2, start_angle = 90 },
    {},
    5,
    true,
    { coils = { { slots = { 1, 2 }, ampere_turns = 100 } } }
)
f.solve(a)
f.solve(b)
for _, point in ipairs({ { 29, 0 }, { 0, 29 }, { 20, 20 }, { 0, 0 } }) do
    local ax, ay = f.sample(a, point[1], point[2])
    local bx, by = f.sample(b, point[1], point[2])
    near(ax, bx, 1e-6)
    near(ay, by, 1e-6)
end
print("Mesh refinement, closed-circle flux and node-order invariance passed.")

assert(not pcall(C.make, "rotor", { slots = { count = 1, neck_depth = 0.001 } }))
assert(
    not pcall(
        C.make,
        "rotor",
        { slots = { count = 1, width = 40, bottom_width = 40, neck_depth = 1, depth = 2 } }
    )
)
