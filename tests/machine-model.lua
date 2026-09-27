-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Shared numerical fixture; no TeX state and no native dependency.
local f = require("luafemm")
local C = require("luafemm-components")
return function(slots, winding, h, rotor_first, extra)
    local m = f.new({
        mesher = "delaunay",
        xmin = -52,
        xmax = 52,
        ymin = -52,
        ymax = 52,
        h = h or 4,
        depth = 10,
    })
    f.material(m, "iron", { mur = 2000 })
    f.material(m, "air", { mur = 1 })
    f.material(m, "copper", { mur = 1 })
    local o = { core_material = "iron", slots = slots, winding = winding, curve_tolerance = 0.04 }
    for k, v in pairs(extra or {}) do
        o[k] = v
    end
    local stator = C.make("stator", o)
    local rotor = C.make("rotor", { radius = 28, core_material = "iron", curve_tolerance = 0.04 })
    local parts = rotor_first and { rotor, stator } or { stator, rotor }
    for _, c in ipairs(parts) do
        for _, r in ipairs(c.regions) do
            local j = r.turns and C.current(r.points, m.options.unit, r.turns)
                or r.current_density
                or 0
            f.region_contours(
                m,
                r.material,
                r.contours or { r.points },
                j,
                0,
                r.mesh_size,
                "even odd"
            )
        end
    end
    -- A narrow air-only refinement band, wholly inside the polygonal gap.
    local rings = {}
    for _, radius in ipairs({ 28.2, 29.8 }) do
        local p = {}
        for i = 1, 120 do
            local a = 2 * math.pi * i / 120
            p[i] = { radius * math.cos(a), radius * math.sin(a) }
        end
        rings[#rings + 1] = p
    end
    f.region_contours(m, "air", rings, 0, 0, (h or 4) / 4, "even odd")
    return m, stator
end
