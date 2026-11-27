-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- FEMM's axisymmetric Kelvin exterior coefficient, evaluated per element.
-- Coordinates of the exterior are computational, not physical-space samples.
-- A mapped source-free isotropic linear region has nu_eff=nu*rho^2*Ri/Ro^3.
-- The source region and paired/interface geometry remain the caller's model.
local K = {}
local function check(v, s)
    assert(v, "luafemm exterior: " .. s)
end
function K.prepare(m, triangles, nodes)
    local o = m.options
    for _, e in ipairs(triangles) do
        if e.kelvin then
            check(
                o.problem_type == "axisymmetric" and not m.operating and not m.transient,
                "Kelvin regions require an ordinary axisymmetric model"
            )
            local ri, ro, zo = o.exterior_inner or 0, o.exterior_outer or 0, o.exterior_center or 0
            check(
                ri > 0 and ro > 0 and zo == zo and math.abs(zo) < math.huge,
                "exterior physical and mapped radii must be positive"
            )
            check(
                not e.p.bh
                    and e.p.hc == 0
                    and e.p.conductivity == 0
                    and e.current == 0
                    and e.p.mur_x == e.p.mur_y
                    and not require("luafemm-losses").active(e.p),
                "Kelvin region must be linear isotropic, source-free and nonconducting"
            )
            local r, z = 0, 0
            for _, id in ipairs(e.ids) do
                r = r + nodes[id][1] / 3
                z = z + nodes[id][2] / 3
            end
            ri, ro, zo = ri * o.unit, ro * o.unit, zo * o.unit
            e.kelvin_factor = (r * r + (z - zo) ^ 2) * ri / ro ^ 3
            check(
                e.kelvin_factor > 0 and e.kelvin_factor < math.huge,
                "singular exterior coefficient"
            )
        end
    end
end
return K
