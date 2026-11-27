-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Physical descriptors shared by analytical and finite-element providers.
-- Physical dimension is independent of the sampled/displayed dimension. The
-- descriptor is passive data, suitable for a future persistence adapter; it is
-- not an interchange format or a promise that every backend supports every law.
-- @module luafemm-physical
local P = {}
function P.copy(t)
    if type(t) ~= "table" then
        return t
    end
    local r = {}
    for k, v in pairs(t) do
        r[k] = P.copy(v)
    end
    return r
end
--- Return an owned descriptor, never mutable provider metadata.
function P.describe(m)
    assert(type(m) == "table", "luafemm physical: model required")
    if m.descriptor then
        return P.copy(m.descriptor)
    end
    assert(m.options, "luafemm physical: unknown model")
    local o = m.options
    return {
        schema = 1,
        physics = "magnetism",
        backend = o.field_model == "ideal" and "ideal circuit" or "finite elements",
        geometry = o.problem_type or "planar",
        physical_dimension = 3,
        evaluation_dimension = 2,
        length_scale = o.unit or 0.001,
        potential = "magnetic vector potential",
        reference = "boundary conditions",
        units = { length = "m", B = "T", H = "A/m" },
    }
end
--- Evaluate named SI quantities, without manufacturing a value at a singularity.
-- Analytical electric providers return nil, reason outside their valid set.
-- Magnetic providers retain their existing solve and sample semantics.
function P.sample(m, x, y, side)
    if m.descriptor and m.descriptor.physics == "electrostatics" then
        return require("luafemm-electrostatic").sample(m, x, y, side)
    end
    assert(
        m.options
            and m.options.problem_type ~= "axial"
            and m.options.frequency == 0
            and not m.phasor_model
            and m.options.field_model ~= "ideal",
        "luafemm physical: magnetic adapter requires a real planar/axisymmetric FE field"
    )
    local bx, by, a, hx, hy = require("luafemm").sample(m, x, y)
    if not bx then
        return nil, "outside magnetic domain"
    end
    return { Bx = bx, By = by, B = math.sqrt(bx * bx + by * by), Hx = hx, Hy = hy, potential = a }
end
return P
