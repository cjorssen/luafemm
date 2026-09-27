-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- luafemm-materials implementation.
-- @module luafemm-materials
-- Complete FEMM material catalogue, with a planar static interpretation.
local data = require("luafemm-materials-data")
local C = { revision = data.revision }
local index, names = {}, {}
for _, r in ipairs(data.records) do
    assert(not index[r.id])
    index[r.id] = r
    names[r.name] = names[r.name] or {}
    table.insert(names[r.name], r.id)
end
local function copy(t)
    local out = {}
    for k, v in pairs(t) do
        out[k] = type(v) == "table" and copy(v) or v
    end
    return out
end
--- Return an independent catalogue record by identifier or unambiguous FEMM name.
-- @tparam string id Identifier, or exact original name.
-- @treturn table Record with id, name, path, fields and original bh data.
-- @raise Unknown identifiers and ambiguous names are errors.
function C.get(id)
    local r = index[id]
    if not r and names[id] then
        assert(
            #names[id] == 1,
            "luafemm materials: ambiguous name " .. id .. "; use " .. table.concat(names[id], ", ")
        )
        r = index[names[id][1]]
    end
    assert(r, "luafemm materials: unknown material " .. tostring(id))
    return copy(r)
end
--- Return independent copies of all records in upstream order.
-- @treturn table Array of 246 catalogue records.
function C.list()
    return copy(data.records)
end
--- Convert a catalogue record to the supported planar static material law.
-- Conductivity and harmonic/loss parameters remain in library.fields.
-- @tparam string id Catalogue identifier.
-- @treturn table Effective properties accepted by luafemm.material.
function C.properties(id)
    local r = C.get(id)
    local p = r.fields
    local mu0 = 4 * math.pi * 1e-7
    assert(p.Mu_x == p.Mu_y, "luafemm materials: anisotropy not supported")
    assert(
        p.LamType == 0 or p.LamType == 3 or p.LamType == 8,
        "luafemm materials: unsupported lamination type"
    )
    local props = {
        mur = p.LamType > 2 and 1 or p.LamFill * p.Mu_x + (1 - p.LamFill),
        coercivity = p.H_c,
        current = p.J_re * 1e6,
        magnetization_angle = p.H_cAngle,
        library = r,
    }
    if #r.bh > 0 then
        assert(
            p.LamType == 0 and (p.LamFill == 1 or p.H_c == 0),
            "luafemm materials: unsupported nonlinear mixture"
        )
        props.bh = {}
        for i, v in ipairs(r.bh) do
            props.bh[i] = { p.LamFill * v[1] + (1 - p.LamFill) * mu0 * v[2], v[2] }
        end
    end
    return props
end
return C
