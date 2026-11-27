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
local editions = {}
local function index_data(d)
    local e = { data = d, index = {}, names = {} }
    for _, r in ipairs(d.records) do
        assert(not e.index[r.id])
        e.index[r.id] = r
        e.names[r.name] = e.names[r.name] or {}
        table.insert(e.names[r.name], r.id)
    end
    editions[d.edition] = e
end
index_data(data)
--- Load only the pinned legacy library, explicitly supplied by its user.
-- This is process-local passive data; no source path is followed implicitly.
function C.configure(o)
    o.material_edition = o.material_edition or "femm5"
    o.material_file = o.material_file or ""
    assert(
        o.material_edition == "femm5" or o.material_edition == "femm4",
        "luafemm materials: unknown material edition"
    )
    if o.material_edition == "femm5" then
        assert(o.material_file == "", "luafemm materials: material file is only for femm4")
    else
        assert(o.material_file ~= "", "luafemm materials: femm4 requires an explicit material file")
        local parser = require("luafemm-catalogue")
        local records, bytes = parser.read(o.material_file)
        assert(
            parser.fingerprint(bytes) == "2addb599:148358",
            "luafemm materials: material file differs from the pinned FEMM 4 catalogue"
        )
        if not editions.femm4 then
            index_data({
                edition = "femm4",
                revision = "7d9e8ed90772ddb939a3256f976e231bcb5cca97",
                fingerprint = parser.fingerprint(bytes),
                records = records,
            })
        end
    end
end
local function edition(name)
    name = name or "femm5"
    return assert(
        editions[name],
        "luafemm materials: unavailable edition "
            .. tostring(name)
            .. "; configure an explicit material file for femm4"
    )
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
function C.get(id, name)
    local e = edition(name)
    local index, names = e.index, e.names
    local r = index[id]
    if not r and names[id] then
        assert(
            #names[id] == 1,
            "luafemm materials: ambiguous name " .. id .. "; use " .. table.concat(names[id], ", ")
        )
        r = index[names[id][1]]
    end
    assert(r, "luafemm materials: unknown material " .. tostring(id))
    local result = copy(r)
    result.edition, result.revision, result.fingerprint =
        e.data.edition, e.data.revision, e.data.fingerprint
    return result
end
--- Return independent copies of all records in upstream order.
-- @treturn table Independent catalogue records in the requested edition.
function C.list(name)
    return copy(edition(name).data.records)
end
--- Convert a catalogue record to the supported planar static material law.
-- Linear LamType=0 records also expose the audited complex-permeability data.
-- Nonlinear curves keep their static interpretation and original provenance.
-- @tparam string id Catalogue identifier.
-- @treturn table Effective properties accepted by luafemm.material.
function C.properties(id, name)
    local r = C.get(id, name)
    local p = r.fields
    local mu0 = 4 * math.pi * 1e-7
    assert(p.Mu_x == p.Mu_y or #r.bh == 0, "luafemm materials: nonlinear anisotropy not supported")
    assert(
        p.LamType == 0
            or p.LamType == 3
            or p.LamType == 4
            or p.LamType == 5
            or p.LamType == 6
            or p.LamType == 7
            or p.LamType == 8,
        "luafemm materials: unsupported catalogue wire/lamination type " .. tostring(p.LamType)
    )
    local props = {
        mur = p.LamType > 2 and 1 or p.LamFill * p.Mu_x + (1 - p.LamFill),
        coercivity = p.H_c,
        conductivity = p.Sigma * 1e6,
        current = p.J_re * 1e6,
        magnetization_angle = p.H_cAngle,
        library = r,
    }
    if p.LamType >= 3 and p.LamType <= 8 then
        props.wire_model = ({
            [3] = "round",
            [4] = "stranded",
            [5] = "litz",
            [6] = "square",
            [7] = "cca10",
            [8] = "cca15",
        })[p.LamType]
        props.wire_diameter, props.strands = p.WireD, p.NStrands
    end
    if #r.bh == 0 and p.LamType == 0 then
        props.mur = p.Mu_x
        if p.Mu_x ~= p.Mu_y then
            props.mur_x, props.mur_y = p.Mu_x, p.Mu_y
        end
        props.loss_angle_x, props.loss_angle_y = p.Phi_hx or 0, p.Phi_hy or 0
        props.library_harmonic = true
        if p.d_lam > 0 or p.LamFill < 1 then
            props.lamination = "stacked"
            props.lamination_fill = p.LamFill
            props.lamination_thickness = p.d_lam
        end
    end
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
