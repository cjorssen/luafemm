-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Native previous-solution export with an explicit companion DC mesh.
-- FEMM reuses element labels and boundary-property indices from this mesh.
-- Both documents must therefore have identical geometry/label ordering.
local X = {}
local function check(v, s)
    assert(v, "luafemm response exchange: " .. s)
end
function X.validate(m)
    require("luafemm-operating").validate(m)
    check(m.options.problem_type == "planar", "response exchange is planar only")
    check(
        not require("luafemm-airgap").active(m),
        "analytical gaps have no native P1 answer representation"
    )
    for _, p in pairs(m.materials) do
        check(not p.nonlinear_sheet, "nonlinear sheet response exchange is unsupported")
        check(
            not p.bh or p.interpolation == "femm",
            "nonlinear response exchange requires interpolation=femm"
        )
    end
    for _, r in ipairs(m.regions) do
        check(
            not r.magnetization or r.magnetization.type == "uniform",
            "response exchange requires uniform magnetisation in each region"
        )
    end
end
--- Write a companion first and publish the response only after it succeeds.
-- A leaf filename keeps the pair relocatable; native solvers must run in that
-- directory. No imported header filename is opened automatically.
function X.write(m, path, format, biasfile)
    -- Validate the destination before creating even the companion file.
    check(
        type(path) == "string"
            and path:match("%S")
            and not path:find("%c")
            and path:match("[^/\\]+$"),
        "invalid output filename"
    )
    check(biasfile == nil or type(biasfile) == "string", "invalid bias filename")
    X.validate(m)
    local F = require("luafemm-fem")
    F.check_companion(path)
    local d = F.document(m)
    local bias = m.exchange_bias
    check(bias, "this response has no owned DC problem; export from a captured operating point")
    check(
        bias.solution_signature == require("luafemm-cache").signature(bias),
        "changed companion bias"
    )
    local bd = F.document(bias)
    check(#d.labels == #bd.labels and #d.boundaries == #bd.boundaries, "companion topology differs")
    for i, l in ipairs(d.labels) do
        check(
            l[1] == bd.labels[i][1] and l[2] == bd.labels[i][2],
            "companion label ordering differs"
        )
    end
    local dir, leaf = path:match("^(.*[/\\])([^/\\]+)$")
    dir, leaf = dir or "", leaf or path
    biasfile = biasfile or ""
    if biasfile == "" then
        biasfile = leaf:gsub("%.[^.]+$", "") .. "-bias.ans"
    end
    check(
        biasfile:match("%.ans$") and not biasfile:find("[/\\%c]") and biasfile ~= leaf,
        "bias file must be a distinct .ans leaf filename"
    )
    F.check_companion(dir .. biasfile)
    require("luafemm-ans").write(bias, dir .. biasfile, { mesh_document = d })
    if format == "fem" then
        d.header.prevsoln = biasfile
        F.write_text(path, F.serialize(d), m.import_document and m.import_document.filename)
    else
        d = require("luafemm-ans").write(m, path, { previous = biasfile })
    end
    return d
end
--- Import an editable response problem only with an explicitly supplied bias.
-- Match native record ordering before reusing any triangle/label indices. This
-- intentionally rejects a merely similar geometry: no bias interpolation occurs.
function X.import(document, options, biasfile)
    local F = require("luafemm-fem")
    local A = require("luafemm-ans")
    local a = A.read(biasfile)
    check(
        (a.previous or 0) == 0 and (a.document.header.frequency or 0) == 0,
        "a DC non-response bias is required"
    )
    check(
        document.header.problemtype == "planar" and a.document.header.problemtype == "planar",
        "response exchange is planar only"
    )
    check(
        document.header.lengthunits == a.document.header.lengthunits,
        "companion length units differ"
    )
    local function equal(a, b)
        if type(a) ~= type(b) then
            return false
        end
        if type(a) ~= "table" then
            return a == b
        end
        for k, v in pairs(a) do
            if not equal(v, b[k]) then
                return false
            end
        end
        for k in pairs(b) do
            if a[k] == nil then
                return false
            end
        end
        return true
    end
    for _, key in ipairs({ "points", "segments", "arcs", "holes", "labels" }) do
        check(
            equal(document[key], a.document[key]),
            "companion " .. key .. " differ; exact native ordering is required"
        )
    end
    check(
        document.header.prevtype == 1 or document.header.prevtype == 2,
        "expected an incremental or frozen FEM problem"
    )
    local bias = A.model(a)
    local m = F.model(document, options, true)
    local nodes, elements = {}, {}
    for i, p in ipairs(bias.nodes) do
        nodes[i] = { p[1], p[2] }
    end
    for i, e in ipairs(bias.triangles) do
        elements[i] = F.copy(e.ids)
    end
    m.operating = {
        mode = document.header.prevtype == 1 and "incremental" or "frozen",
        A = F.copy(bias.A),
        nodes = F.copy(nodes),
        elements = F.copy(elements),
        currents = {},
        provenance = bias.solution_signature,
    }
    for i, e in ipairs(bias.triangles) do
        m.operating.currents[i] = e.current
    end
    m.exchange_bias = require("luafemm-operating").capture(bias)
    require("luafemm").prepare_mesh(m, nodes, elements, bias.exterior_nodes)
    m.mesh_stats = F.copy(bias.mesh_stats)
    m.operating_contract = require("luafemm-cache").signature(m)
    X.validate(m)
    return m
end
return X
