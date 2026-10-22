-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Static, non-incremental magnetic answer-file output from a solved snapshot.
-- Coordinates use file units, A uses T m, element label references are zero-based.
-- @module luafemm-ans
local A = {}
--- Write the existing mesh and potentials, without solving or remeshing.
-- Linear material post-processing agrees exactly in the represented model;
-- nonlinear H agreement requires interpolation='femm' at material preparation.
-- @tparam table m Converged model, including one restored from a valid cache.
-- @tparam string path Destination including .ans.
-- @treturn table Embedded computational problem document.
function A.write(m, path)
    assert(m and m.stats, "luafemm ANS: a converged solution is required")
    assert(
        m.options.field_model ~= "ideal",
        "luafemm ANS: ideal floating boundaries cannot be exported"
    )
    assert(
        m.solution_signature == require("luafemm-cache").signature(m),
        "luafemm ANS: stale solution; rebuild and solve the changed model"
    )
    for _, e in ipairs(m.triangles) do
        assert(
            not e.p.bh or e.p.interpolation == "femm",
            "luafemm ANS: nonlinear export requires interpolation=femm before solving"
        )
    end
    local F = require("luafemm-fem")
    local d, t, labels = F.document(m)
    local element_labels, label_current = {}, {}
    local spatial = false
    for _, r in ipairs(m.regions) do
        spatial = spatial
            or (
                m.materials[r.material].hc > 0
                and r.magnetization
                and r.magnetization.type ~= "uniform"
            )
    end
    if spatial then
        -- FEMM evaluates angle laws at element centroids; luafemm integrates
        -- their source with three-point quadrature. Preserve the solved weak
        -- source exactly by embedding one constant-source block per triangle.
        local old = d.materials
        d.points, d.segments, d.labels, d.materials, d.circuits, d.holes = {}, {}, {}, {}, {}, {}
        for i, p in ipairs(m.nodes) do
            d.points[i] = { p[1] / m.options.unit, p[2] / m.options.unit, 0, 0 }
        end
        local seen = {}
        for i, e in ipairs(m.triangles) do
            local x, y = 0, 0
            for _, id in ipairs(e.ids) do
                x = x + m.nodes[id][1] / (3 * m.options.unit)
                y = y + m.nodes[id][2] / (3 * m.options.unit)
            end
            local face = require("luafemm-topology").locate(t, x, y)
            local material = F.copy(old[labels[face.id] + 1])
            local magnitude = math.sqrt(e.mx ^ 2 + e.my ^ 2)
            material.h_c = e.p.hc * magnitude
            material.blockname = "snapshot element " .. i
            material.j_re = e.current / 1e6
            d.materials[i] = material
            d.labels[i] = { x, y, i, 0, 0, math.deg(math.atan(e.my, e.mx)), 0, 1, 0 }
            element_labels[i] = i - 1
            label_current[i] = 0
            for j = 1, 3 do
                local a, b = e.ids[j], e.ids[j % 3 + 1]
                local key = math.min(a, b) .. ":" .. math.max(a, b)
                if not seen[key] then
                    seen[key] = true
                    d.segments[#d.segments + 1] = { a - 1, b - 1, -1, 0, 0, 0 }
                end
            end
        end
        -- This is an answer snapshot, not an editable reconstruction of the
        -- source model. Use export fem to retain its original angle law.
        d.boundaries = {}
        d.header.comment =
            "luafemm answer snapshot: elementwise quadrature-averaged magnetisation; use the separate .fem export to edit the source model"
    else
        for _, face in ipairs(t.faces) do
            if not face.hole then
                local label = labels[face.id] + 1
                label_current[label] = d.labels[label][5] ~= 0 and face.current / 1e6 or 0
            end
        end
    end
    local out = { F.serialize(d), "[Solution]", tostring(#m.nodes) }
    for i, p in ipairs(m.nodes) do
        out[#out + 1] = string.format(
            "%.17g\t%.17g\t%.17g\t-1",
            p[1] / m.options.unit,
            p[2] / m.options.unit,
            m.A[i]
        )
    end
    out[#out + 1] = tostring(#m.triangles)
    for index, e in ipairs(m.triangles) do
        local p, q, r = m.nodes[e.ids[1]], m.nodes[e.ids[2]], m.nodes[e.ids[3]]
        local f = require("luafemm-topology").locate(
            t,
            (p[1] + q[1] + r[1]) / (3 * m.options.unit),
            (p[2] + q[2] + r[2]) / (3 * m.options.unit)
        )
        assert(f and labels[f.id], "luafemm ANS: missing element label")
        out[#out + 1] = string.format(
            "%d\t%d\t%d\t%d",
            e.ids[1] - 1,
            e.ids[2] - 1,
            e.ids[3] - 1,
            element_labels[index] or labels[f.id]
        )
    end
    out[#out + 1] = tostring(#d.labels)
    for i = 1, #d.labels do
        out[#out + 1] = string.format("1\t%.17g", label_current[i] or 0)
    end
    out[#out + 1] = tostring(#(m.periodic_nodes or {}))
    for _, pair in ipairs(m.periodic_nodes or {}) do
        out[#out + 1] =
            string.format("%d\t%d\t%d", pair[1] - 1, pair[2] - 1, pair[3] == 1 and 0 or 1)
    end
    out[#out + 1] = "0"
    F.write_text(
        path,
        table.concat(out, "\n") .. "\n",
        m.import_document and m.import_document.filename
    )
    return d
end
return A
