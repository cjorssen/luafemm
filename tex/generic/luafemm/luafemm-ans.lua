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
    for _, e in ipairs(m.triangles) do
        local p, q, r = m.nodes[e.ids[1]], m.nodes[e.ids[2]], m.nodes[e.ids[3]]
        local f = require("luafemm-topology").locate(
            t,
            (p[1] + q[1] + r[1]) / (3 * m.options.unit),
            (p[2] + q[2] + r[2]) / (3 * m.options.unit)
        )
        assert(f and labels[f.id], "luafemm ANS: missing element label")
        out[#out + 1] =
            string.format("%d\t%d\t%d\t%d", e.ids[1] - 1, e.ids[2] - 1, e.ids[3] - 1, labels[f.id])
    end
    out[#out + 1] = tostring(#d.labels)
    for _ = 1, #d.labels do
        out[#out + 1] = "1\t0"
    end
    out[#out + 1] = "0"
    out[#out + 1] = "0"
    F.write_text(
        path,
        table.concat(out, "\n") .. "\n",
        m.import_document and m.import_document.filename
    )
    return d
end
return A
