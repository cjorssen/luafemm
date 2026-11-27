-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Static, harmonic and supported time answer snapshots; passive import/export.
-- Coordinates use file units; nodal values are T m planar, Wb axisymmetric.
-- Element label references are zero-based.
-- @module luafemm-ans
local A = {}
--- Preserve supported native answer records, including FEMM 4 AGE quad records
-- and optional FEMM 5 raw ring lists. No field conversion or recalculation is
-- performed. Omitting raw rings disables subsequent angle re-pairing in FEMM 5.
-- Unknown trailing data fail during parsing rather than being silently dropped.
function A.copy(input, output, rings)
    local F = require("luafemm-fem")
    F.check_companion(input)
    F.check_companion(output)
    assert(input ~= output, "luafemm ANS: input and output must differ")
    rings = rings or "preserve"
    assert(rings == "preserve" or rings == "omit", "luafemm ANS: invalid raw-ring policy")
    local reader = require("luafemm-answer-reader")
    local a = reader.read(input)
    local text = A.serialize(a, rings)
    reader.parse(text, output)
    local file = assert(io.open(output, "wb"), "luafemm ANS: cannot write " .. output)
    assert(file:write(text))
    assert(file:close())
    return a
end
--- Serialize already validated passive records without changing their meaning.
function A.serialize(a, rings)
    local F = require("luafemm-fem")
    rings = rings or "preserve"
    local out = { F.serialize(a.document), "[Solution]" }
    local function row(values)
        local tokens = {}
        for _, v in ipairs(values) do
            tokens[#tokens + 1] = string.format("%.17g", v)
        end
        out[#out + 1] = table.concat(tokens, " ")
    end
    local function section(values)
        row({ #values })
        for _, v in ipairs(values) do
            row(v)
        end
    end
    section(a.nodes)
    section(a.element_records)
    section(a.states)
    section(a.pairs)
    row({ #a.ages })
    local raw = #a.ages > 0
    for _, age in ipairs(a.ages) do
        out[#out + 1] = '"' .. age.name .. '"'
        row(age.header)
        for _, q in ipairs(age.quads) do
            row(q)
        end
        raw = raw and age.rings ~= nil
    end
    if raw and rings == "preserve" then
        out[#out + 1] = "[AgeRings]"
        for _, age in ipairs(a.ages) do
            section(age.rings)
        end
    end
    return table.concat(out, "\n") .. "\n"
end
--- Write the existing mesh and potentials, without solving or remeshing.
-- Planar linear material post-processing agrees exactly in the represented model;
-- nonlinear H agreement requires interpolation='femm' at material preparation.
-- Axisymmetric export preserves nodal flux; FEMM interpolation changes element fields.
-- @tparam table m Converged model, including one restored from a valid cache.
-- @tparam string path Destination including .ans.
-- @treturn table Embedded computational problem document.
function A.write(m, path, internal)
    if m.transient then
        require("luafemm-time-exchange").validate(m)
    end
    assert(
        not m.initial_time,
        "luafemm ANS: an initial time field is not an equilibrium answer; save its time state"
    )
    if m.operating and not (internal and internal.previous) then
        return require("luafemm-response-exchange").write(
            m,
            path,
            "ans",
            internal and internal.bias_file
        )
    end
    require("luafemm-fem").check_companion(path)
    if m.sliding then
        require("luafemm-sliding").validate(m)
    end
    if m.operating then
        require("luafemm-response-exchange").validate(m)
    end
    assert(
        not require("luafemm-airgap").active(m) or m.native_age,
        "luafemm ANS: analytical gaps cannot be represented by P1 answer records"
    )
    local ac = m.options.frequency > 0
    assert(m and m.stats, "luafemm ANS: a converged solution is required")
    assert(
        m.options.field_model ~= "ideal",
        "luafemm ANS: ideal floating boundaries cannot be exported"
    )
    assert(
        m.solution_signature == require("luafemm-cache").signature(m),
        "luafemm ANS: stale solution; rebuild and solve the changed model"
    )
    if m.import_answer and m.import_answer.native then
        local F = require("luafemm-fem")
        F.write_text(
            path,
            A.serialize(m.import_answer.native),
            m.import_document and m.import_document.filename
        )
        return F.copy(m.import_answer.native.document)
    end
    for _, e in ipairs(m.triangles) do
        assert(
            not e.p.bh or e.p.interpolation == "femm",
            "luafemm ANS: nonlinear export requires interpolation=femm before solving"
        )
    end
    local F = require("luafemm-fem")
    local d, t, labels = F.document(m, true)
    if m.transient then
        d.header.deltat = m.transient.dt
    end
    if internal and internal.previous then
        d.header.prevsoln = internal.previous
    end
    -- The previous mesh carries boundary indices of the response problem.
    -- Geometry is unchanged, but incremental loads are homogeneous.
    local meshdoc = internal and internal.mesh_document or d
    local top = require("luafemm-topology")
    local function edge_marker(a, b)
        local u = m.options.unit
        a, b = { a[1] / u, a[2] / u }, { b[1] / u, b[2] / u }
        for _, segment in ipairs(meshdoc.segments) do
            if
                segment[4] ~= 0
                and top.distance(a, meshdoc.points[segment[1] + 1], meshdoc.points[segment[2] + 1]) <= t.tolerance * 30
                and top.distance(
                        b,
                        meshdoc.points[segment[1] + 1],
                        meshdoc.points[segment[2] + 1]
                    )
                    <= t.tolerance * 30
            then
                return segment[4] - 1
            end
        end
        return -1
    end
    local element_labels, label_current, label_states = {}, {}, {}
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
            if not face.hole and labels[face.id] ~= nil then
                local label = labels[face.id] + 1
                label_current[label] = d.labels[label][5] ~= 0 and face.current / 1e6 or 0
                if ac then
                    label_states[label] = {
                        require("luafemm-harmonic-exchange").answer_state(m, face, d.labels[label]),
                    }
                end
            end
        end
    end
    local out = { F.serialize(d), "[Solution]", tostring(#m.nodes) }
    for i, p in ipairs(m.nodes) do
        local factor = m.options.problem_type == "axisymmetric" and 2 * math.pi * p[1] or 1
        local row = string.format(
            "%.17g\t%.17g\t%.17g",
            p[1] / m.options.unit,
            p[2] / m.options.unit,
            factor * m.A[i]
        )
        if ac then
            row = row .. string.format("\t%.17g", factor * m.A_im[i])
        end
        row = row .. "\t-1"
        if m.operating then
            row = row .. string.format("\t%.17g", factor * m.operating.A[i])
        end
        if m.transient then
            row = row .. string.format("\t%.17g", factor * m.transient.previous[i])
        end
        out[#out + 1] = row
    end
    out[#out + 1] = tostring(m.native_age and m.sliding.bulk_count or #m.triangles)
    for index, e in ipairs(m.triangles) do
        if not e.analytical_gap then
            local p, q, r = m.nodes[e.ids[1]], m.nodes[e.ids[2]], m.nodes[e.ids[3]]
            local f = require("luafemm-topology").locate(
                t,
                (p[1] + q[1] + r[1]) / (3 * m.options.unit),
                (p[2] + q[2] + r[2]) / (3 * m.options.unit)
            )
            assert(f and labels[f.id], "luafemm ANS: missing element label")
            local row = string.format(
                "%d\t%d\t%d\t%d",
                e.ids[1] - 1,
                e.ids[2] - 1,
                e.ids[3] - 1,
                element_labels[index] or labels[f.id]
            )
            if internal or m.operating or m.native_age or m.transient then
                row = row
                    .. string.format(
                        "\t%d\t%d\t%d\t%.17g",
                        edge_marker(p, q),
                        edge_marker(q, r),
                        edge_marker(r, p),
                        (m.operating and m.operating.currents[index] or e.current) / 1e6
                    )
            end
            out[#out + 1] = row
        end
    end
    out[#out + 1] = tostring(#d.labels)
    for i = 1, #d.labels do
        if ac then
            out[#out + 1] = string.format("%d\t%.17g\t%.17g", table.unpack(label_states[i]))
        else
            out[#out + 1] = string.format("1\t%.17g", label_current[i] or 0)
        end
    end
    out[#out + 1] = tostring(#(m.periodic_nodes or {}))
    for _, pair in ipairs(m.periodic_nodes or {}) do
        out[#out + 1] =
            string.format("%d\t%d\t%d", pair[1] - 1, pair[2] - 1, pair[3] == 1 and 0 or 1)
    end
    if m.native_age then
        local function row(v)
            local words = {}
            for _, x in ipairs(v) do
                words[#words + 1] = string.format("%.17g", x)
            end
            out[#out + 1] = table.concat(words, " ")
        end
        row({ #m.native_age.ages })
        for _, a in ipairs(m.native_age.ages) do
            out[#out + 1] = '"' .. a.name .. '"'
            row(a.header)
            for _, q in ipairs(a.quads) do
                row(q)
            end
        end
        if m.native_age.edition == "femm5" then
            out[#out + 1] = "[AgeRings]"
            for _, a in ipairs(m.native_age.ages) do
                row({ #a.rings })
                for _, q in ipairs(a.rings) do
                    row(q)
                end
            end
        end
    else
        out[#out + 1] = "0"
    end
    F.write_text(
        path,
        table.concat(out, "\n") .. "\n",
        m.import_document and m.import_document.filename
    )
    return d
end
--- Read a passive answer, including the embedded FEM document and exact mesh.
function A.read(path)
    return require("luafemm-answer-reader").read(path)
end
--- Parse answer bytes without evaluating imported expressions.
function A.parse(text, name)
    return require("luafemm-answer-reader").parse(text, name)
end
--- Restore a snapshot without solving; unsupported physical records raise.
function A.model(answer, options)
    return require("luafemm-answer-reader").model(answer, options)
end
return A
