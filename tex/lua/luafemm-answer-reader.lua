-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Bounded, passive non-incremental magnetic ANS reader.
-- Reconstruct interfaces from the supplied mesh, never remesh a saved solution.
-- FEMM flux values become A_phi on axisymmetric nodes; that conversion preserves
-- nodal flux, not FEMM's distinct axisymmetric interpolation between nodes.
-- @module luafemm-answer-reader
local R = {}
local F = require("luafemm-fem")
local function check(v, s)
    assert(v, "luafemm ANS: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
local function integer(v, lo, hi)
    return finite(v) and v % 1 == 0 and v >= lo and v <= hi
end
--- Parse bounded numeric sections. No Lua, TeX or FEMM expressions are run.
function R.parse(text, name)
    check(type(text) == "string" and #text <= 64 * 1024 * 1024, "input exceeds 64 MiB")
    text = text:gsub("\r\n", "\n")
    local first, last = text:find("\n%s*%[Solution%]%s*\n")
    check(first, "missing [Solution] section")
    local a = {
        document = F.parse(text:sub(1, first), name),
        nodes = {},
        elements = {},
        states = {},
        pairs = {},
        ages = {},
        element_records = {},
    }
    local lines = {}
    for line in text:sub(last + 1):gmatch("[^\n]+") do
        if line:match("%S") then
            lines[#lines + 1] = line
        end
    end
    local at = 0
    local function row(minimum, maximum)
        at = at + 1
        check(lines[at], "truncated solution at record " .. at)
        local values = {}
        for token in lines[at]:gmatch("%S+") do
            local v = tonumber(token)
            check(finite(v), "non-finite or nonnumeric solution record " .. at)
            values[#values + 1] = v
        end
        check(
            #values >= minimum and #values <= maximum,
            (a.dt and "transient DeltaT: " or "") .. "invalid solution record width at " .. at
        )
        return values
    end
    local function count()
        local n = row(1, 1)[1]
        check(integer(n, 0, 300000), "invalid solution count")
        return n
    end
    local ac = (a.document.header.frequency or 0) > 0
    local dt = a.document.header.deltat
    if dt ~= nil then
        check(
            finite(dt) and dt > 0 and not ac,
            "transient DeltaT must be positive in a real answer"
        )
        check(
            a.document.header.problemtype == "planar",
            "native transient snapshots are planar only"
        )
    end
    a.dt, a.time_previous = dt, {}
    local previous = a.document.header.prevtype or 0
    check(previous ~= 3, "previous-solution type 3 restarts are unsupported")
    check(previous == 0 or previous == 1 or previous == 2, "unknown previous-solution type")
    check(not dt or previous == 0, "transient state cannot also be an incremental bias")
    check(
        dt or previous ~= 0 or (a.document.header.prevsoln or "") == "",
        "non-incremental previous meshes are unsupported"
    )
    check(
        previous == 0 or a.document.header.problemtype == "planar",
        "response exchange is planar only"
    )
    check(not ac or previous ~= 2, "harmonic frozen answers are unsupported")
    a.previous, a.bias, a.currents = previous, {}, {}
    for i = 1, count() do
        local width = (ac and 5 or 4) + ((previous ~= 0 or dt) and 1 or 0)
        local p = row(width, width)
        if previous ~= 0 then
            a.bias[i] = p[width]
        end
        if dt then
            a.time_previous[i] = p[width]
        end
        check(integer(p[ac and 5 or 4], -1, 300000), "invalid node marker")
        a.nodes[i] = p
    end
    check(#a.nodes >= 3, "empty solution mesh")
    for i = 1, count() do
        local e = row((previous ~= 0 or dt) and 8 or 4, 8)
        -- Ordinary DC companion meshes may retain edge markers and a current
        -- column for a later incremental solve; those are not transient state.
        check(#e == 4 or #e == 7 or #e == 8, "unsupported element record")
        if previous ~= 0 then
            a.currents[i] = e[8] * 1e6
        end
        for j = 1, 3 do
            check(integer(e[j], 0, #a.nodes - 1), "invalid element node index")
        end
        check(integer(e[4], 0, #a.document.labels - 1), "invalid element label index")
        for j = 5, math.min(#e, 7) do
            check(integer(e[j], -300000, 300000), "invalid element edge marker")
        end
        check(e[1] ~= e[2] and e[1] ~= e[3] and e[2] ~= e[3], "repeated triangle node")
        local p, q, r = a.nodes[e[1] + 1], a.nodes[e[2] + 1], a.nodes[e[3] + 1]
        check(
            (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1]) > 0,
            "degenerate or clockwise solution element"
        )
        a.elements[i] = { e[1], e[2], e[3], e[4] }
        a.element_records[i] = e
    end
    check(#a.elements > 0, "empty solution elements")
    check(count() == #a.document.labels, "label state count does not match embedded problem")
    for i = 1, #a.document.labels do
        local c = row(ac and 3 or 2, ac and 3 or 2)
        check(c[1] == 0 or c[1] == 1, "unknown label current/voltage case")
        a.states[i] = c
    end
    for i = 1, count() do
        local p = row(3, 3)
        check(
            integer(p[1], 0, #a.nodes - 1)
                and integer(p[2], 0, #a.nodes - 1)
                and (p[3] == 0 or p[3] == 1),
            "invalid periodic node record"
        )
        a.pairs[i] = p
    end
    for i = 1, count() do
        at = at + 1
        local label = lines[at] and lines[at]:match('^"([^"\r\n]+)"$')
        check(label, "invalid AGE name")
        local header = row(11, 11)
        check(header[1] == 0 or header[1] == 1, "invalid AGE periodicity")
        check(
            header[4] > 0 and header[5] > header[4] and header[6] > 0 and header[6] <= 360 + 1e-8,
            "invalid AGE dimensions"
        )
        check(integer(header[9], 1, 100000), "invalid AGE segment count")
        local age = { name = label, header = header, quads = {} }
        for j = 1, header[9] + 1 do
            local q = row(8, 8)
            for k = 1, 7, 2 do
                check(integer(q[k], 0, #a.nodes - 1), "invalid AGE node index")
            end
            age.quads[j] = q
        end
        a.ages[i] = age
    end
    if lines[at + 1] and lines[at + 1]:lower() == "[agerings]" then
        check(#a.ages > 0, "AgeRings requires an AGE block")
        at = at + 1
        for _, age in ipairs(a.ages) do
            check(count() == age.header[9], "AgeRings count mismatch")
            age.rings = {}
            for j = 1, age.header[9] do
                local p = row(2, 2)
                check(
                    integer(p[1], 0, #a.nodes - 1) and integer(p[2], 0, #a.nodes - 1),
                    "invalid raw AGE ring node"
                )
                age.rings[j] = p
            end
        end
    end
    check(at == #lines, "unexpected trailing solution records")
    return a
end
--- Read at most 64 MiB plus one byte so oversize inputs cannot be truncated silently.
function R.read(path)
    local stream, message = io.open(path, "rb")
    check(stream, "cannot open " .. tostring(path) .. ": " .. tostring(message))
    local text = stream:read(64 * 1024 * 1024 + 1)
    stream:close()
    local a = R.parse(text, path)
    -- A time answer owns DeltaT, A_previous and its current terminal state.
    -- Its solver-control sibling is neither needed nor executed. Ordinary
    -- problem/answer exchange retains the pre-existing companion guard.
    if not a.dt then
        F.check_companion(path)
    end
    return a
end
-- Reconstruct interfaces from label changes and exposed edges. Original segment
-- properties are assigned only where both endpoints lie on the same source
-- segment/chord. This also retains internal prescribed-potential interfaces.
local function mesh_document(a)
    local d = F.copy(a.document)
    if #a.ages > 0 then
        -- Snapshot constraints are verified against their stored nodal pairs;
        -- they are not reconstructed as new, potentially different DOFs.
        for id, b in ipairs(d.boundaries) do
            if b.bdrytype == 4 or b.bdrytype == 5 or b.bdrytype == 6 or b.bdrytype == 7 then
                for _, s in ipairs(d.segments) do
                    if s[4] == id then
                        s[4] = 0
                    end
                end
                for _, arc in ipairs(d.arcs) do
                    if arc[5] == id then
                        arc[5] = 0
                    end
                end
                b.bdrytype, b.c0, b.c1 = 2, 0, 0
            end
        end
        for _, age in ipairs(a.ages) do
            local h = age.header
            local theta = math.rad(h[6] / 2)
            local found = false
            local holes = {}
            for _, p in ipairs(d.holes) do
                local r = math.sqrt((p[1] - h[7]) ^ 2 + (p[2] - h[8]) ^ 2)
                local inside = r > h[4] and r < h[5]
                found = found or inside
                -- In a reduced mesh the missing strip opens to the radial
                -- exterior: it is not a bounded arrangement face. Its native
                -- mesher seed must not be reapplied as an interior hole.
                if not inside or h[6] > 360 - 1e-8 then
                    holes[#holes + 1] = p
                end
            end
            d.holes = holes
            if not found and h[6] > 360 - 1e-8 then
                d.holes[#d.holes + 1] = {
                    h[7] + (h[4] + h[5]) / 2 * math.cos(theta),
                    h[8] + (h[4] + h[5]) / 2 * math.sin(theta),
                    0,
                }
            end
        end
    end
    -- Native answers retain the arc subdivision used by the mesher. Rebuild
    -- those chords; refinement nodes may lie inside a chord, not on the circle.
    local native_arcs = #a.ages > 0 and F.copy(d.arcs) or {}
    for _, arc in ipairs(#a.ages > 0 and {} or d.arcs) do
        local a, b = d.points[arc[1] + 1], d.points[arc[2] + 1]
        local theta = math.rad(arc[3])
        local dx, dy = b[1] - a[1], b[2] - a[2]
        local cx, cy =
            (a[1] + b[1]) / 2 - dy / (2 * math.tan(theta / 2)),
            (a[2] + b[2]) / 2 + dx / (2 * math.tan(theta / 2))
        local radius = math.sqrt(dx * dx + dy * dy) / (2 * math.sin(theta / 2))
        local start = math.atan(a[2] - cy, a[1] - cx)
        local n = math.ceil(arc[3] / arc[4])
        check(n <= 40000, "arc subdivision limit exceeded")
        local previous = arc[1]
        for i = 1, n do
            local id = arc[2]
            if i < n then
                id = #d.points
                d.points[id + 1] = {
                    cx + radius * math.cos(start + theta * i / n),
                    cy + radius * math.sin(start + theta * i / n),
                    0,
                    0,
                }
            end
            d.segments[#d.segments + 1] = { previous, id, -1, arc[5], arc[6], arc[7] }
            previous = id
        end
    end
    d.arcs = {}
    local top = require("luafemm-topology")
    local extent = 0
    for _, p in ipairs(a.nodes) do
        extent = math.max(extent, math.abs(p[1]), math.abs(p[2]))
    end
    local tolerance = math.max(extent, 1e-9) * 1e-9
    local edges, used, first = {}, {}, {}
    for _, e in ipairs(a.elements) do
        first[e[4] + 1] = first[e[4] + 1] or e
        for j = 1, 3 do
            local u, v = e[j] + 1, e[j % 3 + 1] + 1
            used[u] = true
            if u > v then
                u, v = v, u
            end
            local key = u .. ":" .. v
            local edge = edges[key]
            if edge then
                check(edge.count == 1, "non-manifold solution edge")
                edge.count, edge.interface = 2, edge.label ~= e[4]
            else
                edges[key] = { u = u, v = v, count = 1, label = e[4] }
            end
        end
    end
    local oldpoints, oldsegments = d.points, d.segments
    d.points, d.segments = {}, {}
    for i, p in ipairs(a.nodes) do
        check(used[i], "unused solution node")
        d.points[i] = { p[1], p[2], 0, 0 }
        for _, q in ipairs(oldpoints) do
            if
                q[3] ~= 0
                and math.abs(q[1] - p[1]) <= tolerance
                and math.abs(q[2] - p[2]) <= tolerance
            then
                d.points[i][3] = q[3]
            end
        end
    end
    local sorted = {}
    for _, e in pairs(edges) do
        sorted[#sorted + 1] = e
    end
    table.sort(sorted, function(x, y)
        return x.u < y.u or x.u == y.u and x.v < y.v
    end)
    local exterior, covered = {}, {}
    local arcs_covered = {}
    for _, e in ipairs(sorted) do
        local boundary, group = 0, 0
        for sid, s in ipairs(oldsegments) do
            if
                s[4] ~= 0
                and top.distance(a.nodes[e.u], oldpoints[s[1] + 1], oldpoints[s[2] + 1]) <= tolerance
                and top.distance(a.nodes[e.v], oldpoints[s[1] + 1], oldpoints[s[2] + 1])
                    <= tolerance
            then
                check(boundary == 0 or boundary == s[4], "ambiguous solution boundary")
                boundary, group = s[4], s[6]
                covered[sid] = true
            end
        end
        for aid, arc in ipairs(native_arcs) do
            if arc[5] ~= 0 then
                local p, q = oldpoints[arc[1] + 1], oldpoints[arc[2] + 1]
                local theta = math.rad(arc[3])
                local cx = (p[1] + q[1]) / 2 - (q[2] - p[2]) / (2 * math.tan(theta / 2))
                local cy = (p[2] + q[2]) / 2 + (q[1] - p[1]) / (2 * math.tan(theta / 2))
                local radius = math.sqrt((p[1] - cx) ^ 2 + (p[2] - cy) ^ 2)
                local start = math.atan(p[2] - cy, p[1] - cx)
                local function on(id)
                    local v = a.nodes[id]
                    local angle = (math.atan(v[2] - cy, v[1] - cx) - start) % (2 * math.pi)
                    if 2 * math.pi - angle < 1e-9 then
                        angle = 0
                    end
                    return math.abs(math.sqrt((v[1] - cx) ^ 2 + (v[2] - cy) ^ 2) - radius)
                            <= tolerance
                        and angle <= theta + 1e-9
                end
                if on(e.u) and on(e.v) then
                    check(boundary == 0 or boundary == arc[5], "ambiguous native arc boundary")
                    boundary, group = arc[5], arc[7]
                    arcs_covered[aid] = true
                end
            end
        end
        if e.count == 1 or e.interface or boundary ~= 0 then
            d.segments[#d.segments + 1] = { e.u - 1, e.v - 1, -1, boundary, 0, group }
        end
        if e.count == 1 then
            exterior[e.u], exterior[e.v] = true, true
        end
    end
    for sid, segment in ipairs(oldsegments) do
        check(
            segment[4] == 0 or covered[sid],
            "a prescribed boundary segment has no matching solution edge; unsupported arc discretisation or incomplete mesh"
        )
    end
    for aid, arc in ipairs(native_arcs) do
        check(arc[5] == 0 or arcs_covered[aid], "native arc boundary has no matching solution edge")
    end
    for i, label in ipairs(d.labels) do
        local e = first[i]
        check(e, "unused block label in answer")
        label[1], label[2] = 0, 0
        for j = 1, 3 do
            local p = a.nodes[e[j] + 1]
            label[1], label[2] = label[1] + p[1] / 3, label[2] + p[2] / 3
        end
    end
    return d, exterior
end
--- Restore owned mesh/field data with physical metadata fixed by the file.
-- Mesh size has no effect on a snapshot; only units, logging, cache policy and
-- contour rendering controls may be overridden. Imported answers bypass caches.
function R.model(a, options)
    if a.dt then
        return require("luafemm-time-exchange").model(a, options)
    end
    options = options or {}
    local native = a.ages and #a.ages > 0
    if native then
        check(
            options.age_interpolation == "femm4" or options.age_interpolation == "femm5",
            "native AGE snapshots require age interpolation=femm4 or femm5"
        )
        check(
            a.document.header.problemtype == "planar" and (a.previous or 0) == 0,
            "native AGE snapshots require ordinary planar static or harmonic solutions"
        )
    end
    for key in pairs(options) do
        check(
            key == "unit"
                or key == "log"
                or key == "cache"
                or key == "contour_subdivisions"
                or key == "age_interpolation",
            "answer snapshots own their settings; cannot override " .. key
        )
    end
    local d, exterior = mesh_document(a)
    local M = require("luafemm")
    local response = a.previous and a.previous ~= 0
    local model_options = F.copy(options)
    model_options.age_interpolation = nil
    local m = F.model(d, model_options, response, #a.ages == 0 and a.pairs or nil)
    local nodes, elements, values, imaginary = {}, {}, {}, {}
    local unit = F.units[tostring(d.header.lengthunits):lower()]
    local axi, ac = m.options.problem_type == "axisymmetric", m.options.frequency > 0
    for i, p in ipairs(a.nodes) do
        nodes[i] = { p[1] * unit, p[2] * unit }
        local factor = axi and 2 * math.pi * nodes[i][1] or 1
        if factor == 0 then
            check(
                math.abs(p[3]) < 1e-14 and (not ac or math.abs(p[4]) < 1e-14),
                "nonzero flux on axis"
            )
            values[i], imaginary[i] = 0, 0
        else
            values[i], imaginary[i] = p[3] / factor, ac and p[4] / factor or 0
        end
    end
    for i, e in ipairs(a.elements) do
        elements[i] = { e[1] + 1, e[2] + 1, e[3] + 1 }
    end
    if response then
        for i = 1, #nodes do
            check(finite(a.bias[i]), "missing or invalid embedded bias potential")
        end
        for i = 1, #elements do
            check(finite(a.currents[i]), "missing or invalid embedded bias current")
        end
        m.operating = {
            mode = a.previous == 1 and "incremental" or "frozen",
            A = F.copy(a.bias),
            nodes = F.copy(nodes),
            elements = F.copy(elements),
            currents = F.copy(a.currents),
            provenance = "embedded ANS bias",
        }
    end
    if native then
        m.import_answer = { snapshot = true }
        m.native_age = { ages = F.copy(a.ages), edition = options.age_interpolation, unit = unit }
    end
    M.prepare_mesh(m, nodes, elements, exterior)
    -- Quality is measured on the restored mesh. Zero work counters describe
    -- this importer (no flips/refinement), not the unknown upstream mesh run.
    local quality = {
        imported = true,
        boundary_nodes = F.copy(exterior),
        segments = {},
        seconds = 0,
        flips = 0,
        min_angle = 180,
        max_area = 0,
        mesh_tolerance = m.topology.tolerance,
        exact_orientation = 0,
        exact_incircle = 0,
        exact_diametral = 0,
        steiner_points = 0,
        segment_splits = 0,
        circumcenters = 0,
        spatial_queries = 0,
        spatial_candidates = 0,
    }
    for i, segment in ipairs(d.segments) do
        quality.segments[i] = { segment[1] + 1, segment[2] + 1 }
    end
    for _, e in ipairs(m.triangles) do
        quality.max_area = math.max(quality.max_area, e.area / m.options.unit ^ 2)
        for j = 1, 3 do
            local p, q, r = nodes[e.ids[j]], nodes[e.ids[j % 3 + 1]], nodes[e.ids[(j + 1) % 3 + 1]]
            local dot = (q[1] - p[1]) * (r[1] - p[1]) + (q[2] - p[2]) * (r[2] - p[2])
            quality.min_angle = math.min(quality.min_angle, math.deg(math.atan(2 * e.area, dot)))
        end
    end
    m.mesh_stats = quality
    for i, e in ipairs(m.triangles) do
        check(
            m.topology.faces[e.face].label == a.elements[i][4] + 1,
            "element label does not agree with reconstructed topology"
        )
    end
    local norm = 1e-14
    for _, v in ipairs(values) do
        norm = math.max(norm, math.abs(v))
    end
    if ac then
        for _, v in ipairs(imaginary) do
            norm = math.max(norm, math.abs(v))
        end
    end
    for id, v in pairs(m.fixed_values) do
        check(
            math.abs(values[id] - v) <= 1e-7 * norm
                and (not ac or math.abs(imaginary[id]) <= 1e-7 * norm),
            "answer violates prescribed nodal potential"
        )
    end
    m.A, m.A_im = values, ac and imaginary or nil
    -- Store exact passive records in the physical signature. A filename alone
    -- must never identify a bias after the source answer has been replaced.
    m.import_answer = {
        nodes = F.copy(a.nodes),
        elements = F.copy(a.elements),
        states = F.copy(a.states),
        previous = a.previous,
        bias = F.copy(a.bias),
        currents = F.copy(a.currents),
        native = native and F.copy(a) or nil,
        age_interpolation = native and options.age_interpolation or nil,
    }
    for _, p in ipairs(a.pairs) do
        local sign = p[3] == 0 and 1 or -1
        local sa, sb = axi and nodes[p[1] + 1][1] or 1, axi and nodes[p[2] + 1][1] or 1
        local scale = math.max(sa, sb, 1e-30)
        check(
            math.abs(sa * values[p[1] + 1] - sign * sb * values[p[2] + 1]) <= 1e-7 * norm * scale,
            "answer violates a periodic nodal constraint"
        )
        if ac then
            check(
                math.abs(sa * imaginary[p[1] + 1] - sign * sb * imaginary[p[2] + 1])
                    <= 1e-7 * norm * scale,
                "answer violates an imaginary periodic nodal constraint"
            )
        end
    end
    m.stats =
        { imported = true, nodes = #nodes, triangles = #elements, seconds = 0, newton = 0, cg = 0 }
    if response then
        m.operating_contract = require("luafemm-cache").signature(m)
    end
    M.update_fields(m)
    -- A nearly cancelling induced voltage needs an absolute term referenced
    -- to omega*A, not to the tiny terminal difference alone. FEMM's iterative
    -- precision is relative to the full field problem. Cap that allowance so
    -- a loose file precision cannot suppress circuit consistency checks.
    local induced = {}
    if ac then
        local omega = 2 * math.pi * m.options.frequency
        for _, e in ipairs(m.triangles) do
            local label = m.topology.faces[e.face].label
            for _, id in ipairs(e.ids) do
                local radius = axi and m.nodes[id][1] or 1
                local value = omega * radius * (math.abs(m.A[id]) + math.abs(m.A_im[id]))
                check(finite(value), "induced-voltage scale exceeds finite range")
                induced[label] = math.max(induced[label] or 0, value)
            end
        end
    end
    for _, face in ipairs(m.topology.faces) do
        if not face.hole then
            local state, label = a.states[face.label], d.labels[face.label]
            local kind, re, im
            if ac then
                kind, re, im = require("luafemm-harmonic-exchange").answer_state(m, face, label)
            else
                kind, re, im = 1, label[5] ~= 0 and face.current / 1e6 or 0, 0
                if state[1] == 0 then
                    local sigma = m.materials[face.material].conductivity
                    check(
                        label[5] ~= 0 and sigma > 0 and not axi,
                        "static voltage-gradient answers require planar conducting circuits"
                    )
                    kind, re = 0, -face.current / sigma
                end
            end
            check(state[1] == kind, "unsupported or inconsistent label current/voltage case")
            local scale = math.max(
                math.abs(re),
                math.abs(im),
                math.abs(state[2]),
                math.abs(state[3] or 0),
                1e-14
            )
            local allowance = ac
                    and kind == 0
                    and 32 * math.min(m.options.tolerance, 1e-8) * (induced[face.label] or 0)
                or 0
            check(
                math.abs(re - state[2]) + math.abs(im - (state[3] or 0)) <= 1e-5 * scale + allowance,
                "label current/voltage data disagree with reconstructed circuit equations"
            )
        end
    end
    return m
end
return R
