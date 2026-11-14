-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Owned static operating points and same-mesh frozen/incremental responses.
-- No field interpolation between meshes is performed. The bias potential,
-- geometry and source declarations are copied; the original model stays intact.
-- @module luafemm-operating
local O = {}
local function check(v, s)
    assert(v, "luafemm operating point: " .. s)
end
local function copy(v, memo)
    if type(v) ~= "table" then
        return v
    end
    memo = memo or {}
    if memo[v] then
        return memo[v]
    end
    local t = {}
    memo[v] = t
    for k, x in pairs(v) do
        t[k] = copy(x, memo)
    end
    return t
end
local function eligible(m)
    check(m and m.stats, "solve the bias problem first")
    check(not m.sliding, "sliding-position operating points are not validated yet")
    check(not m.operating, "a response cannot become a new nonlinear operating point")
    check(
        m.options.frequency == 0 and m.options.problem_type ~= "axial",
        "a static magnetic bias is required"
    )
    check(m.options.field_model == "open", "ideal confined circuits are unsupported")
    check(not m.import_document, "imported operating points are not supported yet")
    check(not m.options.source and not m.options.boundary, "callbacks cannot be captured")
    check(m.solution_signature == require("luafemm-cache").signature(m), "stale bias solution")
    for _, c in pairs(m.circuits or {}) do
        check(
            c.type == "series" and c.drive == "current",
            "bias circuits require series current drive"
        )
    end
end
--- Capture an independent solved DC model, including the ordered mesh.
-- @tparam table m Converged, unchanged native model.
-- @treturn table Owned snapshot suitable for repeated O.new calls.
function O.capture(m)
    eligible(m)
    local out = {}
    for _, key in ipairs({
        "options",
        "materials",
        "regions",
        "boundaries",
        "frontiers",
        "domain",
        "pairs",
        "topology",
        "native_topology",
        "circuits",
        "nodes",
        "triangles",
        "exterior_nodes",
        "mesh_stats",
        "xs",
        "ys",
        "A",
        "stats",
        "solution_signature",
        "frame",
        "component_basis",
    }) do
        out[key] = copy(m[key])
    end
    return out
end
--- Derive a static or small-signal harmonic response on the exact bias mesh.
-- Frozen defaults retain all bias sources. Incremental defaults zero every
-- independent source (including coercivity and boundary loads). Explicit
-- currents replace the selected circuit currents, in A, in either mode.
-- @tparam table bias Owned capture or unchanged converged native model.
-- @tparam[opt] table o mode, currents (A), frequency (Hz), phases (degrees), cache.
-- @treturn table Unsolved response; original material declarations are retained.
function O.new(bias, o)
    o = o or {}
    for k in pairs(o) do
        check(
            k == "mode" or k == "currents" or k == "cache" or k == "frequency" or k == "phases",
            "unknown option " .. tostring(k)
        )
    end
    local m = O.capture(bias)
    local mode = o.mode or "incremental"
    check(mode == "frozen" or mode == "incremental", "unknown mode")
    local frequency = o.frequency or 0
    check(
        type(frequency) == "number"
            and frequency == frequency
            and frequency >= 0
            and frequency < math.huge,
        "frequency must be finite and nonnegative"
    )
    check(
        frequency == 0 or mode == "incremental",
        "harmonic responses require permeability=incremental"
    )
    m.options.frequency = frequency
    local nodes, elements = {}, {}
    for i, p in ipairs(m.nodes) do
        nodes[i] = { p[1], p[2] }
    end
    for i, e in ipairs(m.triangles) do
        elements[i] = copy(e.ids)
    end
    m.operating = {
        mode = mode,
        A = copy(m.A),
        nodes = copy(nodes),
        elements = copy(elements),
        provenance = m.solution_signature,
    }
    m.nodes, m.triangles, m.A, m.stats, m.solution_signature = nil, nil, nil, nil, nil
    m.options.cache = copy(o.cache or { mode = "off" })
    require("luafemm-cache").configure(m)
    if mode == "incremental" then
        for _, r in ipairs(m.regions) do
            r.current = 0
        end
        for _, c in pairs(m.circuits or {}) do
            c.current = 0
        end
        local function homogeneous(c)
            c.value, c.dx, c.dy = 0, 0, 0
        end
        for _, c in pairs(m.boundaries) do
            homogeneous(c)
        end
        for _, f in pairs(m.frontiers or {}) do
            homogeneous(f.condition)
        end
    end
    check(type(o.currents or {}) == "table", "currents must be a table")
    for name, value in pairs(o.currents or {}) do
        check(m.circuits and m.circuits[name], "unknown circuit " .. tostring(name))
        check(
            type(value) == "number" and value == value and math.abs(value) < math.huge,
            "invalid circuit current"
        )
        m.circuits[name].current = value
    end
    check(type(o.phases or {}) == "table", "phases must be a table")
    for name, value in pairs(o.phases or {}) do
        check(m.circuits and m.circuits[name], "unknown circuit " .. tostring(name))
        check(
            type(value) == "number" and value == value and math.abs(value) < math.huge,
            "invalid circuit phase"
        )
        check(frequency > 0 or value == 0, "circuit phase requires frequency>0")
        m.circuits[name].phase = value
    end
    -- Recompile loads on the same geometry, then reconstruct all DOFs and
    -- element properties. Do not retain a mesh-owned reference from the bias.
    require("luafemm-frontiers").compile(m)
    if m.topology and not (m.domain or m.frontiers or m.circuits) then
        for _, f in ipairs(m.topology.faces) do
            if f.region then
                f.current = m.regions[f.region].current
            end
        end
    end
    require("luafemm").prepare_mesh(m, nodes, elements, m.exterior_nodes)
    m.operating_contract = require("luafemm-cache").signature(m)
    if frequency > 0 then
        require("luafemm-harmonic").validate(m)
    end
    return m
end
--- Reject changed declarations before either solving or reading a cache.
function O.validate(m)
    local bias = m.operating
    check(
        m.nodes and #m.nodes == #bias.nodes and #m.triangles == #bias.elements,
        "response must retain the bias mesh"
    )
    for i, p in ipairs(m.nodes) do
        check(p[1] == bias.nodes[i][1] and p[2] == bias.nodes[i][2], "response mesh changed")
    end
    for i, e in ipairs(m.triangles) do
        for j, id in ipairs(e.ids) do
            check(id == bias.elements[i][j], "response mesh changed")
        end
    end
    check(
        m.options.frequency == 0 or m.operating.mode == "incremental",
        "harmonic responses require permeability=incremental"
    )
    check(
        m.operating_contract == require("luafemm-cache").signature(m),
        "response declarations changed; create a new response from the bias"
    )
end
--- Effective scaled reluctivity tensor at a point on the original mesh.
-- Axisymmetric tensors use the local bias curl, including A_phi/r; the
-- assembly evaluates them at every quadrature point, not at the centroid.
-- @tparam table m Response model.
-- @tparam table e Response triangle.
-- @tparam[opt] number x SI radius (required for axisymmetric evaluation).
-- @tparam[opt] number y SI axial position.
-- @return xx, xy, yy entries of mu0*dH/dB (or the frozen secant tensor).
function O.tensor(m, e, x, y)
    if not e.p.bh then
        return table.unpack(e.p.tensor)
    end
    local bx, by = 0, 0
    local A = m.operating.A
    if m.options.problem_type == "axisymmetric" then
        local p = m.nodes[e.ids[1]]
        local v = e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2])
        local w = e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
        local n, a, ar, az = { 1 - v - w, v, w }, 0, 0, 0
        for i, id in ipairs(e.ids) do
            a, ar, az = a + n[i] * A[id], ar + e.gx[i] * A[id], az + e.gy[i] * A[id]
        end
        bx, by = x == 0 and 0 or -az, ar + (x > 0 and a / x or ar)
    else
        for i, id in ipairs(e.ids) do
            bx, by = bx + e.gy[i] * A[id], by - e.gx[i] * A[id]
        end
    end
    if m.operating.mode == "frozen" then
        local nu = require("luafemm").constitutive(e.p, math.sqrt(bx * bx + by * by))
        return nu, 0, nu
    end
    local _, _, xx, xy, yy = require("luafemm-static").law(e.p, bx, by)
    return xx, xy, yy
end
return O
