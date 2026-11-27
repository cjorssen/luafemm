-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Owned time states and passive restart persistence, independent of TeX.
-- A state contains only one previous step, never a chain of model references.
-- The passive cache codec supplies bounded parsing and atomic file replacement.
local S = {}
local function check(v, s)
    assert(v, "luafemm time state: " .. s)
end
function S.copy(v, memo)
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
        t[k] = S.copy(x, memo)
    end
    return t
end
local fields = {
    "materials",
    "regions",
    "boundaries",
    "frontiers",
    "domain",
    "pairs",
    "answer_pairing",
    "point_constraints",
    "adaptive_mesh",
    "topology",
    "native_topology",
    "circuits",
    "mesh_stats",
    "xs",
    "ys",
    "sliding",
    "machine_nodes",
    "frame",
    "component_basis",
}
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
--- Capture a converged static initial condition or an unchanged time step.
-- time is seconds; speed is rad/s. A template must be the exact sliding capture
-- used for the initial position. Caller-owned tables are copied recursively.
local function capture(m, o)
    o = o or {}
    local allowed = { time = true, initial = true, template = true, speed = true }
    for k in pairs(o) do
        check(allowed[k], "unknown capture option " .. tostring(k))
    end
    local M, C = require("luafemm"), require("luafemm-cache")
    check(
        o.initial == nil or o.initial == "equilibrium" or o.initial == "zero",
        "unknown initial condition"
    )
    if o.initial == "zero" then
        check(not m.transient and not m.stats, "zero initialisation requires an unsolved model")
        require("luafemm-frontiers").compile(m)
        if not m.nodes then
            M.mesh(m)
        end
        for _, v in pairs(m.fixed_values) do
            check(v == 0, "zero initialisation requires homogeneous Dirichlet data")
        end
        for i in ipairs(m.A) do
            m.A[i] = 0
        end
        m.initial_time = o.time or 0
        m.stats = {
            newton = 0,
            cg = 0,
            residual = 0,
            seconds = 0,
            assembly_seconds = 0,
            nodes = #m.nodes,
            triangles = #m.triangles,
        }
        M.update_fields(m)
        if m.sliding then
            m.sliding_contract = C.signature(m)
        end
    else
        M.solve(m)
    end
    check(m.solution_signature == C.signature(m), "stale solution")
    local d = {
        model = {},
        nodes = {},
        elements = {},
        boundary = S.copy(m.exterior_nodes),
        A = S.copy(m.A),
        time = m.transient and m.transient.time or o.time or m.initial_time or 0,
        transient = S.copy(m.transient),
        motion = S.copy(m.time_motion),
        sources = S.copy(m.time_sources),
        stats = S.copy(m.stats),
    }
    check(finite(d.time) and d.time >= 0, "time must be finite and nonnegative")
    check(
        not m.transient or o.time == nil or o.time == d.time,
        "cannot change a solved step's time"
    )
    d.model.options = S.copy(m.options)
    d.model.options.log, d.model.options.cache = nil, { mode = "off" }
    for _, k in ipairs(fields) do
        d.model[k] = S.copy(m[k])
    end
    for i, p in ipairs(m.nodes) do
        d.nodes[i] = { p[1], p[2] }
    end
    for i, e in ipairs(m.triangles) do
        d.elements[i] = S.copy(e.ids)
    end
    if not d.sources then
        d.sources =
            { currents = {}, boundaries = S.copy(m.boundaries), frontiers = S.copy(m.frontiers) }
        for i, r in ipairs(m.regions) do
            d.sources.currents[i] = r.current
        end
    end
    if o.template then
        check(not m.transient, "supply a motion template only for the initial state")
        check(
            m.sliding and o.template.provenance == m.sliding.provenance,
            "unrelated motion template"
        )
        check(
            not o.template.sector and o.template.options.method == "p1",
            "time motion requires a full-circle P1 template"
        )
        d.motion = {
            template = S.copy(o.template),
            angle = m.sliding.angle,
            offset_x = m.sliding.offset_x,
            offset_y = m.sliding.offset_y,
            speed = o.speed or 0,
        }
        d.motion.template.model.options.log = nil
    end
    check(o.speed == nil or d.motion, "speed requires a motion template")
    if o.speed ~= nil then
        check(not m.transient, "set initial speed only on an initial state")
        check(finite(o.speed), "invalid angular speed")
    end
    d.linkage = {}
    for name in pairs(m.circuits or {}) do
        local r = require("luafemm-static").result(m, name)
        if r.linkage then
            d.linkage[name] = r.linkage
        end
    end
    -- Energy is used for the discrete dissipation diagnostic, not as an
    -- initial guess or a substitute for the previous nodal potential.
    d.energy = require("luafemm-post").measure(
        m,
        { name = "_time", quantity = "integrals", register = false }
    ).energy
    local state = { format = "luafemm-time-1", version = M.version, data = d }
    state.signature = C.encode(d)
    return state
end
--- Capture under a bounded initial-equilibrium solve context.
-- A source-free isolated conductor carries no current in DC equilibrium.
-- Permit that solve and its measurements only after checking the time contract.
-- The guard is not physical data and must not survive either success or error.
function S.capture(m, o)
    require("luafemm-time").eligible(m)
    local previous = m.time_initialising
    m.time_initialising = not m.transient
    local ok, result = pcall(capture, m, o)
    m.time_initialising = previous
    if not ok then
        error(result, 0)
    end
    return result
end
--- Validate an owned state before any derived solve can observe its data.
function S.validate(state)
    local C = require("luafemm-cache")
    check(type(state) == "table" and state.format == "luafemm-time-1", "invalid state format")
    check(state.version == require("luafemm").version, "incompatible state version")
    check(
        type(state.data) == "table" and state.signature == C.encode(state.data),
        "state changed after capture"
    )
    local d = state.data
    check(finite(d.time) and d.time >= 0 and type(d.model) == "table", "invalid state metadata")
    check(#d.nodes >= 3 and #d.nodes <= 1000000 and #d.A == #d.nodes, "invalid nodal arrays")
    check(#d.elements > 0 and #d.elements <= 2000000, "invalid element array")
    for i, p in ipairs(d.nodes) do
        check(finite(p[1]) and finite(p[2]) and finite(d.A[i]), "invalid nodal value")
    end
    for _, e in ipairs(d.elements) do
        check(#e == 3, "invalid element")
        for _, id in ipairs(e) do
            check(finite(id) and id % 1 == 0 and d.nodes[id], "invalid node index")
        end
    end
    return d
end
--- Rebuild element/DOF metadata from the state's ordered mesh.
-- Physical declarations were prepared at capture; reapplying material mixing
-- here would incorrectly apply lamination fill for a second time.
function S.model(state)
    local d = S.validate(state)
    local m = require("luafemm").new(S.copy(d.model.options))
    for _, k in ipairs(fields) do
        m[k] = S.copy(d.model[k])
    end
    if m.topology then
        require("luafemm-static").excite(m, m.topology)
    end
    require("luafemm").prepare_mesh(m, S.copy(d.nodes), S.copy(d.elements), S.copy(d.boundary))
    m.A = S.copy(d.A)
    m.transient, m.time_motion, m.time_sources =
        S.copy(d.transient), S.copy(d.motion), S.copy(d.sources)
    m.time_mesh = { nodes = S.copy(d.nodes), elements = S.copy(d.elements) }
    return m
end
--- Restore a solved snapshot for viewing or measurements, without taking a step.
function S.restore(state)
    local d = S.validate(state)
    local m = S.model(state)
    m.stats = S.copy(d.stats)
    if not m.transient then
        m.initial_time = d.time
    end
    require("luafemm").update_fields(m)
    m.transient_contract = require("luafemm-cache").signature(m)
    if m.sliding then
        m.sliding_contract = m.transient_contract
    end
    return m
end
function S.write(state, file)
    S.validate(state)
    require("luafemm-cache").write_data(file, state)
end
function S.read(file)
    local state = require("luafemm-cache").read_data(file)
    S.validate(state)
    local m = S.model(state)
    require("luafemm-time").eligible(m)
    return state
end
return S
