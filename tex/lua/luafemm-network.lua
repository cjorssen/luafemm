-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Lumped electrical graphs coupled to magnetic winding ports.
-- All voltages use V(p)-V(n), currents flow p to n, and complex values are
-- peak phasors with exp(+i omega t). Drawing coordinates never enter this module.
-- Modified nodal equations retain an explicit current for every branch. This
-- makes DC open capacitors and short inductors exact, without artificial shunts.
local E = {}
local function check(v, s)
    assert(v, "luafemm network: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
local function name(v)
    check(type(v) == "string" and v:match("^[%w_-]+$"), "invalid identifier")
    return v
end
local function keys(t)
    local out = {}
    for k in pairs(t) do
        out[#out + 1] = k
    end
    table.sort(out)
    return out
end
local kinds = {
    resistor = true,
    inductor = true,
    capacitor = true,
    ["current source"] = true,
    ["voltage source"] = true,
    winding = true,
}
--- Declare one graph per magnetic problem, before any mesh exists.
function E.declare(m, o)
    check(m and not m.nodes, "declare the network before meshing")
    check(not m.import_document, "declare networks on a native magnetic problem")
    check(not m.network, "duplicate network")
    for k in pairs(o) do
        check(k == "name" or k == "ground", "unknown network option " .. k)
    end
    m.network = { name = name(o.name), ground = name(o.ground or "0"), elements = {} }
end
--- Register an explicitly oriented element; values are SI, phase is in degrees.
function E.element(m, o)
    check(m and m.network and not m.nodes, "declare a network and its elements before meshing")
    for k in pairs(o) do
        check(
            k == "name"
                or k == "kind"
                or k == "p"
                or k == "n"
                or k == "value"
                or k == "phase"
                or k == "circuit",
            "unknown element option " .. k
        )
    end
    local e = {
        name = name(o.name),
        kind = o.kind,
        p = name(o.p),
        n = name(o.n),
        value = o.value or 0,
        phase = o.phase or 0,
        circuit = o.circuit or "",
    }
    check(kinds[e.kind], "unsupported element kind")
    check(e.p ~= e.n, "element terminals must differ")
    check(finite(e.value) and finite(e.phase), "nonfinite element value or phase")
    check(not m.network.elements[e.name], "duplicate element " .. e.name)
    if e.kind == "winding" then
        name(e.circuit)
        check(e.value == 0 and e.phase == 0, "winding values come from the field model")
    elseif e.kind == "resistor" or e.kind == "inductor" or e.kind == "capacitor" then
        check(
            e.value > 0 and e.phase == 0 and e.circuit == "",
            "passive elements need a positive value and no phase/circuit"
        )
    else
        check(e.circuit == "", "source cannot bind a winding")
    end
    m.network.elements[e.name] = e
end
--- Allocate graph unknowns after any pre-existing field/circuit unknowns.
-- bindings maps winding names to current callbacks, voltage constraints are
-- supplied by the magnetic solver; the returned graph contains no owned model.
function E.prepare(m, offset)
    if not m.network then
        return
    end
    local n, g = m.network, { nodes = {}, elements = {}, ports = {}, size = offset or 0 }
    local nodes, adjacency = {}, {}
    for _, id in ipairs(keys(n.elements)) do
        local e = n.elements[id]
        check(m.options.frequency > 0 or e.phase == 0, "DC sources require zero phase")
        for _, a in ipairs({ e.p, e.n }) do
            nodes[a] = true
            adjacency[a] = adjacency[a] or {}
        end
        adjacency[e.p][e.n], adjacency[e.n][e.p] = true, true
        if e.kind == "winding" then
            check(
                m.circuits and m.circuits[e.circuit] and m.circuits[e.circuit].drive == "network",
                "winding requires a circuit with drive=network: " .. e.circuit
            )
            check(not g.ports[e.circuit], "winding port connected more than once")
            g.ports[e.circuit] = e
        end
        g.elements[#g.elements + 1] = e
    end
    check(nodes[n.ground], "reference node is not connected")
    check(#g.elements <= 256 and #keys(nodes) <= 257, "network exceeds 256 unknowns")
    local seen, queue = { [n.ground] = true }, { n.ground }
    local head = 1
    while queue[head] do
        for b in pairs(adjacency[queue[head]]) do
            if not seen[b] then
                seen[b] = true
                queue[#queue + 1] = b
            end
        end
        head = head + 1
    end
    for _, a in ipairs(keys(nodes)) do
        check(seen[a], "disconnected node " .. a)
        if a ~= n.ground then
            g.size = g.size + 1
            g.nodes[a] = g.size
        end
    end
    g.currents = {}
    for _, e in ipairs(g.elements) do
        if e.kind ~= "winding" then
            g.size = g.size + 1
            g.currents[e.name] = g.size
        end
    end
    check(g.size <= 256, "network exceeds 256 unknowns")
    for id, c in pairs(m.circuits or {}) do
        check(c.drive ~= "network" or g.ports[id], "unconnected winding port " .. id)
    end
    return g
end
local function voltage(g, e, u)
    return (u[g.nodes[e.p]] or 0) - (u[g.nodes[e.n]] or 0)
end
E.voltage = voltage
--- Append KCL and component laws; row accepts residual real and imaginary parts.
-- port(name) returns the actual field-port current in the trial solution.
function E.equations(m, g, ur, ui, targets, port, row)
    local kr, ki = {}, {}
    for a in pairs(g.nodes) do
        kr[a], ki[a] = 0, 0
    end
    local w = 2 * math.pi * m.options.frequency
    for _, e in ipairs(g.elements) do
        local ir, ii
        if e.kind == "winding" then
            ir, ii = port(e.circuit)
        else
            ir, ii = ur[g.currents[e.name]], ui[g.currents[e.name]]
        end
        for _, pair in ipairs({ { e.p, 1 }, { e.n, -1 } }) do
            if g.nodes[pair[1]] then
                kr[pair[1]] = kr[pair[1]] + pair[2] * ir
                ki[pair[1]] = ki[pair[1]] + pair[2] * ii
            end
        end
        if e.kind ~= "winding" then
            local vr, vi = voltage(g, e, ur), voltage(g, e, ui)
            local r, i = vr, vi
            if e.kind == "resistor" then
                r, i = vr - e.value * ir, vi - e.value * ii
            elseif e.kind == "inductor" then
                r, i = vr + w * e.value * ii, vi - w * e.value * ir
            elseif e.kind == "capacitor" then
                r, i = ir + w * e.value * vi, ii - w * e.value * vr
            else
                if e.kind == "current source" then
                    r, i = ir, ii
                end
                if targets then
                    r = r - e.value * math.cos(math.rad(e.phase))
                    i = i - e.value * math.sin(math.rad(e.phase))
                end
            end
            row(r, i, "branch")
        end
    end
    for _, a in ipairs(keys(g.nodes)) do
        row(kr[a], ki[a], "kcl")
    end
end
--- Solve DC electrical constraints using the winding resistances already
-- obtained from the painted geometry. DC magnetic nonlinearity does not alter
-- the terminal resistance; linkage constraints are handled by the field solver.
function E.dc(m, resolved)
    if not m.network or m.options.frequency > 0 then
        return
    end
    local ports = {}
    for _, id in ipairs(keys(m.circuits or {})) do
        if m.circuits[id].drive == "network" then
            ports[#ports + 1] = id
        end
    end
    local g = E.prepare(m, #ports)
    local index = {}
    for j, id in ipairs(ports) do
        index[id] = j
    end
    local function equations(u, targets)
        local r, z = {}, {}
        for j = 1, g.size do
            z[j] = 0
        end
        for _, id in ipairs(ports) do
            local resistance = resolved[id].resistance
            check(resistance, "network winding requires conductivity and conductor length")
            r[#r + 1] = resistance * u[index[id]] - voltage(g, g.ports[id], u)
        end
        E.equations(m, g, u, z, targets, function(id)
            return u[index[id]], 0
        end, function(a)
            r[#r + 1] = a
        end)
        return r
    end
    local z = {}
    for j = 1, g.size do
        z[j] = 0
    end
    local b = equations(z, true)
    local ar, ai, bi = {}, {}, {}
    for k = 1, g.size do
        ar[k], ai[k], bi[k] = {}, {}, 0
        b[k] = -b[k]
    end
    for j = 1, g.size do
        local u = {}
        for k = 1, g.size do
            u[k] = k == j and 1 or 0
        end
        local r = equations(u, false)
        for k = 1, g.size do
            ar[k][j], ai[k][j] = r[k], 0
        end
    end
    local u = require("luafemm-ac-circuits").dense(ar, ai, b, bi)
    for _, id in ipairs(ports) do
        local c = resolved[id]
        c.current, c.voltage = u[index[id]], c.resistance * u[index[id]]
    end
    E.results(m, g, u, z, function(id)
        return u[index[id]], 0
    end)
end
--- Snapshot branch/node values and independent absolute KCL/branch residuals.
function E.results(m, g, ur, ui, port)
    local out = {
        nodes = {},
        elements = {},
        kcl_residual = 0,
        branch_residual = 0,
        power = 0,
        reactive_power = 0,
    }
    local decorate = require("luafemm-ac-circuits").decorate
    out.nodes[m.network.ground] = decorate({ voltage_real = 0, voltage_imaginary = 0 })
    if m.options.frequency == 0 then
        out.nodes[m.network.ground].voltage = 0
    end
    for a, j in pairs(g.nodes) do
        out.nodes[a] = decorate({ voltage_real = ur[j], voltage_imaginary = ui[j] })
        if m.options.frequency == 0 then
            out.nodes[a].voltage = ur[j]
        end
    end
    for _, e in ipairs(g.elements) do
        local ir, ii
        if e.kind == "winding" then
            ir, ii = port(e.circuit)
        else
            ir, ii = ur[g.currents[e.name]], ui[g.currents[e.name]]
        end
        local vr, vi = voltage(g, e, ur), voltage(g, e, ui)
        local factor = m.options.frequency > 0 and 0.5 or 1
        local r = decorate({
            current_real = ir,
            current_imaginary = ii,
            voltage_real = vr,
            voltage_imaginary = vi,
            power = factor * (vr * ir + vi * ii),
            reactive_power = factor * (vi * ir - vr * ii),
        })
        if m.options.frequency == 0 then
            r.current, r.voltage = ir, vr
        end
        out.elements[e.name] = r
        out.power = out.power + r.power
        out.reactive_power = out.reactive_power + r.reactive_power
    end
    E.equations(m, g, ur, ui, true, port, function(r, i, kind)
        local key = kind .. "_residual"
        out[key] = math.max(out[key], math.sqrt(r * r + i * i))
    end)
    m.network_result = out
    return out
end
return E
