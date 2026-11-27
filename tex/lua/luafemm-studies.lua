-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Ordered angle sweeps and unit-current response matrices.
-- Studies own numeric results, not the transient meshes. Registration is atomic:
-- a failed solve cannot leave a plausible truncated curve under the requested name.
local W = {}
local results = {}
local function check(v, s)
    assert(v, "luafemm study: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
local function begin(o, allowed)
    check(type(o.name) == "string" and o.name:match("%S"), "name is required")
    check(not results[o.name], "duplicate study " .. o.name)
    for k in pairs(o) do
        check(allowed[k], "unknown option " .. tostring(k))
    end
end
local function cache(o, i)
    local prefix = o.cache_prefix or ""
    check(type(prefix) == "string", "invalid cache prefix")
    local mode = o.cache_mode or "auto"
    check(
        mode == "off" or mode == "auto" or mode == "refresh" or mode == "frozen",
        "invalid cache mode"
    )
    check(mode ~= "frozen" or prefix ~= "", "frozen studies require a cache prefix")
    return prefix == "" and { mode = "off" }
        or { mode = mode, file = prefix .. "-" .. i .. "-cache" }
end
--- Solve independent positions. Angles are degrees, in the caller's order.
-- Sources and frequency are inherited from the captured template.
function W.angles(template, o)
    begin(o, {
        name = true,
        angles = true,
        quantity = true,
        circuit = true,
        cache_prefix = true,
        cache_mode = true,
    })
    check(
        type(o.angles) == "table" and #o.angles > 0 and #o.angles <= 721,
        "provide 1 to 721 angles"
    )
    local quantity = o.quantity or "force"
    check(
        quantity == "force" or quantity == "circuit" or quantity == "integrals",
        "unknown quantity"
    )
    for _, a in ipairs(o.angles) do
        check(finite(a), "angles must be finite")
    end
    local out = { kind = "angles", rows = {} }
    for i, a in ipairs(o.angles) do
        local m = require("luafemm-sliding").new(template, { angle = a, cache = cache(o, i) })
        local r = require("luafemm-post").measure(m, {
            name = "_study",
            register = false,
            quantity = quantity,
            method = quantity == "force" and "band" or "weighted",
            circuit = o.circuit,
        })
        local row = { angle = a, index = i }
        for k, v in pairs(r) do
            if finite(v) then
                row[k] = v
            end
        end
        out.rows[i] = row
    end
    results[o.name] = out
    return out
end
--- Differential linkage (H) and, in AC, terminal impedance (ohm).
-- Column j is a unit-current, zero-phase incremental solve about the SAME DC
-- bias; every other independent source is zero. No reciprocity is imposed.
function W.matrix(bias, o)
    begin(
        o,
        { name = true, circuits = true, frequency = true, cache_prefix = true, cache_mode = true }
    )
    check(
        type(o.circuits) == "table" and #o.circuits > 0 and #o.circuits <= 32,
        "provide 1 to 32 distinct circuits"
    )
    local seen = {}
    for _, name in ipairs(o.circuits) do
        check(
            type(name) == "string" and bias.circuits and bias.circuits[name],
            "unknown circuit " .. tostring(name)
        )
        check(not seen[name], "duplicate circuit " .. name)
        seen[name] = true
    end
    local frequency = o.frequency or 0
    check(finite(frequency) and frequency >= 0, "frequency must be finite and nonnegative")
    local out = { kind = "matrix", circuits = {}, cells = {}, frequency = frequency }
    for i, name in ipairs(o.circuits) do
        out.circuits[i] = name
        out.cells[name] = {}
    end
    for j, source in ipairs(o.circuits) do
        local phases = {}
        for name in pairs(bias.circuits) do
            phases[name] = 0
        end
        local m = require("luafemm-operating").new(bias, {
            mode = "incremental",
            frequency = frequency,
            currents = { [source] = 1 },
            phases = phases,
            cache = cache(o, j),
        })
        require("luafemm").solve(m)
        for _, target in ipairs(o.circuits) do
            local r = require("luafemm-static").result(m, target)
            local cell = {
                linkage_real = r.linkage_real or r.linkage,
                linkage_imaginary = r.linkage_imaginary or 0,
            }
            check(finite(cell.linkage_real), "matrix requires closed signed-turn windings")
            if frequency > 0 then
                cell.impedance_real, cell.impedance_imaginary = r.voltage_real, r.voltage_imaginary
            end
            out.cells[target][source] = cell
        end
    end
    results[o.name] = out
    return out
end
function W.get(name)
    return assert(results[name], "luafemm study: unknown study " .. tostring(name))
end
function W.value(name, index, component)
    local r = W.get(name)
    check(r.kind == "angles", "angle study required")
    local v = r.rows[index] and r.rows[index][component]
    check(finite(v), "unavailable study component " .. tostring(component))
    return v
end
function W.entry(name, row, column, component)
    local r = W.get(name)
    check(r.kind == "matrix", "response matrix required")
    local v = r.cells[row] and r.cells[row][column] and r.cells[row][column][component]
    check(finite(v), "unavailable matrix entry")
    return v
end
function W.coordinates(name, component)
    local r = W.get(name)
    check(r.kind == "angles", "angle study required")
    local out = {}
    for i, row in ipairs(r.rows) do
        out[i] = string.format("(%.17g,%.17g)", row.angle, W.value(name, i, component))
    end
    return table.concat(out, " ")
end
return W
