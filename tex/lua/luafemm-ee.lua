-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Passive IEC views of physical declarations and named result snapshots.
-- Output contains fixed TeX templates, validated ASCII identifiers and numbers.
-- No label is parsed as physics and no view triggers a solve.
local V = {}
local function check(v, s)
    assert(v, "luafemm IEC: " .. s)
end
local function id(s)
    check(type(s) == "string" and s:match("^[%w_-]+$"), "invalid view identifier")
    return (s:gsub("_", "\\_"))
end
local function model(o)
    local M = require("luafemm")
    local r = o.result ~= "" and require("luafemm-post").get(o.result) or nil
    local m = o.model ~= "" and M.tex.models[o.model] or r and r.model or M.current
    check(m, "unknown view model")
    check(o.model == "" or M.tex.models[o.model] == m, "unknown view model")
    check(not r or r.model == m, "result belongs to another model")
    return m, r
end
local function resolve(o)
    local m, r = model(o)
    check(o.circuit == "" or o.element == "", "select a circuit or an element")
    local c = o.circuit ~= "" and m.circuits and m.circuits[o.circuit]
    local e = o.element ~= "" and m.network and m.network.elements[o.element]
    check(c or e, "bind a known circuit or element")
    check(not (c and e), "select a circuit or an element")
    return m, r, c, e
end
function V.label(o)
    local m, r, c, e = resolve(o)
    local output = {}
    for key in o.show:gmatch("[^, ]+") do
        if key == "name" then
            output[#output + 1] = "{\\rm " .. id(c and o.circuit or o.element) .. "}"
        elseif key ~= "none" then
            check(
                key == "current" or key == "voltage" or key == "power" or key == "linkage",
                "unsupported displayed quantity"
            )
            local value, component
            if o.values == "solution" then
                check(r, "solution values require a named result")
                if c then
                    check(
                        r.quantity == "circuit" and r.circuit == o.circuit,
                        "result is not for the selected circuit"
                    )
                else
                    check(
                        r.quantity == "network" and r.element == o.element,
                        "result is not for the selected element"
                    )
                end
                component = m.options.frequency > 0
                        and key ~= "power"
                        and key .. "_" .. o.representation
                    or key
                value = r[component]
            else
                check(o.values == "input", "unknown value source")
                if c and c.drive == key then
                    value = c[key]
                elseif e and e.kind == key .. " source" then
                    value = e.value
                end
                component = key
            end
            check(type(value) == "number", "quantity is not available: " .. key)
            if o.reference == "reverse" and key ~= "power" then
                if component:match("_phase$") then
                    value = (value + 360) % 360 - 180
                elseif not component:match("_amplitude$") and not component:match("_rms$") then
                    value = -value
                end
            end
            local symbol = ({ current = "I", voltage = "U", power = "P", linkage = "\\lambda" })[key]
            local unit = require("luafemm-post").unit(component)
            output[#output + 1] = symbol
                .. "=\\pgfmathprintnumber{"
                .. string.format("%.12g", value)
                .. "}\\,{\\rm "
                .. unit
                .. "}"
        end
    end
    return "$" .. table.concat(output, "\\,;\\quad ") .. "$"
end
--- Canonical port loop or explicit branch-list view of a whole graph.
-- Repeated terminal identifiers in a branch-list denote the same electrical node.
-- Layout is intentionally simple and can be replaced with native hand placement.
function V.picture(o)
    local m = model(o)
    local out = {}
    if o.circuit ~= "" then
        local c = m.circuits and m.circuits[o.circuit]
        check(c, "unknown circuit")
        local source = c.drive == "current" and "current source" or "voltage source"
        check(
            c.drive == "current" or c.drive == "voltage",
            "automatic port loop needs current/voltage drive; use a bound winding view"
        )
        out[1] = "\\draw[femm ee wire] (0,0) to["
            .. source
            .. "] (0,2) -- (3,2) to[inductor={femm/ee={}}] (3,0) -- (0,0);"
    else
        check(m.network, "no network to exhibit")
        check(
            o.values == "input",
            "whole-network view uses input names; bind measured elements for solution annotations"
        )
        local ids = {}
        for name in pairs(m.network.elements) do
            ids[#ids + 1] = name
        end
        table.sort(ids)
        for j, name in ipairs(ids) do
            local e = m.network.elements[name]
            local symbol = e.kind == "winding" and "inductor" or e.kind
            out[#out + 1] = string.format(
                "\\draw[femm ee wire] (0,%g) node[left]{%s} to[%s={info={${\\rm %s}$}}] (3,%g) node[right]{%s};",
                -(j - 1) * 1.2,
                id(e.p),
                symbol,
                id(name),
                -(j - 1) * 1.2,
                id(e.n)
            )
        end
    end
    return table.concat(out, " ")
end
return V
