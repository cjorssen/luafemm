-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Passive FEMM 4/5 electrostatic format-1 documents and raw solution records.
-- No imported text is evaluated. Indices stay native (zero-based geometry,
-- one-based property references) until the physical adapter validates them.
-- @module luafemm-fee
local I = {}
local E = require("luafemm-electrostatic")
local function check(v, s)
    assert(v, "luafemm electric file: " .. s)
end
local sections = {
    pointprops = { "point", "pointprops" },
    bdryprops = { "bdry", "boundaries" },
    blockprops = { "block", "materials" },
    conductorprops = { "conductor", "conductors" },
    numpoints = { 5, "points" },
    numsegments = { 7, "segments" },
    numarcsegments = { 8, "arcs" },
    numholes = { 2, "holes" },
    numblocklabels = { 6, "labels" },
}
local fields = {
    point = { pointname = true, vp = true, qp = true },
    bdry = { bdryname = true, bdrytype = true, vs = true, qs = true, c0 = true, c1 = true },
    block = { blockname = true, ex = true, ey = true, qv = true },
    conductor = { conductorname = true, vc = true, qc = true, conductortype = true },
}
local headers = {
    format = true,
    precision = true,
    minangle = true,
    depth = true,
    lengthunits = true,
    problemtype = true,
    coordinates = true,
    comment = true,
    extzo = true,
    extro = true,
    extri = true,
    dosmartmesh = true,
    smartmesh = true,
}
local function scalar(s)
    s = s:match("^%s*(.-)%s*$")
    if s:sub(1, 1) == '"' then
        check(s:sub(-1) == '"', "unterminated quoted string")
        return s:sub(2, -2):gsub("\\n", "\n"):gsub('\\"', '"')
    end
    return tonumber(s) or s
end
local function row(s, lo, hi)
    local r = {}
    for v in s:gmatch("%S+") do
        r[#r + 1] = E.number(v, "file row")
    end
    check(#r >= lo and #r <= (hi or lo), "wrong column count")
    return r
end
local function count(s)
    return E.integer(s, "record count", 0, 1000000)
end
function I.parse(text)
    check(type(text) == "string" and #text <= 128 * 1024 * 1024, "file exceeds 128 MiB")
    local lines = {}
    for s in (text .. "\n"):gmatch("(.-)\n") do
        s = s:gsub("\r$", ""):match("^%s*(.-)%s*$")
        if s ~= "" then
            lines[#lines + 1] = s
        end
    end
    local d = { header = {} }
    for _, v in pairs(sections) do
        d[v[2]] = {}
    end
    local i, seen = 1, {}
    while i <= #lines do
        local name, value = lines[i]:match("^%[([^%]]+)%]%s*=%s*(.-)%s*$")
        if lines[i]:lower() == "[solution]" then
            check(not d.solution, "duplicate solution")
            local s = { nodes = {}, triangles = {}, conductors = {} }
            d.solution = s
            i = i + 1
            local n = count(lines[i])
            i = i + 1
            for j = 1, n do
                s.nodes[j] = row(assert(lines[i], "truncated solution"), 4)
                i = i + 1
            end
            n = count(lines[i])
            i = i + 1
            for j = 1, n do
                s.triangles[j] = row(assert(lines[i], "truncated solution"), 4)
                i = i + 1
            end
            n = count(lines[i])
            i = i + 1
            check(n == #d.conductors, "solution conductor count mismatch")
            for j = 1, n do
                s.conductors[j] = row(assert(lines[i], "truncated solution"), 2)
                i = i + 1
            end
            check(i > #lines, "trailing solution data are unsupported")
        else
            check(name, "expected header or section at line " .. i)
            name = name:lower()
            check(not seen[name], "duplicate section " .. name)
            seen[name] = true
            local spec = sections[name]
            i = i + 1
            if spec then
                local n = count(value)
                for j = 1, n do
                    local r
                    if type(spec[1]) == "number" then
                        r = row(
                            assert(lines[i], "truncated geometry"),
                            spec[1],
                            name == "numarcsegments" and 9
                                or (name == "numpoints" and 6 or (name == "numholes" and 3 or nil))
                        )
                        i = i + 1
                    else
                        check(
                            lines[i] and lines[i]:lower() == "<begin" .. spec[1] .. ">",
                            "missing property begin"
                        )
                        i = i + 1
                        r = {}
                        while lines[i] and lines[i]:lower() ~= "<end" .. spec[1] .. ">" do
                            local k, v = lines[i]:match("^<([^>]+)>%s*=%s*(.-)%s*$")
                            check(k, "invalid property record")
                            k = k:lower()
                            check(
                                fields[spec[1]][k] and r[k] == nil,
                                "unknown or duplicate property " .. k
                            )
                            r[k] = scalar(v)
                            i = i + 1
                        end
                        check(lines[i], "missing property end")
                        i = i + 1
                    end
                    d[spec[2]][j] = r
                end
            else
                check(headers[name], "unsupported header " .. name)
                d.header[name] = scalar(value)
            end
        end
    end
    check(d.header.format == 1, "only electrostatic format 1 is supported")
    check(
        d.header.problemtype == "planar" or d.header.problemtype == "axisymmetric",
        "unsupported geometry"
    )
    check(
        d.header.coordinates == nil
            or d.header.coordinates == "cartesian"
            or d.header.coordinates == "polar",
        "unknown coordinate display"
    )
    check(
        #d.points >= 3 and #d.materials > 0 and #d.labels > 0,
        "missing geometry or dielectric labels"
    )
    return d
end
function I.read(path)
    check(type(path) == "string" and not path:find("%c"), "invalid input filename")
    local f = assert(io.open(path, "rb"))
    local text = f:read(128 * 1024 * 1024 + 1)
    f:close()
    local d = I.parse(text)
    d.source = path
    return d
end
local function scalar_out(v)
    if type(v) == "number" then
        E.number(v, "output number")
        return string.format("%.17g", v)
    end
    check(type(v) == "string" and not v:find('[\r\n"\\]'), "unrepresentable native string")
    return '"' .. v .. '"'
end
function I.serialize(d)
    local out = {}
    local function put(s)
        out[#out + 1] = s
    end
    for _, k in ipairs({
        "format",
        "precision",
        "minangle",
        "depth",
        "lengthunits",
        "problemtype",
        "coordinates",
        "comment",
    }) do
        local v = d.header[k]
        if v ~= nil then
            put("[" .. k .. "] = " .. (k == "comment" and scalar_out(v) or tostring(v)))
        end
    end
    for _, k in ipairs({
        "pointprops",
        "bdryprops",
        "blockprops",
        "conductorprops",
        "numpoints",
        "numsegments",
        "numarcsegments",
        "numholes",
        "numblocklabels",
    }) do
        local spec = sections[k]
        local list = d[spec[2]]
        put("[" .. k .. "] = " .. #list)
        for _, r in ipairs(list) do
            if type(spec[1]) == "number" then
                local a = {}
                for _, v in ipairs(r) do
                    a[#a + 1] = scalar_out(v)
                end
                put(table.concat(a, "\t"))
            else
                put("<Begin" .. spec[1]:gsub("^%l", string.upper) .. ">")
                local keys = {}
                for name in pairs(r) do
                    keys[#keys + 1] = name
                end
                table.sort(keys)
                for _, name in ipairs(keys) do
                    put("<" .. name .. "> = " .. scalar_out(r[name]))
                end
                put("<End" .. spec[1]:gsub("^%l", string.upper) .. ">")
            end
        end
    end
    if d.solution then
        put("[Solution]")
        for _, list in ipairs({ d.solution.nodes, d.solution.triangles, d.solution.conductors }) do
            put(tostring(#list))
            for _, r in ipairs(list) do
                local a = {}
                for _, v in ipairs(r) do
                    a[#a + 1] = scalar_out(v)
                end
                put(table.concat(a, "\t"))
            end
        end
    end
    return table.concat(out, "\n") .. "\n"
end
return I
