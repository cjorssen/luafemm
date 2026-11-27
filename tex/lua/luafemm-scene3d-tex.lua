-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Process-local named scene registry and generic TikZ serialization.
-- Labels/styles are explicitly supplied TeX tokens; coordinates are finite numbers.
-- Scene data are independent copies and can be released after all desired views.
local T = {}
local S = require("luafemm-scene3d")
local scenes = {}
local function name(s)
    S.check(type(s) == "string" and s:match("^[%w][%w:_%-]*$"), "invalid scene name")
    return s
end
local function get(s)
    return assert(scenes[name(s)], "luafemm 3D view: unknown scene")
end
function T.store(id, scene)
    scenes[name(id)] = S.copy(scene)
end
function T.clear(id)
    get(id)
    scenes[id] = nil
end
function T.plane(id, o)
    S.add_plane(get(id), o)
end
local function coordinate(p)
    for k = 1, 2 do
        S.check(p[k] == p[k] and math.abs(p[k]) < 400, "annotation exceeds TeX coordinate range")
    end
    return string.format("(%.8fcm,%.8fcm)", p[1], p[2])
end
function T.picture(id, o)
    local layout = S.layout(get(id), o)
    local b = layout.bounds
    local clipping = "\\clip "
        .. coordinate({ b[1], b[2] })
        .. " rectangle "
        .. coordinate({ b[3], b[4] })
        .. ";"
    local text = {}
    local active
    for _, p in ipairs(layout.primitives) do
        local clipped = p.clip == true
        if active == nil or active ~= clipped then
            if active ~= nil then
                text[#text + 1] = "\\endscope"
            end
            text[#text + 1] = "\\scope"
            if clipped then
                text[#text + 1] = clipping
            end
            active = clipped
        end
        if p.kind == "label" then
            text[#text + 1] = "\\node["
                .. p.style
                .. "] at "
                .. coordinate(p.position)
                .. " {"
                .. p.text
                .. "};"
        else
            local points = {}
            for _, q in ipairs(p.points) do
                points[#points + 1] = coordinate(q)
            end
            local body = table.concat(points, "--")
            if p.controls then
                body = points[1]
                    .. " .. controls "
                    .. coordinate(p.controls[1])
                    .. " and "
                    .. coordinate(p.controls[2])
                    .. " .. "
                    .. points[#points]
            end
            text[#text + 1] = "\\draw["
                .. p.style
                .. "] "
                .. body
                .. (p.closed and "--cycle" or "")
                .. ";"
        end
    end
    if active ~= nil then
        text[#text + 1] = "\\endscope"
    end
    -- tex.sprint does not turn embedded newline characters into input lines.
    return table.concat(text, " ")
end
return T
