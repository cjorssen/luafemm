-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Display-only arrow supports shared by canvas paths and projected scenes.
-- All distances use the caller's drawing units. Original curves and arclength
-- positions remain unchanged. A bounded local Hermite fit is used only for the
-- small carrier; adaptive rigid tips are reduced when the curve turns sharply.
local A = {}
local function check(ok, s)
    assert(ok, "luafemm arrows: " .. s)
end
local function distance(a, b)
    return math.sqrt((a[1] - b[1]) ^ 2 + (a[2] - b[2]) ^ 2)
end
function A.path(points)
    local p = { points = { points[1] }, s = { 0 }, length = 0 }
    for i = 2, #points do
        local d = distance(points[i], p.points[#p.points])
        if d > 1e-12 then
            p.length = p.length + d
            p.points[#p.points + 1] = points[i]
            p.s[#p.s + 1] = p.length
        end
    end
    return p
end
function A.at(p, s)
    s = math.max(0, math.min(p.length, s))
    local lo, hi = 2, #p.points
    while lo < hi do
        local mid = math.floor((lo + hi) / 2)
        if p.s[mid] < s then
            lo = mid + 1
        else
            hi = mid
        end
    end
    if not p.points[lo] then
        return
    end
    local a, b = p.points[lo - 1], p.points[lo]
    local d = p.s[lo] - p.s[lo - 1]
    local t = (s - p.s[lo - 1]) / d
    return {
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
        (a[3] or 0) + t * ((b[3] or 0) - (a[3] or 0)),
    }, {
        (b[1] - a[1]) / d,
        (b[2] - a[2]) / d,
    }
end
local function tangent(p, s, h)
    local a, b = A.at(p, s - h), A.at(p, s + h)
    local d = distance(a, b)
    if d < 1e-12 then
        return
    end
    return { (b[1] - a[1]) / d, (b[2] - a[2]) / d }
end
--- Return carrier and size factor, or nil when a stable tip cannot fit.
-- available is the unoccluded backward arclength, supplied by the renderer.
function A.mark(p, s, o, available)
    local length, tol, minimum = o.length, o.tolerance or 0.002, o.minimum or 0.25
    local mode = o.mode or "bend"
    check(length and length > 0 and length < math.huge, "positive finite arrow length required")
    check(
        tol > 0 and tol < math.huge and minimum > 0 and minimum <= 1,
        "invalid arrow fit controls"
    )
    check(mode == "adaptive" or mode == "bend" or mode == "tangent", "unknown arrow mode")
    if p.length == 0 then
        return
    end
    local factor = math.min(1, s / (2 * length), (available or s) / (2 * length))
    while factor >= minimum do
        local span = 2 * length * factor
        local a, b = A.at(p, s - span), A.at(p, s)
        local ta, tb = tangent(p, s - span, span * 0.08), tangent(p, s, span * 0.08)
        if ta and tb then
            local c = { a[1] + ta[1] * span / 3, a[2] + ta[2] * span / 3 }
            local d = { b[1] - tb[1] * span / 3, b[2] - tb[2] * span / 3 }
            local error = 0
            for i = 1, 15 do
                local t = i / 16
                local u = 1 - t
                local q = A.at(p, s - span + t * span)
                local fit = {
                    u ^ 3 * a[1] + 3 * u * u * t * c[1] + 3 * u * t * t * d[1] + t ^ 3 * b[1],
                    u ^ 3 * a[2] + 3 * u * u * t * c[2] + 3 * u * t * t * d[2] + t ^ 3 * b[2],
                }
                error = math.max(error, distance(fit, q))
                if mode == "adaptive" and t >= 0.5 then
                    local offset = (1 - t) * span
                    error = math.max(
                        error,
                        distance(q, { b[1] - tb[1] * offset, b[2] - tb[2] * offset })
                    )
                end
            end
            if mode == "tangent" or error <= tol then
                return {
                    points = { a, b },
                    controls = { c, d },
                    factor = factor,
                    start = s - span,
                    finish = s,
                    error = error,
                }
            end
        end
        factor = factor * 0.75
    end
end
--- Generic PGF postaction adapter. Soft paths are already in canvas points;
-- matrices are reset by the TeX caller, so scaling is applied exactly once.
function A.tex(soft, o)
    local P = require("luafemm-path")
    for _, key in ipairs({ "length", "tolerance", "minimum", "position", "spacing" }) do
        o[key] = tonumber(o[key])
        check(
            o[key] and o[key] == o[key] and math.abs(o[key]) < math.huge,
            "finite " .. key .. " required"
        )
    end
    check(
        o.position >= 0 and o.position <= 1 and o.spacing >= 0,
        "invalid mark position or spacing"
    )
    local result, count = {}, 0
    -- Split independent moves before invoking the established single-contour
    -- reader. Isolated moves emitted after closed PGF shapes are harmless.
    local starts = {}
    for pos in soft:gmatch("()\\pgfsyssoftpath@movetotoken") do
        starts[#starts + 1] = pos
    end
    if #starts == 0 then
        starts = { 1 }
    elseif starts[1] > 1 then
        table.insert(starts, 1, 1)
    end
    starts[#starts + 1] = #soft + 1
    for i = 1, #starts - 1 do
        local part = soft:sub(starts[i], starts[i + 1] - 1)
        if
            part:find("linetotoken", 1, true)
            or part:find("curvetotoken", 1, true)
            or part:find("rectsizetoken", 1, true)
        then
            local points = P.read(part, { 1, 0, 0, 1, 0, 0 }, o.tolerance / 4, true)
            local p = A.path(points)
            local function emit(s)
                count = count + 1
                check(count <= 20000, "arrow budget exceeded")
                local m = A.mark(p, s, o)
                if m then
                    local a, b = m.points[1], m.points[2]
                    local c, d = m.controls[1], m.controls[2]
                    result[#result + 1] = string.format(
                        "\\csname luafemm@curvearrow\\endcsname{%.10f}{%.10f}{%.10f}{%.10f}{%.10f}{%.10f}{%.10f}{%.10f}{%.10f}",
                        m.factor,
                        a[1],
                        a[2],
                        c[1],
                        c[2],
                        d[1],
                        d[2],
                        b[1],
                        b[2]
                    )
                end
            end
            if o.spacing == 0 then
                emit(o.position * p.length)
            else
                for s = o.position * p.length, 0.95 * p.length, o.spacing do
                    emit(s)
                end
            end
        end
    end
    return table.concat(result, " ")
end
return A
