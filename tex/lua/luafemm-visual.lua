-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Quantitative magnetic colour cells and vectors, separated from physics.
-- Raw fields use the actual P1 basis. Recovered vectors are area-weighted at
-- shared vertices within one material only, then barycentrically interpolated.
-- Recovery is a display operation, never the source for forces or loss integrals.
local V = {}
local M, H = require("luafemm"), require("luafemm-harmonic")
local function check(v, s)
    assert(v, "luafemm visual: " .. s)
end
local function number(v, name)
    local n = tonumber(v)
    check(n and n == n and math.abs(n) < math.huge, "finite " .. name .. " required")
    return n
end
local units = { B = "T", Bx = "T", By = "T", H = "A/m", Hx = "A/m", Hy = "A/m", A = "T m" }
function V.prepare(m, o)
    o = o or {}
    check(
        m.options.field_model == "open" and m.options.problem_type ~= "axial" and not m.sliding,
        "requires an ordinary planar or axisymmetric P1 field"
    )
    M.solve(m)
    local v = {
        model = m,
        component = o.component or "B",
        values = o.values or "raw",
        representation = o.representation or "snapshot",
        phase = number(o.phase or 0, "phase"),
    }
    check(units[v.component], "unknown scalar component")
    check(v.values == "raw" or v.values == "recovered", "unknown value policy")
    H.value(1, 1, v.representation, v.phase)
    check(v.representation ~= "phase", "phase maps are unavailable")
    check(
        m.options.frequency > 0
            or v.representation == "real"
            or (v.representation == "snapshot" and v.phase == 0),
        "static maps require real values or a zero-phase snapshot"
    )
    v.unit = units[v.component]
    v.index = require("luafemm-spatial").new(m.options.h * m.options.unit)
    v.recovery = {}
    for i, e in ipairs(m.triangles) do
        local x0, y0, x1, y1 = math.huge, math.huge, -math.huge, -math.huge
        local x, y = 0, 0
        for _, id in ipairs(e.ids) do
            local p = m.nodes[id]
            x0, y0, x1, y1 =
                math.min(x0, p[1]), math.min(y0, p[2]), math.max(x1, p[1]), math.max(y1, p[2])
            x, y = x + p[1] / 3, y + p[2] / 3
        end
        v.index.insert(i, x0, y0, x1, y1)
        if v.values == "recovered" then
            local f = V.raw(m, e, x, y)
            local material = e.material
            local group = v.recovery[material] or {}
            v.recovery[material] = group
            for _, id in ipairs(e.ids) do
                local n = group[id] or { weight = 0 }
                group[id] = n
                n.weight = n.weight + e.area
                for _, key in ipairs({ "Bx", "By", "Hx", "Hy" }) do
                    n[key] = n[key] or { 0, 0 }
                    for part = 1, 2 do
                        n[key][part] = n[key][part] + e.area * f[key][part]
                    end
                end
            end
        end
    end
    return v
end
function V.raw(m, e, x, y)
    if m.options.frequency > 0 then
        return H.field(m, e, x, y)
    end
    local bx, by, a, hx, hy
    if m.options.problem_type == "axisymmetric" then
        bx, by, a, hx, hy = require("luafemm-axisymmetric").field(m, e, x, y)
    else
        bx, by, hx, hy = e.bx, e.by, e.hx, e.hy
        local p = m.nodes[e.ids[1]]
        local b, c =
            e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2]), e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
        a = (1 - b - c) * m.A[e.ids[1]] + b * m.A[e.ids[2]] + c * m.A[e.ids[3]]
    end
    return { Bx = { bx, 0 }, By = { by, 0 }, Hx = { hx, 0 }, Hy = { hy, 0 }, A = { a, 0 } }
end
function V.field(v, i, x, y)
    local e, m = v.model.triangles[i], v.model
    local f = V.raw(m, e, x, y)
    if v.values == "recovered" then
        local p = m.nodes[e.ids[1]]
        local b, c =
            e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2]), e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
        local w = { 1 - b - c, b, c }
        for _, key in ipairs({ "Bx", "By", "Hx", "Hy" }) do
            f[key] = { 0, 0 }
            for j, id in ipairs(e.ids) do
                local n = v.recovery[e.material][id]
                for part = 1, 2 do
                    f[key][part] = f[key][part] + w[j] * n[key][part] / n.weight
                end
            end
        end
    end
    return f
end
function V.scalar(v, f)
    local key = v.component
    if key == "B" or key == "H" then
        local x, y = f[key .. "x"], f[key .. "y"]
        if v.representation == "amplitude" or v.representation == "rms" then
            local s = x[1] ^ 2 + x[2] ^ 2 + y[1] ^ 2 + y[2] ^ 2
            return math.sqrt(s / (v.representation == "rms" and 2 or 1))
        end
        return math.sqrt(
            H.value(x[1], x[2], v.representation, v.phase) ^ 2
                + H.value(y[1], y[2], v.representation, v.phase) ^ 2
        )
    end
    return H.value(f[key][1], f[key][2], v.representation, v.phase)
end
function V.locate(v, x, y)
    for _, i in ipairs(v.index.query(x, y, x, y)) do
        local e = v.model.triangles[i]
        local p = v.model.nodes[e.ids[1]]
        local b, c =
            e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2]), e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
        if b >= -1e-10 and c >= -1e-10 and b + c <= 1 + 1e-10 then
            return i
        end
    end
end
function V.cells(v, o)
    o = o or {}
    local m = v.model
    local cells, lo, hi = {}, math.huge, -math.huge
    for i, e in ipairs(m.triangles) do
        local x, y = 0, 0
        for _, id in ipairs(e.ids) do
            x = x + m.nodes[id][1] / 3
            y = y + m.nodes[id][2] / 3
        end
        local q = V.scalar(v, V.field(v, i, x, y))
        check(q == q and math.abs(q) < math.huge, "non-finite map sample")
        cells[i] = q
        lo, hi = math.min(lo, q), math.max(hi, q)
    end
    if o.minimum and o.minimum ~= "auto" then
        lo = number(o.minimum, "minimum")
    end
    if o.maximum and o.maximum ~= "auto" then
        hi = number(o.maximum, "maximum")
    end
    check(lo <= hi, "reversed colour range")
    -- Numerical roundoff in an affine patch must not create an effectively
    -- empty PGFPlots colour-bar range. Respect explicitly distinct limits.
    if
        lo == hi
        or (
            (not o.minimum or o.minimum == "auto")
            and (not o.maximum or o.maximum == "auto")
            and hi - lo < 1e-4 * math.max(math.abs(lo), math.abs(hi))
        )
    then
        local d = math.max(math.abs(lo) * 0.01, 1e-12)
        lo, hi = lo - d, hi + d
    end
    return cells, lo, hi
end
function V.vectors(v, o)
    o = o or {}
    check(v.component == "B" or v.component == "H", "vectors require component=B or H")
    check(
        v.representation == "snapshot"
            or v.representation == "real"
            or v.representation == "imaginary",
        "vectors require a signed instantaneous representation"
    )
    local m = v.model
    local span = math.max(m.options.xmax - m.options.xmin, m.options.ymax - m.options.ymin)
    local step = (o.spacing == nil or o.spacing == "auto") and span / 18
        or number(o.spacing, "spacing")
    local length = (o.length == nil or o.length == "auto") and 0.65 * step
        or number(o.length, "length")
    check(step > 0 and length > 0, "spacing and length must be positive")
    local mode = o.mode or "normalized"
    check(mode == "normalized" or mode == "scaled", "unknown vector mode")
    local scale = number(o.scale or 1, "vector scale")
    check(scale > 0, "vector scale must be positive")
    local nx, ny =
        math.floor((m.options.xmax - m.options.xmin) / step),
        math.floor((m.options.ymax - m.options.ymin) / step)
    check((nx + 1) * (ny + 1) <= 100000, "vector sample budget exceeded")
    local out = {}
    for j = 0, ny do
        for i = 0, nx do
            local x, y = m.options.xmin + (i + 0.5) * step, m.options.ymin + (j + 0.5) * step
            local id = V.locate(v, x * m.options.unit, y * m.options.unit)
            if id then
                local f = V.field(v, id, x * m.options.unit, y * m.options.unit)
                local fx, fy = f[v.component .. "x"], f[v.component .. "y"]
                local dx, dy =
                    H.value(fx[1], fx[2], v.representation, v.phase),
                    H.value(fy[1], fy[2], v.representation, v.phase)
                local mag = math.sqrt(dx * dx + dy * dy)
                if mag > 1e-30 then
                    local l = mode == "normalized" and length or scale * mag
                    out[#out + 1] = { x, y, x + dx * l / mag, y + dy * l / mag }
                end
            end
        end
    end
    return out
end
return V
