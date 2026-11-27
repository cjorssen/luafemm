-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Physical eddy-current views: transverse streamlines, normal glyphs and rings.
-- Returns geometry only. The TeX adapter owns styles and picture transforms.
-- @module luafemm-currents
local C = {}
local function check(v, s)
    assert(v, "luafemm currents: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
local function chosen(e, o)
    return o.material == "" or e.material == o.material
end
local function element(m, x, y)
    x, y = x * m.options.unit, y * m.options.unit
    for _, e in ipairs(m.triangles) do
        local p = m.nodes[e.ids[1]]
        local b = e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2])
        local c = e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
        if b >= -1e-10 and c >= -1e-10 and b + c <= 1 + 1e-10 then
            return e, x, y
        end
    end
end
--- Produce an oriented, phase-specific view of induced current only.
-- @tparam table m Harmonic model, solved on demand.
-- @tparam[opt] table o view, phase, material, lines, spacing, section,
-- reference, minimum_fraction and centre coordinates, in model units.
-- @treturn table View with paths or signed symbols; original solution unchanged.
function C.geometry(m, o)
    o = o or {}
    local p = {
        view = o.view or "auto",
        phase = o.phase or 90,
        material = o.material or "",
        lines = o.lines or 7,
        spacing = o.spacing or 2,
        section = o.section or 0,
        reference = o.reference or 0,
        minimum_fraction = o.minimum_fraction or 0.02,
        center_x = o.center_x or 0,
        center_y = o.center_y or 0,
    }
    check(
        m and (m.options.frequency > 0 or m.transient),
        "eddy currents require frequency>0 or a time step"
    )
    for _, k in ipairs({
        "phase",
        "lines",
        "spacing",
        "section",
        "reference",
        "minimum_fraction",
        "center_x",
        "center_y",
    }) do
        check(finite(p[k]), "invalid " .. k)
    end
    check(
        p.lines >= 1 and p.lines <= 200 and p.lines % 1 == 0,
        "lines must be an integer from 1 to 200"
    )
    check(
        p.spacing > 0 and p.reference >= 0 and p.minimum_fraction >= 0 and p.minimum_fraction < 1,
        "invalid spacing or threshold"
    )
    check(p.material == "" or m.materials[p.material], "unknown material " .. p.material)
    local axial = m.options.problem_type == "axial"
    local axi = m.options.problem_type == "axisymmetric"
    if p.view == "auto" then
        p.view = axial and "lines" or "symbols"
    end
    check(p.view == "lines" or p.view == "symbols" or p.view == "loops", "unknown current view")
    check(
        p.view ~= "lines" or axial,
        "in-plane current lines require type=axial; use symbols or axisymmetric loops"
    )
    check(p.view ~= "symbols" or not axial, "axial currents lie in the plane; use view=lines")
    check(p.view ~= "loops" or axi, "circular current loops require type=axisymmetric")
    require("luafemm").solve(m)
    check(m.solution_signature == require("luafemm-cache").signature(m), "stale current solution")
    local harmonic = require("luafemm-harmonic")
    local out = { view = p.view, paths = {}, symbols = {}, reference = 0 }
    if p.view == "lines" then
        local s = {
            options = { unit = m.options.unit },
            nodes = m.nodes,
            triangles = {},
            A = {},
            stats = {},
        }
        local angle = math.rad(m.options.applied_phase + p.phase)
        local h = m.options.applied_field / (4 * math.pi * 1e-7)
        for i, a in ipairs(m.A) do
            s.A[i] = h * (a * math.cos(angle) - m.A_im[i] * math.sin(angle))
        end
        local low, high = math.huge, -math.huge
        for _, e in ipairs(m.triangles) do
            if e.p.conductivity > 0 and chosen(e, p) then
                local jx, jy = 0, 0
                for i, id in ipairs(e.ids) do
                    jx = jx + e.gy[i] * s.A[id]
                    jy = jy - e.gx[i] * s.A[id]
                    low = math.min(low, s.A[id])
                    high = math.max(high, s.A[id])
                end
                s.triangles[#s.triangles + 1] = { ids = e.ids, bx = jx, by = jy }
            end
        end
        if high - low > math.max(math.abs(h) * 1e-13, 1e-30) then
            local levels = {}
            for i = 1, p.lines do
                levels[i] = low + (high - low) * i / (p.lines + 1)
            end
            out.paths = require("luafemm").contours(s, p.lines, levels)
        end
        return out
    end
    local samples = {}
    local function add(x, y)
        local e, sx, sy = element(m, x, y)
        if e and e.sigma > 0 and chosen(e, p) then
            local j
            if m.transient then
                j = require("luafemm-time").field(m, e, sx, sy).Jeddy
            else
                local f = harmonic.field(m, e, sx, sy)
                j = harmonic.value(f.Jeddy[1], f.Jeddy[2], "snapshot", p.phase)
            end
            samples[#samples + 1] = { x = x, y = y, j = j }
            out.reference = math.max(out.reference, math.abs(j))
        end
    end
    local o = m.options
    if p.view == "symbols" then
        local nx, ny =
            math.ceil((o.xmax - o.xmin) / p.spacing), math.ceil((o.ymax - o.ymin) / p.spacing)
        check(nx * ny <= 10000, "too many current symbols; increase spacing")
        for i = 0, nx - 1 do
            for j = 0, ny - 1 do
                add(
                    o.xmin + (i + 0.5) * (o.xmax - o.xmin) / nx,
                    o.ymin + (j + 0.5) * (o.ymax - o.ymin) / ny
                )
            end
        end
    else
        check(p.section >= o.ymin and p.section <= o.ymax, "section outside axial bounds")
        local n = math.ceil((o.xmax - o.xmin) / p.spacing)
        check(n <= 1000, "too many current loops; increase spacing")
        for i = 0, n - 1 do
            add(o.xmin + (i + 0.5) * (o.xmax - o.xmin) / n, p.section)
        end
    end
    out.reference = p.reference > 0 and p.reference or out.reference
    if out.reference == 0 then
        return out
    end
    for _, sample in ipairs(samples) do
        local fraction = math.min(1, math.abs(sample.j) / out.reference)
        if fraction > p.minimum_fraction then
            if p.view == "symbols" then
                sample.positive = sample.j > 0
                -- In a right-handed r,z picture, +phi points INTO the page.
                sample.dot = axi and sample.j < 0 or not axi and sample.j > 0
                sample.opacity = 0.2 + 0.8 * fraction
                out.symbols[#out.symbols + 1] = sample
            else
                local path = { opacity = 0.2 + 0.8 * fraction }
                for k = 0, 96 do
                    local a = (sample.j > 0 and 1 or -1) * 2 * math.pi * k / 96
                    path[#path + 1] = {
                        x = p.center_x + sample.x * math.cos(a),
                        y = p.center_y + sample.x * math.sin(a),
                    }
                end
                out.paths[#out.paths + 1] = path
            end
        end
    end
    return out
end
return C
