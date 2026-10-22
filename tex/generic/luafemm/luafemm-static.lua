-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Static tensor materials, spatial magnetisation and named DC windings.
-- All orientations are measured in the physical model frame, in degrees.
-- Public geometry uses model units; electrical properties use SI.
-- @module luafemm-static
local S = {}
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm static: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
local function number(v, default, name)
    if v == nil then
        v = default
    end
    check(finite(v), "invalid " .. name)
    return v
end

-- Deterministic accumulation keeps passive cache descriptors reproducible
-- across LuaTeX formats and independently seeded Lua string hash tables.
local function keys(t)
    local list = {}
    for key in pairs(t) do
        list[#list + 1] = key
    end
    table.sort(list)
    return list
end

--- Prepare a symmetric positive tensor from principal permeabilities.
-- Linear lamination: arithmetic mixing along sheets, harmonic across them.
-- Nonlinear anisotropy/lamination is deliberately rejected, not approximated
-- by independent scalar B-H curves. Conductivity is in S/m.
-- @tparam table p Owned material with scalar defaults already validated.
function S.material(p)
    local x = number(p.mur_x, p.mur, "mu x")
    local y = number(p.mur_y, p.mur, "mu y")
    local angle = number(p.material_angle, 0, "material angle")
    local fill = number(p.lamination_fill, 1, "lamination fill")
    local lam = p.lamination or "none"
    check(x > 0 and y > 0, "principal permeabilities must be positive")
    check(fill > 0 and fill <= 1, "lamination fill must be in (0,1]")
    check(lam == "none" or lam == "parallel x" or lam == "parallel y", "unknown lamination")
    check(lam ~= "none" or fill == 1, "lamination fill requires a sheet direction")
    check(
        not p.bh or (x == y and x == p.mur and fill == 1),
        "nonlinear anisotropy/lamination is unsupported"
    )
    check(fill == 1 or p.hc == 0, "laminated permanent magnets are unsupported")
    if lam == "parallel x" then
        x, y = fill * x + 1 - fill, 1 / (fill / y + 1 - fill)
    elseif lam == "parallel y" then
        x, y = 1 / (fill / x + 1 - fill), fill * y + 1 - fill
    end
    check(p.remanence == 0 or x == y, "anisotropic magnets require coercivity, not remanence")
    if p.remanence > 0 then
        p.hc = p.remanence / (mu0 * x)
    end
    check(finite(p.hc), "non-finite coercivity")
    p.mur_x, p.mur_y, p.material_angle = x, y, angle
    p.lamination, p.lamination_fill = lam, fill
    p.conductivity = number(p.conductivity, 0, "conductivity")
    check(p.conductivity >= 0, "conductivity must be nonnegative")
    local c, s = math.cos(math.rad(angle)), math.sin(math.rad(angle))
    p.tensor = { c * c / x + s * s / y, c * s * (1 / x - 1 / y), s * s / x + c * c / y }
end

--- Scaled H and its symmetric tangent, excluding the coercive source.
-- @tparam table p Prepared material.
-- @tparam number bx Model-frame Bx in tesla.
-- @tparam number by Model-frame By in tesla.
-- @return mu0*Hx, mu0*Hy, tangent xx, xy, yy, before coercivity.
function S.law(p, bx, by)
    if not p.bh then
        local t = p.tensor
        return t[1] * bx + t[2] * by, t[2] * bx + t[3] * by, t[1], t[2], t[3]
    end
    local b = math.sqrt(bx * bx + by * by)
    local nu, d = require("luafemm").constitutive(p, b)
    local k = b > 1e-15 and (d - nu) / (b * b) or 0
    return nu * bx, nu * by, nu + k * bx * bx, k * bx * by, nu + k * by * by
end

--- Configure a region's spatial magnetisation without executing source text.
-- angle = offset + order*atan2(y-cy,x-cx); radial/tangential are presets.
-- This bounded family includes cylindrical Halbach angle laws.
-- @tparam table r Owned region returned by luafemm.region, before meshing.
-- @tparam[opt] table o Magnetisation and winding data; lengths in model units.
function S.region(r, o)
    o = o or {}
    local kind = o.magnetization or "uniform"
    check(
        kind == "uniform" or kind == "radial" or kind == "tangential" or kind == "angular",
        "unknown magnetization"
    )
    r.magnetization = {
        type = kind,
        x = number(o.center_x, 0, "magnetization center x"),
        y = number(o.center_y, 0, "magnetization center y"),
        order = number(o.order, 1, "magnetization order"),
    }
    if o.circuit and o.circuit ~= "" then
        check(type(o.circuit) == "string" and o.circuit:match("%S"), "invalid circuit name")
        r.circuit = o.circuit
        r.winding = o.winding and o.winding ~= "" and o.winding or "main"
        check(type(r.winding) == "string" and r.winding:match("%S"), "invalid winding name")
        r.turns = number(o.turns, 1, "turns")
        r.fill = number(o.fill, 1, "winding fill")
        r.length = number(o.length, 0, "conductor length")
        check(
            r.turns ~= 0 and r.fill > 0 and r.fill <= 1 and r.length >= 0,
            "invalid winding geometry"
        )
        check(r.current == 0, "circuit regions cannot also prescribe current density")
    end
end

--- Average the coercive direction with symmetric three-point area quadrature.
-- P1 test curls are constant; the averaged source is also used in energy/H.
-- @tparam table r Region with a prepared angular law.
-- @tparam table vertices Three coordinate pairs in model units.
-- @return Quadrature-averaged x and y components of the unit direction.
function S.direction(r, vertices)
    local law = r.magnetization
    if not law or law.type == "uniform" then
        return r.mx, r.my
    end
    local mx, my = 0, 0
    for k = 1, 3 do
        local x, y = 0, 0
        for j = 1, 3 do
            local w = j == k and 2 / 3 or 1 / 6
            x, y = x + w * vertices[j][1], y + w * vertices[j][2]
        end
        x, y = x - law.x, y - law.y
        check(x * x + y * y > 1e-28, "magnetization is undefined at its centre")
        local order = law.type == "angular" and law.order or 1
        local a = math.rad(r.angle + (law.type == "tangential" and 90 or 0))
            + order * math.atan(y, x)
        mx, my = mx + math.cos(a) / 3, my + math.sin(a) / 3
    end
    return mx, my
end

--- Declare a named DC series or parallel circuit before meshing.
-- Voltage drive requires conductivity and a conductor length for every leg.
-- @tparam table m Unmeshed model.
-- @tparam string name Problem-local unique circuit name.
-- @tparam[opt] table o type, drive, current (A), voltage (V).
function S.circuit(m, name, o)
    check(not m.nodes, "declare circuits before meshing")
    check(type(name) == "string" and name:match("%S"), "circuit needs a name")
    o = o or {}
    local c = {
        type = o.type or "series",
        drive = o.drive or "current",
        current = number(o.current, 0, "circuit current"),
        voltage = number(o.voltage, 0, "voltage"),
    }
    check(c.type == "series" or c.type == "parallel", "unknown circuit type")
    check(c.drive == "current" or c.drive == "voltage", "unknown circuit drive")
    check(
        c.drive == "current" and c.voltage == 0 or c.drive == "voltage" and c.current == 0,
        "specify only the driven quantity"
    )
    m.circuits = m.circuits or {}
    check(not m.circuits[name], "duplicate circuit " .. name)
    m.circuits[name] = c
end

--- Resolve active region areas and DC branch currents on the final arrangement.
-- A winding groups its outgoing and returning legs. Paint order is respected.
-- @tparam table m Model with input regions, materials and circuit declarations.
-- @tparam table t Final arrangement; face currents are updated in place.
function S.excite(m, t)
    local areas, groups = {}, {}
    for _, f in ipairs(t.faces) do
        if not f.hole and f.region then
            areas[f.region] = (areas[f.region] or 0) + f.area * m.options.unit ^ 2
        end
    end
    for i, r in ipairs(m.regions) do
        if r.circuit then
            check(m.circuits and m.circuits[r.circuit], "unknown circuit " .. r.circuit)
            check(areas[i] and areas[i] > 0, "winding region has no active area")
            local cg = groups[r.circuit] or {}
            groups[r.circuit] = cg
            local g = cg[r.winding] or { resistance = 0, legs = {}, turns = 0 }
            cg[r.winding] = g
            g.legs[#g.legs + 1] = i
            g.turns = g.turns + r.turns
            local sigma = m.materials[r.material].conductivity
            local length = r.length > 0 and r.length * m.options.unit
                or (m.options.depth or 0) * m.options.unit
            local resistance = sigma > 0
                    and length > 0
                    and r.turns ^ 2 * length / (sigma * r.fill * areas[i])
                or nil
            if not resistance then
                g.missing = true
            else
                check(
                    finite(resistance) and resistance > 0,
                    "non-finite or underflowed leg resistance"
                )
                g.resistance = g.resistance + resistance
            end
        end
    end
    local currents, resolved = {}, {}
    for _, name in ipairs(keys(m.circuits or {})) do
        local c = m.circuits[name]
        local gs = groups[name]
        check(gs, "circuit has no winding: " .. name)
        local sum, missing = 0, false
        for _, winding in ipairs(keys(gs)) do
            local g = gs[winding]
            missing = missing or g.missing
            if not g.missing then
                sum = sum + (c.type == "series" and g.resistance or 1 / g.resistance)
            end
        end
        check(
            not missing or (c.type == "series" and c.drive == "current"),
            "voltage/parallel circuits require conductivity and conductor length"
        )
        local resistance = not missing and (c.type == "series" and sum or 1 / sum) or nil
        check(
            not resistance or (finite(resistance) and resistance > 0),
            "invalid circuit resistance"
        )
        local current = c.drive == "current" and c.current or c.voltage / resistance
        local voltage = c.drive == "voltage" and c.voltage or resistance and current * resistance
        check(
            finite(current) and (not voltage or finite(voltage)),
            "non-finite circuit operating point"
        )
        for _, winding in ipairs(keys(gs)) do
            local g = gs[winding]
            g.current = c.type == "series" and current or voltage / g.resistance
            for _, i in ipairs(g.legs) do
                currents[i] = m.regions[i].turns * g.current / areas[i]
                check(finite(currents[i]), "non-finite circuit current density")
            end
        end
        resolved[name] = {
            current = current,
            voltage = voltage,
            resistance = resistance,
            windings = gs,
            areas = areas,
        }
    end
    for _, f in ipairs(t.faces) do
        if currents[f.region] then
            f.current = currents[f.region]
        end
    end
    m.static_circuits = resolved
end

--- Return SI terminal quantities and linkage for closed signed-turn windings.
-- Linkage is depth*sum(N/area*integral(A dS)); no gauge-dependent open leg
-- is reported as an inductance. Secant inductance excludes permanent magnets.
-- @tparam table m Converged model with an unchanged physical signature.
-- @tparam string name Circuit identifier.
-- @treturn table SI terminal quantities; unavailable observables are absent.
function S.result(m, name)
    check(m.stats, "solve before measuring a circuit")
    check(m.solution_signature == require("luafemm-cache").signature(m), "stale circuit solution")
    local c = m.static_circuits and m.static_circuits[name]
    check(c, "unknown or unresolved circuit " .. tostring(name))
    check(m.options.depth, "set problem depth for flux linkage")
    local linkage = {}
    for _, e in ipairs(m.triangles) do
        local face = m.topology.faces[e.face]
        local r = face.region and m.regions[face.region]
        if r and r.circuit == name then
            local a = (m.A[e.ids[1]] + m.A[e.ids[2]] + m.A[e.ids[3]]) / 3
            linkage[r.winding] = (linkage[r.winding] or 0)
                + r.turns * a * e.area / c.areas[face.region] * m.options.depth * m.options.unit
        end
    end
    local total, closed = 0, true
    local series = m.circuits[name].type == "series"
    for _, key in ipairs(keys(c.windings)) do
        local g = c.windings[key]
        closed = closed and math.abs(g.turns) < 1e-10
        local weight = series and 1 or c.resistance / g.resistance
        total = total + weight * (linkage[key] or 0)
    end
    local r = {
        current = c.current,
        voltage = c.voltage,
        resistance = c.resistance,
        power = c.voltage and c.voltage * c.current,
        windings = linkage,
    }
    if closed then
        r.linkage = total
        local magnets = false
        for _, e in ipairs(m.triangles) do
            magnets = magnets or e.p.hc > 0
        end
        if c.current ~= 0 and not magnets then
            r.inductance = total / c.current
        end
    end
    return r
end
return S
