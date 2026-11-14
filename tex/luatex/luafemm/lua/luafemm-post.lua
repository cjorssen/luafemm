-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Named static integrals and mechanical forces on the unsmoothed field.
-- Selections never alter physical geometry. Results are document-wide snapshots.
-- All integrated outputs use SI; geometry and torque origins use model units.
-- @module luafemm-post
local P = {}
local f = require("luafemm")
local geometry = require("luafemm-geometry")
local selections, results = {}, {}
local function check(v, message)
    assert(v, "luafemm post: " .. message)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
local function name(s)
    check(type(s) == "string" and s:match("^[%w_-]+$"), "invalid name")
end
local function options_copy(options)
    local out = {}
    for key, value in pairs(options or {}) do
        out[key] = value
    end
    return out
end
local function center(m, e)
    local x, y = 0, 0
    for _, id in ipairs(e.ids) do
        x, y = x + m.nodes[id][1] / 3, y + m.nodes[id][2] / 3
    end
    return x, y
end
local function air(e)
    return not e.p.bh and e.p.mur_x == 1 and e.p.mur_y == 1 and e.p.hc == 0 and e.current == 0
end

--- Register a material filter, optionally intersected with closed contours.
-- Options: material (empty = any), contours (model units), fill_rule.
-- Whole elements are selected by centroid; no partial-element clipping occurs.
-- @tparam table m Source model, not necessarily solved.
-- @tparam string id Document-wide selection identifier.
-- @tparam[opt] table o Material and geometric filters; contours are copied.
function P.select(m, id, o)
    name(id)
    check(m, "declare a problem before a selection")
    o = o or {}
    check(
        o.fill_rule == nil or o.fill_rule == "nonzero" or o.fill_rule == "even odd",
        "invalid fill rule"
    )
    local contours
    if o.contours then
        contours = {}
        check(#o.contours > 0, "empty selection path")
        for i, c in ipairs(o.contours) do
            require("luafemm-mesh").validate(c)
            contours[i] = {}
            for j, p in ipairs(c) do
                check(finite(p[1]) and finite(p[2]), "invalid selection coordinate")
                contours[i][j] = { p[1], p[2] }
            end
        end
    end
    selections[id] = {
        model = m,
        material = o.material or "",
        contours = contours,
        fill_rule = o.fill_rule or "nonzero",
    }
end

--- Capture a closed TikZ selection; no mesh or solve is triggered.
function P.capture(m, id, softpath, o)
    check(m and m.frame, "selection requires a femm/problem picture")
    o = options_copy(o)
    o.contours = require("luafemm-path").read_contours(softpath, m.frame, o.tolerance or 0.03)
    P.select(m, id, o)
end

local function selected(m, id)
    local s = id and id ~= "" and selections[id] or nil
    check(not id or id == "" or s, "unknown selection " .. tostring(id))
    check(not s or s.model == m, "selection belongs to another problem")
    check(not s or s.material == "" or m.materials[s.material], "unknown selection material")
    local chosen, count = {}, 0
    for i, e in ipairs(m.triangles) do
        local x, y = center(m, e)
        if
            not s
            or (
                (
                    s.material == ""
                    or e.material == s.material
                    or (m.sliding and m.sliding.aliases[s.material] == e.material)
                )
                and (not s.contours or geometry.contains(s, x / m.options.unit, y / m.options.unit))
            )
        then
            chosen[i], count = true, count + 1
        end
    end
    check(count > 0, "selection contains no elements")
    return chosen, count
end

--- Primitive of H_table(B), exactly integrated per linear/cubic interval.
-- Two-point Gauss integration is exact for the cubic interpolation in B.
-- Coercivity is applied separately in the element energy below.
-- @tparam table p Prepared isotropic material.
-- @tparam number b Nonnegative B magnitude in tesla.
-- @treturn number Integral in J/m^3, before the coercive contribution.
function P.primitive(p, b)
    check(finite(b) and b >= 0, "invalid flux density")
    if not p.bh then
        local nu = f.constitutive(p, b)
        return b * b * nu / (2 * f.mu0)
    end
    local knots = p.curve and p.curve.points or p.bh
    local cuts = { 0 }
    for _, v in ipairs(knots) do
        if v[1] > 0 and v[1] < b then
            cuts[#cuts + 1] = v[1]
        end
    end
    cuts[#cuts + 1] = b
    local value = 0
    for i = 2, #cuts do
        local half = (cuts[i] - cuts[i - 1]) / 2
        local mid = (cuts[i] + cuts[i - 1]) / 2
        for _, sign in ipairs({ -1, 1 }) do
            local x = mid + sign * half / math.sqrt(3)
            local nu = f.constitutive(p, x)
            value = value + half * nu * x / f.mu0
        end
    end
    return value
end

local function tensor(e)
    return (e.bx * e.bx - e.by * e.by) / (2 * f.mu0), e.bx * e.by / f.mu0
end
local function addforce(r, fx, fy, x, y)
    r.Fx, r.Fy = r.Fx + fx, r.Fy + fy
    r.torque = r.torque + (x - r.origin_x) * fy - (y - r.origin_y) * fx
end

-- Integrate -T grad(w); P1 weights give a constant gradient per element.
-- The moving band prescribes these weights; body selection solves for them.
local function stress_weight(m, chosen, r, values)
    local axi = m.options.problem_type == "axisymmetric"
    for i, e in ipairs(m.triangles) do
        if not chosen[i] and air(e) then
            local wx, wy = 0, 0
            for j, id in ipairs(e.ids) do
                wx, wy = wx + values[id] * e.gx[j], wy + values[id] * e.gy[j]
            end
            if axi then
                for _, q in ipairs(e.quadrature) do
                    local br, bz = require("luafemm-axisymmetric").field(m, e, q.r, q.z)
                    r.Fy = r.Fy
                        - 2
                            * math.pi
                            * q.weight
                            * (br * bz * wx + (bz * bz - br * br) * wy / 2)
                            / f.mu0
                end
            else
                local xx, xy = tensor(e)
                local v = e.area * r.depth
                local x, y = center(m, e)
                addforce(r, -v * (xx * wx + xy * wy), -v * (xy * wx - xx * wy), x, y)
            end
        end
    end
end

-- Dirichlet Laplace problem in source-free air, w=1 on the selected body,
-- w=0 on other bodies and exposed boundaries. No weight crosses a body.
local function weighted(m, chosen, r, o)
    local axi = m.options.problem_type == "axisymmetric"
    local one, zero, active = {}, {}, {}
    for id in pairs(m.exterior_nodes) do
        if not axi or m.nodes[id][1] > 0 then
            zero[id] = true
        end
    end
    for i, e in ipairs(m.triangles) do
        for _, id in ipairs(e.ids) do
            if chosen[i] then
                one[id] = true
            elseif not air(e) then
                zero[id] = true
            else
                active[id] = true
            end
        end
    end
    for id in pairs(one) do
        check(
            not zero[id],
            "selected body touches another source/material or domain boundary; surround it with air"
        )
    end
    local dof, values, n = {}, {}, 0
    for id in ipairs(m.nodes) do
        if one[id] then
            values[id] = 1
        elseif zero[id] or not active[id] then
            values[id] = 0
        else
            n = n + 1
            dof[id] = n
        end
    end
    local K, rhs = {}, {}
    for i = 1, n do
        K[i], rhs[i] = {}, 0
    end
    for i, e in ipairs(m.triangles) do
        if not chosen[i] and air(e) then
            -- Area-weighted diffusion, optionally rounded at half-height, is
            -- useful for avoiding tiny corner elements; harmonic is unweighted.
            local rho = o.weighting == "harmonic" and 1 or math.sqrt(e.area)
            for a, ia in ipairs(e.ids) do
                local row = dof[ia]
                if row then
                    for b, ib in ipairs(e.ids) do
                        local k = rho * e.area * (e.gx[a] * e.gx[b] + e.gy[a] * e.gy[b])
                        local col = dof[ib]
                        if col then
                            K[row][col] = (K[row][col] or 0) + k
                        else
                            rhs[row] = rhs[row] - k * values[ib]
                        end
                    end
                end
            end
        end
    end
    local sol, iterations = f.solve_spd(K, rhs, o.tolerance)
    r.weight_iterations = iterations
    for id in ipairs(m.nodes) do
        if dof[id] then
            values[id] = sol[dof[id]]
        end
        if o.weighting == "midgap" then
            values[id] = values[id] >= 0.5 and 1 or 0
        end
    end
    stress_weight(m, chosen, r, values)
end

-- Split a stress-contour segment at every element crossing. Integration is
-- exact for the P1 field and linear torque lever on this computational polyline.
local function stress_segment(m, a, b, orientation, r)
    local u = m.options.unit
    local x, y, dx, dy = a[1] * u, a[2] * u, (b[1] - a[1]) * u, (b[2] - a[2]) * u
    if m.options.problem_type == "axisymmetric" and x == 0 and dx == 0 then
        return
    end
    if dx == 0 and dy == 0 then
        return
    end
    local cuts, intervals = { 0, 1 }, {}
    for _, e in ipairs(m.triangles) do
        local p = m.nodes[e.ids[1]]
        local v = {
            0,
            e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2]),
            e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2]),
        }
        v[1] = 1 - v[2] - v[3]
        local lo, hi = 0, 1
        for j = 1, 3 do
            local slope = e.gx[j] * dx + e.gy[j] * dy
            if math.abs(slope) < 1e-12 then
                if v[j] < -1e-10 then
                    hi = -1
                end
            elseif slope > 0 then
                lo = math.max(lo, -v[j] / slope)
            else
                hi = math.min(hi, -v[j] / slope)
            end
        end
        if hi - lo > 1e-12 then
            intervals[#intervals + 1] = { lo, hi, e }
            cuts[#cuts + 1], cuts[#cuts + 2] = lo, hi
        end
    end
    table.sort(cuts)
    for i = 2, #cuts do
        local lo, hi = cuts[i - 1], cuts[i]
        if hi - lo > 1e-12 then
            local mid, element = (lo + hi) / 2
            for _, v in ipairs(intervals) do
                if mid >= v[1] - 1e-12 and mid <= v[2] + 1e-12 then
                    check(
                        air(v[3]),
                        "Maxwell contour must stay in source-free air, away from material interfaces"
                    )
                    element = element or v[3]
                end
            end
            check(element, "Maxwell contour leaves the mesh")
            if m.options.problem_type == "axisymmetric" then
                local nx, ny = orientation * dy, -orientation * dx
                for _, q in ipairs({
                    { -math.sqrt(3 / 5), 5 / 9 },
                    { 0, 8 / 9 },
                    {
                        math.sqrt(3 / 5),
                        5 / 9,
                    },
                }) do
                    local t = mid + (hi - lo) * q[1] / 2
                    local radius = x + t * dx
                    local br, bz =
                        require("luafemm-axisymmetric").field(m, element, radius, y + t * dy)
                    r.Fy = r.Fy
                        + 2
                            * math.pi
                            * radius
                            * (hi - lo)
                            * q[2]
                            / 2
                            * (br * bz * nx + (bz * bz - br * br) * ny / 2)
                            / f.mu0
                end
            else
                local xx, xy = tensor(element)
                local nx, ny = orientation * dy, -orientation * dx
                local scale = (hi - lo) * r.depth
                addforce(
                    r,
                    scale * (xx * nx + xy * ny),
                    scale * (xy * nx - xx * ny),
                    x + mid * dx,
                    y + mid * dy
                )
            end
        end
    end
end

local function maxwell(m, r, profile)
    local p = require("luafemm-profiles").get(profile)
    check(
        p.model == m and p.offset == 0 and not p.circle,
        "use a zero-offset closed TikZ profile on this model"
    )
    local points = p.points
    local a, b = points[1], points[#points]
    check(a[1] == b[1] and a[2] == b[2], "Maxwell contour must be closed")
    local polygon, area = {}, 0
    for i = 1, #points - 1 do
        polygon[i] = points[i]
        area = area + points[i][1] * points[i + 1][2] - points[i + 1][1] * points[i][2]
    end
    require("luafemm-mesh").validate(polygon)
    check(area ~= 0, "zero-area Maxwell contour")
    for i = 2, #points do
        stress_segment(m, points[i - 1], points[i], area > 0 and 1 or -1, r)
    end
end

--- Compute regional integrals and, optionally, force/torque; solve on demand.
-- Options: name, selection (empty=all), quantity=integrals|force,
-- method=weighted|maxwell|lorentz|band, profile, weighting=midgap|area|harmonic,
-- origin_x/y in model units, tolerance=1e-10, parameter=0 for result plots.
-- Positive depth is required. Permanent-magnet energy uses W(B=0)=0.
-- @tparam table m Source model, solved on demand.
-- @tparam table o Measurement options; name is required.
-- @treturn table Snapshot of SI integrals, fields and optional forces.
function P.measure(m, o)
    o = options_copy(o)
    name(o.name)
    check(m, "declare a problem before measuring")
    local axi = m.options.problem_type == "axisymmetric"
    check(axi or m.options.depth, "set a positive problem depth before measuring")
    local quantity, method = o.quantity or "integrals", o.method or "weighted"
    if method == "band" then
        check(quantity == "force", "method=band requires quantity=force")
        check(m.sliding, "method=band requires a sliding position")
        check(not o.profile or o.profile == "", "method=band does not accept a profile")
        check(not o.selection or o.selection == "", "method=band measures the whole rotor")
    end
    check(
        quantity == "integrals" or quantity == "force" or quantity == "circuit",
        "unknown quantity"
    )
    check(
        method == "weighted" or method == "maxwell" or method == "lorentz" or method == "band",
        "unknown force method"
    )
    o.weighting, o.tolerance = o.weighting or "midgap", o.tolerance or 1e-10
    check(
        o.weighting == "midgap" or o.weighting == "area" or o.weighting == "harmonic",
        "unknown weighting"
    )
    check(finite(o.tolerance) and o.tolerance > 0 and o.tolerance < 1, "invalid weight tolerance")
    for _, k in ipairs({ "origin_x", "origin_y", "parameter" }) do
        check(finite(o[k] or 0), "invalid " .. k)
    end
    if quantity == "force" then
        check(not m.operating, "operating-point force derivatives are not implemented")
        check(m.options.field_model == "open", "mechanical forces require field model=open")
        check(
            method == "maxwell" or method == "band" or (o.selection and o.selection ~= ""),
            "force requires an explicit body selection"
        )
    end
    f.solve(m)
    check(m.solution_signature == require("luafemm-cache").signature(m), "stale solution")
    if quantity == "circuit" then
        local r = require("luafemm-static").result(m, o.circuit, o.winding)
        r.parameter, r.model, r.quantity = o.parameter or 0, m, "circuit"
        results[o.name] = r
        return r
    end
    local chosen, count = selected(m, o.selection)
    if m.options.frequency > 0 then
        check(quantity == "integrals", "harmonic force/torque observables are not implemented")
        local r = require("luafemm-harmonic").integrals(m, chosen)
        r.model, r.parameter, r.quantity, r.elements = m, o.parameter or 0, quantity, count
        results[o.name] = r
        return r
    end
    local r = {
        model = m,
        parameter = o.parameter or 0,
        method = method,
        quantity = quantity,
        elements = count,
        depth = (m.options.depth or 1) * m.options.unit,
        volume = 0,
        origin_x = (o.origin_x or 0) * m.options.unit,
        origin_y = (o.origin_y or 0) * m.options.unit,
        area = 0,
        current = 0,
        energy = 0,
        coenergy = 0,
        Bx = 0,
        By = 0,
        B = 0,
        Hx = 0,
        Hy = 0,
        H = 0,
        x = 0,
        y = 0,
    }
    if quantity == "force" then
        r.Fx, r.Fy, r.torque = 0, 0, 0
    end
    for i, e in ipairs(m.triangles) do
        if chosen[i] then
            if axi then
                require("luafemm-axisymmetric").regional(m, e, r, quantity, method)
            else
                local x, y = center(m, e)
                r.area, r.current = r.area + e.area, r.current + e.area * e.current
                r.x, r.y = r.x + e.area * x, r.y + e.area * y
                local fields = {
                    Bx = e.bx,
                    By = e.by,
                    B = e.b,
                    Hx = e.hx,
                    Hy = e.hy,
                    H = math.sqrt(e.hx ^ 2 + e.hy ^ 2),
                }
                for key, value in pairs(fields) do
                    r[key] = r[key] + e.area * value
                end
                local energy
                if e.p.bh and not m.operating then
                    energy = P.primitive(e.p, e.b)
                else
                    local hx, hy = require("luafemm-static").element_law(m, e, e.bx, e.by)
                    energy = (e.bx * hx + e.by * hy) / (2 * require("luafemm").mu0)
                end
                energy = energy
                    - require("luafemm-static").coercivity(m, e) * (e.mx * e.bx + e.my * e.by)
                r.energy = r.energy + energy * e.area * r.depth
                r.coenergy = r.coenergy + (e.bx * e.hx + e.by * e.hy - energy) * e.area * r.depth
                if quantity == "force" and method == "lorentz" then
                    check(
                        not e.p.bh and e.p.mur_x == 1 and e.p.mur_y == 1 and e.p.hc == 0,
                        "Lorentz body must be nonmagnetic; use stress for iron/magnets"
                    )
                    addforce(
                        r,
                        -e.current * e.by * e.area * r.depth,
                        e.current * e.bx * e.area * r.depth,
                        x,
                        y
                    )
                end
            end
        end
    end
    for _, key in ipairs({ "Bx", "By", "B", "Hx", "Hy", "H", "x", "y" }) do
        r[key] = r[key] / (axi and r.volume or r.area)
    end
    r.x, r.y = r.x / m.options.unit, r.y / m.options.unit
    if not axi then
        r.volume = r.area * r.depth
    end
    if quantity == "force" then
        if method == "weighted" then
            weighted(m, chosen, r, o)
        elseif method == "band" then
            stress_weight(m, {}, r, m.sliding.weights)
        elseif method == "maxwell" then
            maxwell(m, r, o.profile or "")
        end
        r.force = math.sqrt(r.Fx ^ 2 + r.Fy ^ 2)
    end
    if axi then
        r.Br, r.Bz, r.Hr, r.Hz = r.Bx, r.By, r.Hx, r.Hy
        r.Fz = r.Fy
    end
    results[o.name] = r
    return r
end

--- Read a snapshot component; force components are absent for integral-only data.
function P.value(id, component)
    local r = P.get(id)
    check(finite(r[component]), "unavailable result component " .. tostring(component))
    return r[component]
end
--- Format a validated result as an expandable TeX-compatible numeric literal.
function P.formatted(id, component)
    return string.format("%.17g", P.value(id, component))
end
--- Return a stored snapshot (treat as read-only).
function P.get(id)
    return assert(results[id], "luafemm post: unknown result " .. tostring(id))
end
--- TeX-safe SI unit of a public scalar; empty for counts/model parameters.
function P.unit(component)
    local units = {
        area = "m^2",
        volume = "m^3",
        Br = "T",
        Bz = "T",
        Hr = "A/m",
        Hz = "A/m",
        Fz = "N",
        current = "A",
        energy = "J",
        linkage_real = "Wb",
        linkage_imaginary = "Wb",
        voltage_real = "V",
        voltage_imaginary = "V",
        impedance_real = "\\Omega",
        impedance_imaginary = "\\Omega",
        reactive_power = "var",
        joule_loss = "W",
        magnetic_loss = "W",
        hysteresis_loss = "W",
        lamination_loss = "W",
        total_loss = "W",
        eddy_loss = "W",
        current_real = "A",
        current_imaginary = "A",
        coenergy = "J",
        B = "T",
        Bx = "T",
        By = "T",
        H = "A/m",
        Hx = "A/m",
        Hy = "A/m",
        linkage = "Wb",
        inductance = "H",
        incremental_inductance = "H",
        incremental_inductance_real = "H",
        incremental_inductance_imaginary = "H",
        resistance = "\\Omega",
        voltage = "V",
        power = "W",
        Fx = "N",
        Fy = "N",
        force = "N",
        torque = "N\\,m",
    }
    local base, view = component:match("^(%a+)_(%a+)$")
    if
        (view == "amplitude" or view == "rms" or view == "phase")
        and (base == "current" or base == "voltage" or base == "linkage" or base == "impedance")
    then
        return view == "phase" and "{}^\\circ" or units[base .. "_real"]
    end
    return units[component] or ""
end
--- PGFPlots coordinates in the explicit order of the comma-separated names.
function P.coordinates(ids, component, abscissa)
    local out = {}
    for id in ids:gmatch("[^,%s]+") do
        out[#out + 1] =
            string.format("(%.17g,%.17g)", P.value(id, abscissa), P.value(id, component))
    end
    check(#out > 0, "result plot needs at least one result")
    return table.concat(out, " ")
end

--- Emit a force glyph in the saved model frame; length is illustrative unless scale>0.
function P.arrow(m, id, length, scale, x, y, label)
    local r = P.get(id)
    check(m == r.model, "draw force in its source picture")
    check(r.force, "force glyph needs a force result")
    check(
        finite(length) and length > 0 and finite(scale) and scale >= 0,
        "invalid force glyph size"
    )
    if r.force < 1e-30 then
        return
    end
    local function canvas(a, b)
        local z = m.frame
        return string.format(
            "(%.9fpt,%.9fpt)",
            z[1] * a + z[3] * b + z[5],
            z[2] * a + z[4] * b + z[6]
        )
    end
    x, y = x or r.x, y or r.y
    check(finite(x) and finite(y), "invalid force glyph origin")
    local size = scale > 0 and scale or length / r.force
    local text = label == "none" and ""
        or " node[femm force label]{\\luafemmForceLabel{"
            .. string.format("%.12g", r.force)
            .. "}}"
    tex.sprint(
        "\\draw[femm force arrow] "
            .. canvas(x, y)
            .. "--"
            .. canvas(x + size * r.Fx, y + size * r.Fy)
            .. text
            .. ";"
    )
end
return P
