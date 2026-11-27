-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Named static and harmonic integrals and mechanical forces on the unsmoothed field.
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
    return not e.p.bh
        and e.p.mur_x == 1
        and e.p.mur_y == 1
        and e.p.hc == 0
        and e.current == 0
        and (e.sigma or 0) == 0
        and (e.jre or 0) == 0
        and (e.jim or 0) == 0
        and not e.stranded
        and not e.ac_material
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

-- For peak phasors a,b: <ab>=Re(a conj(b))/2 and (ab)_2=ab/2.
-- A signed imaginary coefficient multiplies -sin(2 omega t), never +sin.
local function product(a, b, part)
    if part == "mean" then
        return (a[1] * b[1] + a[2] * b[2]) / 2
    elseif part == "2f_real" then
        return (a[1] * b[1] - a[2] * b[2]) / 2
    end
    return (a[1] * b[2] + a[2] * b[1]) / 2
end
-- First variation at a captured DC bias: no peak-phasor factor 1/2.
-- The response is at f; the ordinary unbiased quadratic load is at 0 and 2f.
local function response_field(m, e, x, y)
    if m.options.frequency > 0 then
        return require("luafemm-harmonic").field(m, e, x, y)
    end
    local bx, by = e.bx, e.by
    if e.analytical_gap then
        bx, by = require("luafemm-airgap").field(m, m.A, x, y)
    end
    if m.options.problem_type == "axisymmetric" then
        bx, by = require("luafemm-axisymmetric").field(m, e, x, y)
    end
    return { Bx = { bx, 0 }, By = { by, 0 }, J = { e.current, 0 } }
end
local function tensor(m, e, r, x, y)
    if r.response then
        local h = response_field(m, e, x, y)
        local bx, by = require("luafemm-operating").field(m, e, x, y)
        local k = r.part == "1f_imaginary" and 2 or 1
        return (bx * h.Bx[k] - by * h.By[k]) / f.mu0, (bx * h.By[k] + by * h.Bx[k]) / f.mu0
    end
    if r.part then
        local h = require("luafemm-harmonic").field(m, e, x, y)
        return (product(h.Bx, h.Bx, r.part) - product(h.By, h.By, r.part)) / (2 * f.mu0),
            product(h.Bx, h.By, r.part) / f.mu0
    end
    local bx, by = e.bx, e.by
    if e.analytical_gap then
        bx, by = require("luafemm-airgap").field(m, m.A, x, y)
    end
    if m.options.problem_type == "axisymmetric" then
        bx, by = require("luafemm-axisymmetric").field(m, e, x, y)
    end
    return (bx * bx - by * by) / (2 * f.mu0), bx * by / f.mu0
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
                    local xx, xy = tensor(m, e, r, q.r, q.z)
                    r.Fy = r.Fy - 2 * math.pi * q.weight * (xy * wx - xx * wy)
                end
            else
                local x, y = center(m, e)
                local xx, xy = tensor(m, e, r, x, y)
                local v = e.area * r.depth
                addforce(r, -v * (xx * wx + xy * wy), -v * (xy * wx - xx * wy), x, y)
            end
        end
    end
end

-- Analytical bands use a circular stress contour and sufficiently many angular
-- points to integrate every product of retained Fourier modes without aliasing.
local function band_force(m, r)
    local G = require("luafemm-airgap")
    if not G.active(m) then
        return stress_weight(m, {}, r, m.sliding.weights)
    end
    G.contour(m, function(x, y, nx, ny)
        local xx, xy = tensor(m, { analytical_gap = true }, r, x, y)
        addforce(r, r.depth * (xx * nx + xy * ny), r.depth * (xy * nx - xx * ny), x, y)
    end)
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
    local sol, iterations = require("luafemm-linear").real(m, K, rhs, o.tolerance, true)
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
    return values
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
                    local xx, xy = tensor(m, element, r, radius, y + t * dy)
                    r.Fy = r.Fy + 2 * math.pi * radius * (hi - lo) * q[2] / 2 * (xy * nx - xx * ny)
                end
            else
                local xx, xy = tensor(m, element, r, x + mid * dx, y + mid * dy)
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

-- Evaluate the quadratic coefficients or the two first-variation coefficients
-- on the same mesh/selection. Static variations have a zero imaginary part.
-- The real auxiliary weight is solved once, independent of phasor phase.
-- No field, material, potential or cache signature is modified.
local function harmonic_force(m, chosen, r, o, method)
    local H = require("luafemm-harmonic")
    local axi = m.options.problem_type == "axisymmetric"
    check(m.options.problem_type ~= "axial", "axial out-of-plane field forces are not supported")
    r.depth = (m.options.depth or 1) * m.options.unit
    r.origin_x, r.origin_y = (o.origin_x or 0) * m.options.unit, (o.origin_y or 0) * m.options.unit
    r.x, r.y = 0, 0
    for i, e in ipairs(m.triangles) do
        if chosen[i] then
            for _, q in ipairs(H.quadrature(m, e)) do
                local v = q.weight * (axi and 2 * math.pi or r.depth)
                r.x, r.y = r.x + v * q.r, r.y + v * q.z
            end
        end
    end
    r.x, r.y = r.x / r.volume / m.options.unit, r.y / r.volume / m.options.unit
    local weights
    local parts = r.response and { "1f_real", "1f_imaginary" }
        or { "mean", "2f_real", "2f_imaginary" }
    for _, part in ipairs(parts) do
        r.part, r.Fx, r.Fy, r.torque = part, 0, 0, 0
        if method == "weighted" then
            if weights then
                stress_weight(m, chosen, r, weights)
            else
                weights = weighted(m, chosen, r, o)
            end
        elseif method == "band" then
            band_force(m, r)
        elseif method == "maxwell" then
            maxwell(m, r, o.profile or "")
        else
            for i, e in ipairs(m.triangles) do
                if chosen[i] then
                    check(
                        not e.p.bh
                            and e.p.mur_x == 1
                            and e.p.mur_y == 1
                            and e.p.hc == 0
                            and not e.ac_material,
                        "Lorentz body must be nonmagnetic; use stress for iron/magnets"
                    )
                    for _, q in ipairs(H.quadrature(m, e)) do
                        local h = response_field(m, e, q.r, q.z)
                        local function load(b, axis)
                            if not r.response then
                                return product(h.J, b, part)
                            end
                            local bx, by = require("luafemm-operating").field(m, e, q.r, q.z)
                            local k = part == "1f_real" and 1 or 2
                            return m.operating.currents[i] * b[k]
                                + h.J[k] * (axis == "x" and bx or by)
                        end
                        local v = q.weight * (axi and 2 * math.pi or r.depth)
                        if axi then
                            -- e_phi cross e_r = -e_z. Radial resultant cancels
                            -- over a full revolution, as does axial torque.
                            r.Fy = r.Fy - v * load(h.Bx, "x")
                        else
                            addforce(r, -v * load(h.By, "y"), v * load(h.Bx, "x"), q.r, q.z)
                        end
                    end
                end
            end
        end
        for _, k in ipairs({ "Fx", "Fy", "torque" }) do
            r[k .. "_" .. part] = r[k]
        end
        if axi then
            r["Fz_" .. part] = r.Fy
        end
    end
    r.part, r.harmonic_force = nil, m.options.frequency > 0
    r.harmonic_order = r.response and 1 or 2
    r.force_component, r.phase = o.force_component, o.phase
    local prefix = r.response and "_1f_" or "_2f_"
    for _, k in ipairs({ "Fx", "Fy", "torque" }) do
        if r.response then
            r[k .. "_mean"] = 0
        end
        if o.force_component == "instantaneous" then
            r[k] = r[k .. "_mean"]
                + H.value(
                    r[k .. prefix .. "real"],
                    r[k .. prefix .. "imaginary"],
                    "snapshot",
                    r.harmonic_order * o.phase
                )
        elseif o.force_component == "mean" then
            r[k] = r[k .. "_mean"]
        else
            r[k] = r[k .. prefix .. o.force_component]
        end
        r[k .. prefix .. "amplitude"] =
            math.sqrt(r[k .. prefix .. "real"] ^ 2 + r[k .. prefix .. "imaginary"] ^ 2)
    end
    if axi then
        r.Fz, r.Fz_mean = r.Fy, r.Fy_mean
        r["Fz" .. prefix .. "amplitude"] = r["Fy" .. prefix .. "amplitude"]
    end
    r.force = math.sqrt(r.Fx ^ 2 + r.Fy ^ 2)
end

--- Compute regional integrals and, optionally, force/torque; solve on demand.
-- Options: name, selection (empty=all), quantity=integrals|force,
-- method=weighted|maxwell|lorentz|band, profile, weighting=midgap|area|harmonic,
-- origin_x/y in model units, tolerance=1e-10, parameter=0 for result plots.
-- force_component=auto|mean|instantaneous|real|imaginary; phase is in drive degrees.
-- Auto selects the mean for ordinary AC and the real first variation for a
-- captured incremental response. Response harmonics are at f, ordinary loads at 2f.
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
    o.force_component, o.phase = o.force_component or "auto", o.phase or 0
    check(
        o.force_component == "auto"
            or o.force_component == "mean"
            or o.force_component == "instantaneous"
            or o.force_component == "real"
            or o.force_component == "imaginary",
        "unknown force component"
    )
    check(finite(o.phase), "invalid force phase")
    check(
        quantity == "force"
            or ((o.force_component == "mean" or o.force_component == "auto") and o.phase == 0),
        "force component and phase require quantity=force"
    )
    check(
        m.options.frequency > 0
            or ((o.force_component == "mean" or o.force_component == "auto") and o.phase == 0),
        "force component and phase require a harmonic model"
    )
    if method == "band" then
        check(quantity == "force", "method=band requires quantity=force")
        check(
            m.sliding or m.native_age,
            "method=band requires a sliding position or native AGE snapshot"
        )
        check(not o.profile or o.profile == "", "method=band does not accept a profile")
        check(not o.selection or o.selection == "", "method=band measures the whole rotor")
    end
    check(
        quantity == "integrals"
            or quantity == "force"
            or quantity == "circuit"
            or quantity == "network",
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
        check(
            not m.operating or m.operating.mode == "incremental",
            "force derivatives require permeability=incremental"
        )
        check(
            not require("luafemm-fundamental").active(m)
                and not require("luafemm-equivalent").active(m),
            "forces from a fundamental nonlinear approximation are not supported"
        )
        check(m.options.field_model == "open", "mechanical forces require field model=open")
        check(
            method == "maxwell" or method == "band" or (o.selection and o.selection ~= ""),
            "force requires an explicit body selection"
        )
    end
    local response = quantity == "force" and m.operating and m.operating.mode == "incremental"
    if o.force_component == "auto" or (response and m.options.frequency == 0) then
        o.force_component = response and "real" or "mean"
    end
    f.solve(m)
    check(m.solution_signature == require("luafemm-cache").signature(m), "stale solution")
    if quantity == "network" then
        local n = m.network_result
        check(n, "no solved electrical network")
        local element, node = o.element or "", o.node or ""
        check(element == "" or node == "", "select an element or a node, not both")
        local source
        if element ~= "" then
            source = n.elements[element]
        elseif node ~= "" then
            source = n.nodes[node]
        else
            source = n
        end
        check(source, "unknown electrical element/node")
        local r = {
            model = m,
            quantity = "network",
            parameter = o.parameter or 0,
            element = element,
            node = node,
        }
        for k, v in pairs(source) do
            if finite(v) then
                r[k] = v
            end
        end
        if o.register ~= false then
            results[o.name] = r
        end
        return r
    end
    if quantity == "circuit" then
        local r = require("luafemm-static").result(m, o.circuit, o.winding)
        r.parameter, r.model, r.quantity, r.circuit = o.parameter or 0, m, "circuit", o.circuit
        if o.register ~= false then
            results[o.name] = r
        end
        return r
    end
    local chosen, count = selected(m, o.selection)
    if m.options.frequency > 0 then
        local r = require("luafemm-harmonic").integrals(m, chosen)
        r.model, r.parameter, r.quantity, r.elements = m, o.parameter or 0, quantity, count
        if quantity == "force" then
            r.method, r.response = method, response
            harmonic_force(m, chosen, r, o, method)
        end
        if o.register ~= false then
            results[o.name] = r
        end
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
                r.area, r.current =
                    r.area + e.area,
                    r.current + (m.transient and e.time_current or e.area * e.current)
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
                if e.p.nonlinear_sheet and not m.operating then
                    energy = require("luafemm-lamination").energy(e.p, e.bx, e.by)
                elseif e.p.bh and not m.operating then
                    energy = P.primitive(e.p, e.b)
                else
                    local hx, hy = require("luafemm-static").element_law(m, e, e.bx, e.by)
                    energy = (e.bx * hx + e.by * hy) / (2 * require("luafemm").mu0)
                end
                energy = energy
                    - require("luafemm-static").coercivity(m, e) * (e.mx * e.bx + e.my * e.by)
                r.energy = r.energy + (e.analytical_gap and 0 or energy * e.area * r.depth)
                r.coenergy = r.coenergy
                    + (
                        e.analytical_gap and 0
                        or (e.bx * e.hx + e.by * e.hy - energy) * e.area * r.depth
                    )
                if quantity == "force" and method == "lorentz" and not response then
                    check(
                        not e.p.bh and e.p.mur_x == 1 and e.p.mur_y == 1 and e.p.hc == 0,
                        "Lorentz body must be nonmagnetic; use stress for iron/magnets"
                    )
                    if m.transient then
                        for _, q in ipairs(require("luafemm-harmonic").quadrature(m, e)) do
                            local j = require("luafemm-time").field(m, e, q.r, q.z).J
                            addforce(
                                r,
                                -j * e.by * q.weight * r.depth,
                                j * e.bx * q.weight * r.depth,
                                q.r,
                                q.z
                            )
                        end
                    else
                        addforce(
                            r,
                            -(m.transient and e.time_current or e.current * e.area) * e.by * r.depth,
                            (m.transient and e.time_current or e.current * e.area) * e.bx * r.depth,
                            x,
                            y
                        )
                    end
                end
            end
        end
    end
    local gap_energy = require("luafemm-airgap").selected_energy(m, chosen)
    if m.transient then
        r.time, r.time_step, r.joule_loss = m.transient.time, m.transient.dt, 0
        for i, e in ipairs(m.triangles) do
            if chosen[i] then
                r.joule_loss = r.joule_loss + e.time_loss
            end
        end
        if m.time_motion then
            r.angle, r.speed = m.time_motion.angle, m.time_motion.speed
        end
    end
    r.energy, r.coenergy = r.energy + gap_energy, r.coenergy + gap_energy
    for _, key in ipairs({ "Bx", "By", "B", "Hx", "Hy", "H", "x", "y" }) do
        r[key] = r[key] / (axi and r.volume or r.area)
    end
    r.x, r.y = r.x / m.options.unit, r.y / m.options.unit
    if not axi then
        r.volume = r.area * r.depth
    end
    if quantity == "force" then
        if response then
            r.response = true
            harmonic_force(m, chosen, r, o, method)
        elseif method == "weighted" then
            weighted(m, chosen, r, o)
        elseif method == "band" then
            band_force(m, r)
        elseif method == "maxwell" then
            maxwell(m, r, o.profile or "")
        end
        r.force = math.sqrt(r.Fx ^ 2 + r.Fy ^ 2)
    end
    if axi then
        r.Br, r.Bz, r.Hr, r.Hz = r.Bx, r.By, r.Hx, r.Hy
        r.Fz = r.Fy
    end
    if o.register ~= false then
        results[o.name] = r
    end
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
        time = "s",
        time_step = "s",
        speed = "rad/s",
        angle = "deg",
        induced_voltage = "V",
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
        wire_transport_loss = "W",
        proximity_loss = "W",
        wire_energy = "J",
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
    local mechanical = component:match("^(%a+)_mean$") or component:match("^(%a+)_[12]f_%a+$")
    if mechanical == "Fx" or mechanical == "Fy" or mechanical == "Fz" or mechanical == "torque" then
        return units[mechanical]
    end
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

--- Reconstruct a harmonic mechanical scalar at temporal phase omega*t (degrees).
-- Stored snapshots are immutable; this does not solve or change the source model.
function P.instantaneous(id, component, phase)
    local r = P.get(id)
    check(r.harmonic_force, "force waveform requires a harmonic force result")
    check(
        component == "Fx"
            or component == "Fy"
            or component == "Fz"
            or component == "torque"
            or component == "force",
        "unknown force waveform component"
    )
    check(finite(phase), "invalid force phase")
    if component == "force" then
        return math.sqrt(
            P.instantaneous(id, "Fx", phase) ^ 2 + P.instantaneous(id, "Fy", phase) ^ 2
        )
    end
    local prefix = r.response and "_1f_" or "_2f_"
    return P.value(id, component .. "_mean")
        + require("luafemm-harmonic").value(
            P.value(id, component .. prefix .. "real"),
            P.value(id, component .. prefix .. "imaginary"),
            "snapshot",
            (r.harmonic_order or 2) * phase
        )
end

--- PGFPlots samples of one saved force waveform; angles use the drive frequency.
function P.force_coordinates(id, component, first, last, samples)
    check(finite(first) and finite(last) and last > first, "invalid force phase interval")
    check(
        finite(samples) and samples == math.floor(samples) and samples >= 2 and samples <= 10001,
        "force plot samples must be an integer in [2,10001]"
    )
    local out = {}
    for i = 0, samples - 1 do
        local phase = first + (last - first) * i / (samples - 1)
        out[#out + 1] = string.format("(%.17g,%.17g)", phase, P.instantaneous(id, component, phase))
    end
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
