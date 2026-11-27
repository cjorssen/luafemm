-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Coupled harmonic winding networks and integral solid-conductor currents.
-- Field response columns reduce the sparse PDE to a small complex Schur system.
-- Currents/voltages are reconstructed from A on cache restore; no hidden state
-- or native library is required. All internal lengths and terminal values use SI.
-- @module luafemm-ac-circuits
local N = {}
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm AC circuits: " .. s)
end
local function keys(t)
    local a = {}
    for k in pairs(t or {}) do
        a[#a + 1] = k
    end
    table.sort(a)
    return a
end
local function add(t, k, v)
    t[k] = (t[k] or 0) + v
end
local function dot(v, a)
    local s = 0
    for i, x in pairs(v) do
        s = s + x * a[i]
    end
    return s
end
local function div(ar, ai, br, bi)
    local scale = math.max(math.abs(br), math.abs(bi))
    check(scale > 0 and scale < math.huge, "singular circuit system")
    br, bi = br / scale, bi / scale
    local d = br * br + bi * bi
    return (ar / scale * br + ai / scale * bi) / d, (ai / scale * br - ar / scale * bi) / d
end
-- Row/column equilibration precedes complex partial pivoting. Circuit rows mix
-- amperes and volts; an unscaled pivot test would depend on arbitrary SI ratios.
local function dense(ar, ai, br, bi)
    local n = #br
    local r, s, x, y, scale = {}, {}, {}, {}, {}
    for j = 1, n do
        scale[j] = 0
        for i = 1, n do
            scale[j] = math.max(scale[j], math.abs(ar[i][j]), math.abs(ai[i][j]))
        end
        check(scale[j] > 0, "unconnected circuit unknown")
    end
    for i = 1, n do
        r[i], s[i] = {}, {}
        local row = 0
        for j = 1, n do
            r[i][j], s[i][j] = ar[i][j] / scale[j], ai[i][j] / scale[j]
            row = math.max(row, math.abs(r[i][j]), math.abs(s[i][j]))
        end
        check(row > 0, "empty circuit equation")
        for j = 1, n do
            r[i][j], s[i][j] = r[i][j] / row, s[i][j] / row
        end
        x[i], y[i] = br[i] / row, bi[i] / row
    end
    for k = 1, n do
        local pivot, best = k, 0
        for i = k, n do
            local v = r[i][k] ^ 2 + s[i][k] ^ 2
            if v > best then
                pivot, best = i, v
            end
        end
        check(best > 1e-28, "singular or ill-conditioned circuit system")
        r[k], r[pivot], s[k], s[pivot] = r[pivot], r[k], s[pivot], s[k]
        x[k], x[pivot], y[k], y[pivot] = x[pivot], x[k], y[pivot], y[k]
        for i = k + 1, n do
            local a, b = div(r[i][k], s[i][k], r[k][k], s[k][k])
            for j = k + 1, n do
                r[i][j], s[i][j] =
                    r[i][j] - a * r[k][j] + b * s[k][j], s[i][j] - a * s[k][j] - b * r[k][j]
            end
            x[i], y[i] = x[i] - a * x[k] + b * y[k], y[i] - a * y[k] - b * x[k]
            r[i][k], s[i][k] = 0, 0
        end
    end
    for i = n, 1, -1 do
        for j = i + 1, n do
            x[i], y[i] =
                x[i] - r[i][j] * x[j] + s[i][j] * y[j], y[i] - r[i][j] * y[j] - s[i][j] * x[j]
        end
        x[i], y[i] = div(x[i], y[i], r[i][i], s[i][i])
    end
    for i = 1, n do
        x[i], y[i] = x[i] / scale[i], y[i] / scale[i]
    end
    -- A homogeneous current row may have an exactly zero solution component.
    -- Use the normwise backward error of the equilibrated system: a componentwise
    -- denominator made only of that rounded component would demand exact zero.
    local solution_norm = 0
    for j = 1, n do
        solution_norm = math.max(solution_norm, (math.abs(x[j]) + math.abs(y[j])) * scale[j])
    end
    for i = 1, n do
        local er, ei, bound = -br[i], -bi[i], math.abs(br[i]) + math.abs(bi[i])
        for j = 1, n do
            er = er + ar[i][j] * x[j] - ai[i][j] * y[j]
            ei = ei + ar[i][j] * y[j] + ai[i][j] * x[j]
            bound = bound + (math.abs(ar[i][j]) + math.abs(ai[i][j])) / scale[j] * solution_norm
        end
        check(
            math.abs(er) + math.abs(ei) <= 1e-9 * math.max(bound, 1e-30),
            "circuit backward residual exceeds tolerance"
        )
    end
    return x, y
end
-- Shared private equilibrated solve for small electrical constraint systems.
N.dense = dense
--- Whether this model needs coupled terminal/conductor equations.
-- @tparam table m Physical model, before or after meshing.
function N.active(m)
    if m.network then
        return true
    end
    for _, c in pairs(m.circuits or {}) do
        if c.type ~= "series" or c.drive ~= "current" then
            return true
        end
    end
    for _, r in ipairs(m.regions) do
        if
            r.current_constraint == "zero"
            or r.circuit and r.conductor == "solid"
            or m.materials[r.material].wire_model ~= "none"
        then
            return true
        end
    end
    return false
end
--- Check unambiguous source, solid-leg and dimensional contracts.
-- @tparam table m Physical harmonic model.
function N.validate(m)
    for _, r in ipairs(m.regions) do
        if r.current_constraint == "zero" then
            check(
                not r.circuit and r.current == 0 and (r.current_phase or 0) == 0,
                "zero current constraint cannot also prescribe a circuit or current density"
            )
            check(r.conductor ~= "stranded", "zero current constraint requires a solid conductor")
        end
        if r.current_constraint == "zero" or r.circuit and r.conductor == "solid" then
            if m.options.problem_type == "axisymmetric" then
                for _, contour in ipairs(r.contours) do
                    for _, point in ipairs(contour) do
                        check(point[1] > 0, "constrained solid rings must stay away from the axis")
                    end
                end
            end
            check(
                m.materials[r.material].conductivity > 0,
                "constrained solids require positive conductivity"
            )
            check(
                (r.fill or 1) == 1 and (r.length or 0) == 0,
                "solid legs require winding fill=1 and no conductor length override"
            )
            check(not r.circuit or math.abs(r.turns) == 1, "solid legs require turns=1 or turns=-1")
        end
    end
    if N.active(m) then
        check(
            m.options.problem_type ~= "axial",
            "axial diffusion has no longitudinal circuit ports"
        )
        check(
            m.options.problem_type == "axisymmetric" or m.options.depth,
            "set positive depth for coupled planar circuits"
        )
    end
end
--- Build quadrature-consistent field loads, linkage and current functionals.
-- @tparam table m Meshed model with harmonic element regions already identified.
-- @treturn table Derived network, also stored as m.ac_network.
function N.prepare(m)
    local p = { branches = {}, solids = {}, by_region = {}, circuits = {}, loads = {} }
    local axi = m.options.problem_type == "axisymmetric"
    local unit = m.options.unit
    for _, name in ipairs(keys(m.circuits)) do
        local c = m.circuits[name]
        local resolved = m.static_circuits[name]
        local list = {}
        p.circuits[name] = list
        for _, winding in ipairs(keys(resolved.windings)) do
            local old = resolved.windings[winding]
            local b = {
                name = winding,
                circuit = name,
                index = #p.branches + 1,
                link = {},
                resistance = 0,
                stranded_resistance = 0,
                extra_real = 0,
                extra_imaginary = 0,
                solids = {},
                closed = axi or math.abs(old.turns) < 1e-10,
            }
            p.branches[#p.branches + 1] = b
            list[#list + 1] = b
            check(
                b.closed or c.type == "series" and c.drive == "current",
                "voltage/parallel planar windings require closed signed turns"
            )
            for _, id in ipairs(old.legs) do
                p.by_region[id] = { branch = b }
            end
        end
    end
    for id, r in ipairs(m.regions) do
        local b = p.by_region[id] and p.by_region[id].branch
        if r.current_constraint == "zero" or b and r.conductor == "solid" then
            local s = {
                index = #p.branches + #p.solids + 1,
                region = id,
                branch = b,
                sign = r.turns or 1,
                conductance = 0,
                current = {},
                elements = {},
            }
            p.solids[#p.solids + 1] = s
            p.by_region[id] = { solid = s, branch = b }
            if b then
                b.solids[#b.solids + 1] = s
            end
        elseif b then
            local area = m.static_circuits[r.circuit].areas[id]
            local moment = m.static_circuits[r.circuit].moments[id]
            local length = r.length > 0 and r.length * unit
                or axi and 2 * math.pi * moment / area
                or m.options.depth * unit
            local sigma = m.materials[r.material].conductivity
            if sigma == 0 then
                b.missing = true
            else
                local t = m.wire_regions and m.wire_regions[id]
                local factor = r.turns ^ 2 * length / area
                local dc = t and t.dc_resistivity or 1 / (sigma * r.fill)
                b.resistance = b.resistance + factor * dc
                b.stranded_resistance = b.stranded_resistance + factor * dc
                if t then
                    b.extra_real = b.extra_real + factor * (t.resistivity_real - dc)
                    b.extra_imaginary = b.extra_imaginary + factor * t.resistivity_imaginary
                end
            end
        end
    end
    p.size = #p.branches + #p.solids
    p.graph = require("luafemm-network").prepare(m, p.size)
    if p.graph then
        p.size = p.graph.size
    end
    for i = 1, p.size do
        p.loads[i] = {}
    end
    for eid, e in ipairs(m.triangles) do
        local entry = p.by_region[e.region_index]
        if entry then
            local r = m.regions[e.region_index]
            local s, b = entry.solid, entry.branch
            local index = s and s.index or b.index
            if s then
                s.elements[#s.elements + 1] = eid
            end
            for _, q in ipairs(require("luafemm-harmonic").quadrature(m, e)) do
                local length = axi and 2 * math.pi * q.r or m.options.depth * unit
                local density = s and e.sigma / length
                    or r.turns / m.static_circuits[r.circuit].areas[e.region_index]
                if s then
                    s.conductance = s.conductance + q.weight / (axi and q.r or 1) * density
                end
                for j, id in ipairs(e.ids) do
                    local basis = q.n[j]
                    local w = q.weight * basis
                    local dof = m.nodes[id].dof
                    if dof then
                        add(p.loads[index], dof, mu0 * w * density * (m.nodes[id].dof_sign or 1))
                    end
                    if s then
                        add(s.current, id, w / (axi and q.r or 1) * e.sigma)
                    else
                        add(
                            b.link,
                            id,
                            w * (axi and 2 * math.pi or m.options.depth * unit) * density
                        )
                    end
                end
            end
        end
    end
    for _, s in ipairs(p.solids) do
        check(s.conductance > 0, "constrained conductor has no active area")
        if s.branch then
            s.branch.resistance = s.branch.resistance + 1 / s.conductance
        end
        -- Each region is one electrical leg: disconnected painted remnants
        -- would otherwise share an unrequested voltage/zero-current constraint.
        local parents, edges = {}, {}
        local function root(i)
            while parents[i] ~= i do
                parents[i] = parents[parents[i]]
                i = parents[i]
            end
            return i
        end
        for _, eid in ipairs(s.elements) do
            parents[eid] = eid
            local e = m.triangles[eid]
            for j, a in ipairs(e.ids) do
                local b = e.ids[j % 3 + 1]
                local key = math.min(a, b) .. ":" .. math.max(a, b)
                if edges[key] then
                    parents[root(eid)] = root(edges[key])
                else
                    edges[key] = eid
                end
            end
        end
        local first = root(s.elements[1])
        for _, id in ipairs(s.elements) do
            check(
                root(id) == first,
                "a constrained solid region must be connected; declare separate legs"
            )
        end
    end
    m.ac_network = p
    return p
end
-- Evaluate the linear circuit residual E(u,A)-target. Omitting targets builds
-- columns; A includes fixed nodes only for the particular field solution.
local function equations(m, p, ur, ui, ar, ai, targets)
    local omega = 2 * math.pi * m.options.frequency
    local vr, vi = {}, {}
    for _, b in ipairs(p.branches) do
        local r, i =
            (b.stranded_resistance + b.extra_real) * ur[b.index]
                - b.extra_imaginary * ui[b.index]
                - omega * dot(b.link, ai),
            (b.stranded_resistance + b.extra_real) * ui[b.index]
                + b.extra_imaginary * ur[b.index]
                + omega * dot(b.link, ar)
        -- The branch resistance above includes solid DC resistance for reports;
        -- actual solid terminal voltages are independent unknowns here.
        for _, s in ipairs(b.solids) do
            r = r + s.sign * ur[s.index]
            i = i + s.sign * ui[s.index]
        end
        vr[b.index], vi[b.index] = r, i
    end
    local rr, ri = {}, {}
    local function row(r, i)
        rr[#rr + 1] = r
        ri[#ri + 1] = i
    end
    for _, s in ipairs(p.solids) do
        local r = s.conductance * ur[s.index] + omega * dot(s.current, ai)
        local i = s.conductance * ui[s.index] - omega * dot(s.current, ar)
        if s.branch then
            r = r - s.sign * ur[s.branch.index]
            i = i - s.sign * ui[s.branch.index]
        end
        row(r, i)
    end
    for _, name in ipairs(keys(p.circuits)) do
        local c, list = m.circuits[name], p.circuits[name]
        local first = list[1].index
        for j = 2, #list do
            local id = list[j].index
            if c.type == "series" then
                row(ur[id] - ur[first], ui[id] - ui[first])
            else
                row(vr[id] - vr[first], vi[id] - vi[first])
            end
        end
        local r, i = 0, 0
        if c.drive == "network" then
            if c.type == "parallel" then
                r, i = vr[first], vi[first]
            else
                for _, b in ipairs(list) do
                    r, i = r + vr[b.index], i + vi[b.index]
                end
            end
            local graph = require("luafemm-network")
            local e = p.graph.ports[name]
            r = r - graph.voltage(p.graph, e, ur)
            i = i - graph.voltage(p.graph, e, ui)
        elseif c.drive == "linkage" or c.drive == "flux" then
            check(c.type == "series", "linkage/flux drive requires series windings")
            for _, b in ipairs(list) do
                r, i = r + dot(b.link, ar), i + dot(b.link, ai)
                for _, s in ipairs(b.solids) do
                    r = r + s.sign * dot(s.current, ar) / s.conductance
                    i = i + s.sign * dot(s.current, ai) / s.conductance
                end
            end
        elseif c.drive == "current" then
            if c.type == "series" then
                r, i = ur[first], ui[first]
            else
                for _, b in ipairs(list) do
                    r = r + ur[b.index]
                    i = i + ui[b.index]
                end
            end
        elseif c.type == "parallel" then
            r, i = vr[first], vi[first]
        else
            for _, b in ipairs(list) do
                r = r + vr[b.index]
                i = i + vi[b.index]
            end
        end
        if targets and c.drive ~= "network" then
            local amplitude = c.drive == "current" and c.current
                or c.drive == "voltage" and c.voltage
                or c.drive == "flux" and c.flux * c.flux_turns
                or c.linkage
            local phase = math.rad(c.phase)
            r = r - amplitude * math.cos(phase)
            i = i - amplitude * math.sin(phase)
        end
        row(r, i)
    end
    if p.graph then
        require("luafemm-network").equations(m, p.graph, ur, ui, targets, function(name)
            local list = p.circuits[name]
            if m.circuits[name].type == "series" then
                return ur[list[1].index], ui[list[1].index]
            end
            local r, i = 0, 0
            for _, b in ipairs(list) do
                r, i = r + ur[b.index], i + ui[b.index]
            end
            return r, i
        end, row)
    end
    return rr, ri, vr, vi
end
local function zeros(n)
    local t = {}
    for i = 1, n do
        t[i] = 0
    end
    return t
end
local function nodal(m, ar, ai, fixed)
    local r, i = {}, {}
    for id, p in ipairs(m.nodes) do
        r[id] = p.dof and (p.dof_sign or 1) * ar[p.dof] or fixed and m.A[id] or 0
        i[id] = p.dof and (p.dof_sign or 1) * ai[p.dof] or 0
    end
    return r, i
end
--- Solve the coupled field/circuit system from its assembled sparse field block.
-- @tparam table m Model with an ac_network.
-- @tparam table kr Sparse real field matrix.
-- @tparam table ki Sparse imaginary field matrix.
-- @tparam table br Real field load, excluding controlled circuit sources.
-- @tparam table bi Imaginary field load.
-- @return Free-node real/imaginary fields, total GMRES iterations, true residual.
function N.solve(m, kr, ki, br, bi)
    local p = m.ac_network
    local linear = require("luafemm-linear")
    local _, target, options = linear.policy(m.options, "general", m.options.tolerance)
    local tol = math.min(target * 0.01, 1e-12)
    local function solve(a, ai, b, bi, t)
        return require("luafemm-complex").solve(a, ai, b, bi, t, options)
    end
    local ar, ai, it, _, telemetry = solve(kr, ki, br, bi, tol)
    local statistics = linear.accumulate(nil, telemetry)
    local a0, b0 = nodal(m, ar, ai, true)
    local zero = zeros(p.size)
    local rr, ri = equations(m, p, zero, zero, a0, b0, true)
    local mr, mi, xr, xi = {}, {}, {}, {}
    for k = 1, p.size do
        mr[k] = {}
        mi[k] = {}
        rr[k] = -rr[k]
        ri[k] = -ri[k]
    end
    for j = 1, p.size do
        local load = zeros(#br)
        for k, v in pairs(p.loads[j]) do
            load[k] = v
        end
        local steps, residual, telemetry
        xr[j], xi[j], steps, residual, telemetry = solve(kr, ki, load, zeros(#br), tol)
        it = it + steps
        statistics = linear.accumulate(statistics, telemetry)
        local a, b = nodal(m, xr[j], xi[j], false)
        local u = zeros(p.size)
        u[j] = 1
        local r, i = equations(m, p, u, zero, a, b, false)
        for k = 1, p.size do
            mr[k][j], mi[k][j] = r[k], i[k]
        end
    end
    local ur, ui = dense(mr, mi, rr, ri)
    m.circuit_state = { real = ur, imaginary = ui }
    for j = 1, p.size do
        for k = 1, #br do
            ar[k] = ar[k] + xr[j][k] * ur[j] - xi[j][k] * ui[j]
            ai[k] = ai[k] + xr[j][k] * ui[j] + xi[j][k] * ur[j]
        end
    end
    local combined_real, combined_imaginary = {}, {}
    for k = 1, #br do
        local r, i = br[k], bi[k]
        for j = 1, p.size do
            local v = p.loads[j][k] or 0
            r, i = r + v * ur[j], i + v * ui[j]
        end
        combined_real[k], combined_imaginary[k] = r, i
    end
    local operator = linear.prepare(kr, ki, combined_real, combined_imaginary, target)
    local residual, backward = linear.original_residual(operator, ar, ai)
    check(linear.accept(operator, residual, backward), "coupled field residual exceeds tolerance")
    statistics.relative_residual = math.max(statistics.relative_residual, residual)
    statistics.backward_error = math.max(statistics.backward_error, backward)
    return ar, ai, it, residual, statistics
end
--- Reconstruct terminal unknowns from the solved field, including restored caches.
-- @tparam table m Model with A/A_im and a prepared network; updates derived state.
function N.recover(m)
    local p = m.ac_network
    local z = zeros(p.size)
    local rr, ri = equations(m, p, z, z, m.A, m.A_im, true)
    local mr, mi = {}, {}
    for k = 1, p.size do
        mr[k] = {}
        mi[k] = {}
        rr[k] = -rr[k]
        ri[k] = -ri[k]
    end
    local nodes = zeros(#m.nodes)
    for j = 1, p.size do
        local u = zeros(p.size)
        u[j] = 1
        local r, i = equations(m, p, u, z, nodes, nodes, false)
        for k = 1, p.size do
            mr[k][j], mi[k][j] = r[k], i[k]
        end
    end
    local constrained = false
    for _, c in pairs(m.circuits or {}) do
        constrained = constrained or c.drive == "linkage" or c.drive == "flux"
    end
    if constrained then
        check(
            m.circuit_state
                and m.circuit_state.real
                and m.circuit_state.imaginary
                and #m.circuit_state.real == p.size
                and #m.circuit_state.imaginary == p.size,
            "missing flux-driven circuit state"
        )
        p.real, p.imaginary = m.circuit_state.real, m.circuit_state.imaginary
    else
        p.real, p.imaginary = dense(mr, mi, rr, ri)
    end
    local er, ei, vr, vi = equations(m, p, p.real, p.imaginary, m.A, m.A_im, true)
    -- Certify circuit/flux rows separately from the field residual, including
    -- restored unknowns. Scale each row by its actual field, source and terminal
    -- contributions; a small linkage target must not hide behind a volt-scale row.
    local tr, ti = equations(m, p, z, z, nodes, nodes, true)
    local residual = 0
    for k = 1, p.size do
        local scale = math.abs(tr[k])
            + math.abs(ti[k])
            + math.abs(rr[k] + tr[k])
            + math.abs(ri[k] + ti[k])
        for j = 1, p.size do
            scale = scale
                + (math.abs(mr[k][j]) + math.abs(mi[k][j]))
                    * (math.abs(p.real[j]) + math.abs(p.imaginary[j]))
        end
        residual = math.max(residual, math.sqrt(er[k] ^ 2 + ei[k] ^ 2) / math.max(scale, 1e-30))
    end
    check(
        residual <= math.max(1e-9, m.options.tolerance),
        "circuit constraint residual exceeds tolerance"
    )
    m.stats = m.stats or {}
    m.stats.constraint_residual = residual
    for _, b in ipairs(p.branches) do
        b.ir, b.ii = p.real[b.index], p.imaginary[b.index]
        b.vr, b.vi = vr[b.index], vi[b.index]
        b.lr, b.li = dot(b.link, m.A), dot(b.link, m.A_im)
        for _, s in ipairs(b.solids) do
            b.lr = b.lr + s.sign * dot(s.current, m.A) / s.conductance
            b.li = b.li + s.sign * dot(s.current, m.A_im) / s.conductance
        end
    end
    if p.graph then
        require("luafemm-network").results(m, p.graph, p.real, p.imaginary, function(name)
            local list = p.circuits[name]
            if m.circuits[name].type == "series" then
                return list[1].ir, list[1].ii
            end
            local r, i = 0, 0
            for _, b in ipairs(list) do
                r, i = r + b.ir, i + b.ii
            end
            return r, i
        end)
    end
end
--- Evaluate controlled transport and total current at an element point.
-- @tparam table m Solved network model.
-- @tparam table e Element with region_index.
-- @tparam number radius SI radius in axisymmetry (ignored in planar models).
-- @tparam table a Potential phasor {real,imaginary}.
-- @return Total J and zero-net-current redistribution phasors, or nil.
function N.field(m, e, radius, a)
    local entry = m.ac_network.by_region[e.region_index]
    if not entry then
        return
    end
    local s, b = entry.solid, entry.branch
    if not s then
        local r = m.regions[e.region_index]
        local density = r.turns / m.static_circuits[r.circuit].areas[e.region_index]
        return { density * b.ir, density * b.ii }, { 0, 0 }
    end
    local length = m.options.problem_type == "axisymmetric" and 2 * math.pi * radius
        or m.options.depth * m.options.unit
    local density = e.sigma / length
    local w = 2 * math.pi * m.options.frequency * e.sigma
    local jr = density * m.ac_network.real[s.index] + w * a[2]
    local ji = density * m.ac_network.imaginary[s.index] - w * a[1]
    local tr = b and density / s.conductance * s.sign * b.ir or 0
    local ti = b and density / s.conductance * s.sign * b.ii or 0
    return { jr, ji }, { jr - tr, ji - ti }
end
--- Terminal or branch observables for a solved coupled network.
-- Parallel linkage is the DC-conductance-weighted equivalent, not a sum.
-- @tparam table m Solved network model.
-- @tparam string name Circuit name.
-- @tparam[opt] string winding Select one branch; empty selects the whole circuit.
-- @treturn table SI peak phasor scalars and cycle-average powers.
function N.result(m, name, winding)
    local list = m.ac_network.circuits[name]
    check(list, "unknown circuit " .. tostring(name))
    local series = m.circuits[name].type == "series"
    if winding and winding ~= "" then
        local found
        for _, b in ipairs(list) do
            if b.name == winding then
                found = b
            end
        end
        check(found, "unknown winding " .. winding)
        list = { found }
        series = true
    end
    local resistance, missing, closed = 0, false, true
    for _, b in ipairs(list) do
        resistance = resistance + (series and b.resistance or 1 / b.resistance)
        missing = missing or b.missing
        closed = closed and b.closed
    end
    if not series then
        resistance = 1 / resistance
    end
    local ir, ii, vr, vi, lr, li = 0, 0, 0, 0, 0, 0
    for _, b in ipairs(list) do
        if series then
            ir, ii = b.ir, b.ii
            vr, vi = vr + b.vr, vi + b.vi
        else
            ir, ii = ir + b.ir, ii + b.ii
            vr, vi = b.vr, b.vi
        end
        local weight = series and 1 or resistance / b.resistance
        lr, li = lr + weight * b.lr, li + weight * b.li
    end
    local r = {
        current_real = ir,
        current_imaginary = ii,
        resistance = not missing and resistance or nil,
    }
    if closed then
        r.linkage_real, r.linkage_imaginary = lr, li
        if not missing then
            r.voltage_real, r.voltage_imaginary = vr, vi
            r.power, r.reactive_power = (vr * ir + vi * ii) / 2, (vi * ir - vr * ii) / 2
            if ir * ir + ii * ii > 0 then
                r.impedance_real, r.impedance_imaginary = div(vr, vi, ir, ii)
            end
        end
    end
    if m.operating and closed and (not winding or winding == "") and ir * ir + ii * ii > 0 then
        local isolated = true
        for other, c in pairs(m.circuits) do
            isolated = isolated and (other == name or c.current == 0)
        end
        if isolated then
            r.incremental_inductance_real, r.incremental_inductance_imaginary = div(lr, li, ir, ii)
        end
    end
    return N.decorate(r)
end
--- Add peak magnitude, RMS and phase to available terminal phasor pairs.
-- @tparam table r SI terminal result, updated in place; zero phases stay absent.
-- @treturn table The same result with display-friendly scalar components.
function N.decorate(r)
    for _, name in ipairs({ "current", "voltage", "linkage", "impedance" }) do
        local a, b = r[name .. "_real"], r[name .. "_imaginary"]
        if a and b then
            local magnitude = math.sqrt(a * a + b * b)
            r[name .. "_amplitude"] = magnitude
            if name ~= "impedance" then
                r[name .. "_rms"] = magnitude / math.sqrt(2)
            end
            if magnitude > 0 then
                r[name .. "_phase"] = math.deg(math.atan(b, a))
            end
        end
    end
    return r
end
return N
