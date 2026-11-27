-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Real backward-Euler magnetic diffusion and material-coordinate motion.
-- Assembly uses SI A (T m), conductivity (S/m), and time (s). The magnetic
-- residual is scaled by mu0, exactly like the static nonlinear solver.
-- No harmonic effective permeability or frequency-domain winding impedance
-- is interpreted as a time-domain constitutive law.
local T = {}
local S = require("luafemm-time-state")
local copy = S.copy
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm transient: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
--- Common physical contract, checked before cache lookup and state capture.
function T.eligible(m)
    check(
        m.options.problem_type == "axisymmetric" or m.options.depth,
        "planar time states require problem depth"
    )
    check(
        m.options.frequency == 0 and m.options.problem_type ~= "axial",
        "requires a real planar or axisymmetric problem"
    )
    check(
        m.options.field_model == "open" and not m.operating,
        "responses and ideal circuits cannot be time states"
    )
    check(
        not m.network and not m.native_age,
        "external networks and native AGE snapshots are unsupported"
    )
    check(
        not m.options.source and not m.options.boundary,
        "callbacks cannot become persistent time states"
    )
    -- The static exterior audit checks impressed current only. Extending that
    -- approximation to time fields requires an audit of induced return currents.
    for _, frontier in pairs(m.frontiers or {}) do
        check(
            not (frontier.condition and frontier.condition.open_radius),
            "open multipole boundaries are unsupported in time; use a finite exterior domain"
        )
    end
    for _, p in pairs(m.materials) do
        check(
            not require("luafemm-losses").active(p)
                and p.wire_model == "none"
                and p.harmonic_law == "standard",
            "harmonic effective losses/wires require a separate time-domain law"
        )
        check(
            p.lamination == "none" or p.conductivity == 0,
            "conducting homogenised sheets need a time-domain conductivity law; resolve the sheets"
        )
    end
    for _, c in pairs(m.circuits or {}) do
        check(
            c.type == "series" and c.drive == "current",
            "time circuits require series current drive"
        )
    end
    for _, r in ipairs(m.regions) do
        check(
            not r.circuit or r.conductor ~= "solid",
            "solid circuit terminals are not supported in time"
        )
        if r.current_constraint == "zero" then
            check(
                m.options.problem_type == "planar" and not r.circuit and r.current == 0,
                "zero-current constraints require source-free planar solid regions"
            )
        end
    end
    if m.sliding then
        check(
            m.sliding.options.method == "p1" and (m.sliding.options.sectors or 1) == 1,
            "time motion requires a full-circle P1 band"
        )
    end
end
local function scale_boundary(c, scale)
    for _, k in ipairs({ "value", "dx", "dy" }) do
        if c[k] then
            c[k] = c[k] * scale
        end
    end
end
--- Construct one step, with an owned previous potential and unchanged bulk IDs.
-- Sources are sampled at the new time. Motion is imposed, or staggered with
-- omega_(n+1)=(I*omega_n+dt*(torque_n-load))/(I+dt*damping), then theta advances
-- with omega_(n+1). This is partitioned first-order mechanics, not an implicit
-- monolithic electromechanical solve. Load is a signed opposing torque.
function T.new(state, o)
    o = o or {}
    local allowed = {
        dt = true,
        currents = true,
        source_scale = true,
        boundary_scale = true,
        cache = true,
        numerical = true,
        angle = true,
        offset_x = true,
        offset_y = true,
        motion = true,
        inertia = true,
        load = true,
        damping = true,
    }
    for k in pairs(o) do
        check(allowed[k], "unknown step option " .. tostring(k))
    end
    local d = S.validate(state)
    check(finite(o.dt) and o.dt > 0, "time step must be positive")
    check(finite(d.time + o.dt) and d.time + o.dt > d.time, "time step does not advance time")
    local mode = o.motion or "fixed"
    check(
        mode == "dynamic" or (o.inertia == nil and o.load == nil and o.damping == nil),
        "mechanical coefficients require dynamic motion"
    )
    check(mode == "fixed" or mode == "prescribed" or mode == "dynamic", "unknown motion mode")
    local motion = copy(d.motion)
    local M = require("luafemm")
    local m
    if mode ~= "fixed" then
        check(motion, "moving steps require an initial motion template")
        local angle = o.angle or motion.angle
        if mode == "dynamic" then
            check(
                o.angle == nil and o.offset_x == nil and o.offset_y == nil,
                "dynamic motion owns the rotation; offsets remain fixed"
            )
            local inertia, load, damping = o.inertia, o.load or 0, o.damping or 0
            check(
                finite(inertia)
                    and inertia > 0
                    and finite(load)
                    and finite(damping)
                    and damping >= 0,
                "dynamic motion requires positive inertia and nonnegative damping"
            )
            local old = S.restore(state)
            local spec = motion.template.options
            local r = require("luafemm-post").measure(old, {
                name = "_motion",
                register = false,
                quantity = "force",
                method = "band",
                origin_x = spec.center_x + motion.offset_x,
                origin_y = spec.center_y + motion.offset_y,
            })
            motion.speed = (inertia * motion.speed + o.dt * (r.torque - load))
                / (inertia + o.dt * damping)
            angle = motion.angle + math.deg(o.dt * motion.speed)
            motion.inertia, motion.load, motion.damping = inertia, load, damping
        else
            check(finite(angle), "invalid prescribed angle")
            motion.speed = math.rad(angle - motion.angle) / o.dt
        end
        motion.angle = angle
        motion.offset_x = o.offset_x or motion.offset_x
        motion.offset_y = o.offset_y or motion.offset_y
        m = require("luafemm-sliding").new(motion.template, {
            angle = angle,
            offset_x = motion.offset_x,
            offset_y = motion.offset_y,
            currents = o.currents,
            cache = o.cache,
            numerical = o.numerical,
        })
        check(#m.nodes == #d.nodes, "moving step changed material node identities")
        for i = 1, m.sliding.bulk_count do
            for j, id in ipairs(m.triangles[i].ids) do
                check(id == d.elements[i][j], "bulk connectivity changed")
            end
        end
    else
        check(
            o.angle == nil and o.offset_x == nil and o.offset_y == nil,
            "position overrides require motion=prescribed"
        )
        m = S.model(state)
        if motion then
            motion.speed = 0
        end
    end
    local source_scale = o.source_scale or (d.transient and d.transient.source_scale) or 1
    local boundary_scale = o.boundary_scale or (d.transient and d.transient.boundary_scale) or 1
    check(finite(source_scale) and finite(boundary_scale), "invalid source/boundary scale")
    for i, r in ipairs(m.regions) do
        r.current = d.sources.currents[i] * source_scale
    end
    m.boundaries = copy(d.sources.boundaries)
    for _, c in pairs(m.boundaries) do
        scale_boundary(c, boundary_scale)
    end
    -- Moving templates forbid attached non-rigid frontiers; fixed frontiers
    -- retain their geometry while affine loads scale relative to initial data.
    if mode == "fixed" then
        m.frontiers = copy(d.sources.frontiers)
    end
    for _, f in pairs(m.frontiers or {}) do
        if f.condition then
            scale_boundary(f.condition, boundary_scale)
        end
    end
    check(type(o.currents or {}) == "table", "currents must be a table")
    for name, c in pairs(m.circuits or {}) do
        c.current = d.model.circuits[name].current
    end
    for name, v in pairs(o.currents or {}) do
        check(
            m.circuits and m.circuits[name] and finite(v),
            "invalid circuit current " .. tostring(name)
        )
        m.circuits[name].current = v
    end
    local nodes, elements, boundary = {}, {}, copy(m.exterior_nodes)
    for i, p in ipairs(m.nodes) do
        nodes[i] = { p[1], p[2] }
    end
    for i, e in ipairs(m.triangles) do
        elements[i] = copy(e.ids)
    end
    m.transient, m.stats, m.solution_signature = nil, nil, nil
    require("luafemm-frontiers").compile(m)
    if m.topology then
        for _, f in ipairs(m.topology.faces) do
            if f.region then
                f.current = m.regions[f.region].current
            end
        end
        require("luafemm-static").excite(m, m.topology)
    end
    M.prepare_mesh(m, nodes, elements, boundary)
    m.transient = {
        dt = o.dt,
        time = d.time + o.dt,
        previous = copy(d.A),
        linkage = copy(d.linkage),
        energy = d.energy,
        source_scale = source_scale,
        boundary_scale = boundary_scale,
        motion = mode,
    }
    m.time_motion, m.time_sources = motion, copy(d.sources)
    m.time_mesh = { nodes = copy(nodes), elements = copy(elements) }
    for _, p in ipairs(m.time_mesh.nodes) do
        p.dof, p.dof_sign = nil, nil
    end
    for i, p in ipairs(m.nodes) do
        if p.dof then
            m.A[i] = d.A[i]
        end
    end
    m.options.cache = copy(o.cache or { mode = "off" })
    require("luafemm-cache").configure(m)
    require("luafemm-linear").override(m, o.numerical)
    T.eligible(m)
    T.prepare(m)
    m.transient_contract = require("luafemm-cache").signature(m)
    if m.sliding then
        m.sliding_contract = m.transient_contract
    end
    return m
end
function T.validate(m)
    T.eligible(m)
    local t = m.transient
    check(
        t and finite(t.dt) and t.dt > 0 and #t.previous == #m.nodes,
        "invalid previous time potential"
    )
    check(
        m.transient_contract == require("luafemm-cache").signature(m),
        "step declarations changed; derive a new step"
    )
    for i, p in ipairs(m.nodes) do
        local old = m.time_mesh.nodes[i]
        check(old and old[1] == p[1] and old[2] == p[2], "time mesh changed")
    end
    for i, e in ipairs(m.triangles) do
        for j, id in ipairs(e.ids) do
            check(m.time_mesh.elements[i][j] == id, "time connectivity changed")
        end
    end
end
--- Identify effective resolved conductors and constraint quadrature vectors.
function T.prepare(m)
    if not m.transient then
        return
    end
    check(
        not m.gauge_node and not (m.gauge_nodes and next(m.gauge_nodes)),
        "time diffusion requires an explicit potential reference; automatic static gauges are unsupported"
    )
    local groups = {}
    for index, e in ipairs(m.triangles) do
        local region = e.face and m.topology.faces[e.face].region
        if not region then
            local x, y = 0, 0
            for _, id in ipairs(e.ids) do
                x = x + m.nodes[id][1] / 3
                y = y + m.nodes[id][2] / 3
            end
            for i, r in ipairs(m.regions) do
                if
                    require("luafemm-geometry").contains(r, x / m.options.unit, y / m.options.unit)
                then
                    region = i
                end
            end
        end
        local r = region and m.regions[region]
        e.time_region = region
        e.sigma = r and (r.conductor == "stranded" or r.circuit) and 0 or e.p.conductivity
        check(
            not m.sliding or index <= m.sliding.bulk_count or e.sigma == 0,
            "moving gap must be nonconducting"
        )
        e.time_group = nil
        if r and r.current_constraint == "zero" then
            check(e.sigma > 0, "zero-current region must be conducting")
            local g = groups[region] or { weights = {}, sum = 0 }
            groups[region] = g
            e.time_group = g
            for _, id in ipairs(e.ids) do
                local w = e.sigma * e.area / 3
                g.weights[id] = (g.weights[id] or 0) + w
                g.sum = g.sum + w
            end
        end
    end
    for _, g in pairs(groups) do
        g.ids = {}
        for id in pairs(g.weights) do
            g.ids[#g.ids + 1] = id
        end
        table.sort(g.ids)
        check(#g.ids <= 2048, "isolated conductor exceeds 2048 nodes")
    end
    m.time_groups = groups
end
local function means(m, A)
    local values = {}
    for id, g in pairs(m.time_groups) do
        local sum = 0
        for _, n in ipairs(g.ids) do
            sum = sum + g.weights[n] * (A[n] - m.transient.previous[n])
        end
        values[id] = sum / g.sum
    end
    return values
end
--- Add consistent diffusion mass and its constrained-conductor Schur terms.
function T.assemble(m, A, R, K)
    if not m.transient then
        return
    end
    if not m.time_groups then
        T.prepare(m)
    end
    local scale = mu0 / m.transient.dt
    for _, e in ipairs(m.triangles) do
        if e.sigma > 0 then
            for _, q in ipairs(require("luafemm-harmonic").quadrature(m, e)) do
                local delta = 0
                for i, id in ipairs(e.ids) do
                    delta = delta + q.n[i] * (A[id] - m.transient.previous[id])
                end
                for i, id in ipairs(e.ids) do
                    local p = m.nodes[id]
                    local row = p.dof
                    if row then
                        local w = scale * e.sigma * q.weight * q.n[i] * (p.dof_sign or 1)
                        R[row] = R[row] + w * delta
                        if K then
                            for j, jd in ipairs(e.ids) do
                                local p2 = m.nodes[jd]
                                if p2.dof then
                                    K[row][p2.dof] = (K[row][p2.dof] or 0)
                                        + w * q.n[j] * (p2.dof_sign or 1)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    local average = means(m, A)
    for id, g in pairs(m.time_groups) do
        for _, i in ipairs(g.ids) do
            local p = m.nodes[i]
            local row = p.dof
            if row then
                local w = scale * g.weights[i] * (p.dof_sign or 1)
                R[row] = R[row] - w * average[id]
                if K then
                    for _, j in ipairs(g.ids) do
                        local p2 = m.nodes[j]
                        if p2.dof then
                            K[row][p2.dof] = (K[row][p2.dof] or 0)
                                - w * g.weights[j] * (p2.dof_sign or 1) / g.sum
                        end
                    end
                end
            end
        end
    end
end
--- Pointwise instantaneous J, induced J, and material derivative of A.
-- Coordinates are SI; moving bulk IDs implement the material derivative,
-- including the planar motional v cross B term without a separate source.
function T.field(m, e, x, y)
    local p = m.nodes[e.ids[1]]
    local b = e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2])
    local c = e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
    local n = { 1 - b - c, b, c }
    local derivative = 0
    for i, id in ipairs(e.ids) do
        derivative = derivative + n[i] * (m.A[id] - m.transient.previous[id]) / m.transient.dt
    end
    local electric = e.time_group and e.time_group.electric or 0
    local induced = e.sigma * (electric - derivative)
    local source = e.current
    return { J = source + induced, Jeddy = induced, dAdt = derivative, E = electric - derivative }
end
function T.sample(m, x, y)
    require("luafemm").solve(m)
    x, y = x * m.options.unit, y * m.options.unit
    for _, e in ipairs(m.triangles) do
        local p = m.nodes[e.ids[1]]
        local b = e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2])
        local c = e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
        if b >= -1e-10 and c >= -1e-10 and b + c <= 1 + 1e-10 then
            return T.field(m, e, x, y)
        end
    end
    error("luafemm transient: sample outside domain")
end
--- Refresh conductor fields and instantaneous (not cycle-average) losses.
function T.update(m)
    if not m.transient then
        return
    end
    T.prepare(m)
    local average = means(m, m.A)
    for id, g in pairs(m.time_groups) do
        g.electric = average[id] / m.transient.dt
    end
    local loss, current = 0, 0
    local axi = m.options.problem_type == "axisymmetric"
    local depth = (m.options.depth or 1) * m.options.unit
    for _, e in ipairs(m.triangles) do
        e.time_loss, e.time_current = 0, 0
        for _, q in ipairs(require("luafemm-harmonic").quadrature(m, e)) do
            local f = T.field(m, e, q.r, q.z)
            local volume = q.weight * (axi and 2 * math.pi or depth)
            local r = e.time_region and m.regions[e.time_region]
            local sigma = e.p.conductivity * (r and r.fill or 1)
            e.time_loss = e.time_loss + (sigma > 0 and volume * f.J * f.J / sigma or 0)
            e.time_current = e.time_current + q.weight * f.J / (axi and q.r or 1)
        end
        loss, current = loss + e.time_loss, current + e.time_current
    end
    m.stats.joule_loss, m.stats.net_current = loss, current
    m.stats.time, m.stats.time_step = m.transient.time, m.transient.dt
    if m.time_motion then
        m.stats.angle, m.stats.speed = m.time_motion.angle, m.time_motion.speed
    end
end
--- Add induced terminal voltage to the existing current-fed winding result.
function T.circuit(m, name, r)
    if not m.transient then
        return r
    end
    r.inductance = nil
    local old = m.transient.linkage[name]
    if old and r.linkage then
        r.induced_voltage = (r.linkage - old) / m.transient.dt
        if r.resistance then
            r.voltage = r.resistance * r.current + r.induced_voltage
        end
        r.power = r.voltage and r.voltage * r.current
    end
    return r
end
return T
