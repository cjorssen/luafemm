-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Single-frequency magnetic diffusion and bounded fundamental approximations with peak exp(+i omega t) phasors.
-- Uncontrolled solids have zero imposed voltage; constrained conductors and
-- winding terminals use the coupled network in luafemm-ac-circuits. Stranded
-- regions carry uniform currents without resolved eddy currents.
-- An operating response uses the captured DC differential tensor for both
-- phasor components. The explicitly selected fundamental law instead solves
-- a lossless, real nonlinear problem at AC amplitude.
-- @module luafemm-harmonic
local H = {}
local mu0 = 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm harmonic: " .. s)
end
local function finite(v)
    return type(v) == "number" and v == v and math.abs(v) < math.huge
end
--- Reject unsupported physics before a cache can return an apparently valid field.
-- @tparam table m Model with physical declarations.
function H.validate(m)
    if m.operating then
        require("luafemm-operating").validate(m)
    end
    if m.options.problem_type == "axial" then
        require("luafemm-axial").validate(m)
    end
    local o = m.options
    check(finite(o.frequency) and o.frequency > 0, "frequency must be positive")
    check(o.field_model == "open", "ideal circuits are static only")
    for _, frontier in pairs(m.frontiers or {}) do
        check(
            not frontier.condition or not frontier.condition.open_radius,
            "harmonic circular open boundaries are not validated yet"
        )
    end
    check(not o.source and not o.boundary, "harmonic source/boundary callbacks are unsupported")
    for _, p in pairs(m.materials) do
        check(
            not p.bh
                or m.operating
                or p.harmonic_law == "fundamental"
                or p.harmonic_law == "equivalent",
            "nonlinear harmonic material laws require a static operating point or harmonic law=fundamental"
        )
        check(p.hc == 0 or m.operating, "permanent-magnet DC bias requires a separate static solve")
        require("luafemm-losses").validate(m, p)
    end
    require("luafemm-windings").validate(m)
    require("luafemm-fundamental").validate(m)
    require("luafemm-equivalent").validate(m)
    require("luafemm-ac-circuits").validate(m)
    for _, r in ipairs(m.regions) do
        check(
            not r.circuit or (r.current_phase or 0) == 0,
            "set circuit phase on the circuit, not its legs"
        )
    end
end
--- Prepare each triangle's source phase and effective eddy conductivity.
-- @tparam table m Meshed model; derived element data are updated in place.
function H.prepare(m)
    local network = require("luafemm-ac-circuits")
    local coupled = network.active(m)
    local geo = require("luafemm-geometry")
    local losses = require("luafemm-losses")
    require("luafemm-windings").prepare(m)
    local tensors = {}
    for _, p in pairs(m.materials) do
        if losses.active(p) and p.harmonic_law ~= "equivalent" then
            tensors[p] = losses.material(p, m.options.frequency)
        end
    end
    for index, e in ipairs(m.triangles) do
        local region
        e.region_index = nil
        if e.face then
            e.region_index = m.topology.faces[e.face].region
            region = m.regions[e.region_index]
        else
            local x, y = 0, 0
            for _, id in ipairs(e.ids) do
                x = x + m.nodes[id][1] / (3 * m.options.unit)
                y = y + m.nodes[id][2] / (3 * m.options.unit)
            end
            for id, r in ipairs(m.regions) do
                if geo.contains(r, x, y) then
                    region = r
                    e.region_index = id
                end
            end
        end
        local c = region and region.circuit and m.circuits[region.circuit]
        if coupled and c then
            e.current = 0
            for _, q in ipairs(e.quadrature or {}) do
                q.current = 0
            end
        end
        local phase = math.rad(c and c.phase or region and region.current_phase or 0)
        e.jre, e.jim = e.current * math.cos(phase), e.current * math.sin(phase)
        e.stranded = region
            and (region.conductor == "stranded" or c and region.conductor ~= "solid")
        e.wire = m.wire_regions[e.region_index]
        e.ac_material = e.wire or tensors[e.p]
        if e.kelvin_factor then
            local t, k = e.p.tensor, e.kelvin_factor
            e.ac_material = {
                real = { k * t[1], k * t[2], k * t[3] },
                imaginary = { 0, 0, 0 },
                hysteresis = { 0, 0, 0 },
                lamination = { 0, 0, 0 },
            }
        end
        if e.p.harmonic_law == "equivalent" then
            e.ac_material = m.equivalent_tensors and m.equivalent_tensors[index]
                or require("luafemm-equivalent").element(m, e)
        end
        e.sigma = (e.stranded or e.p.lamination == "stacked") and 0 or e.p.conductivity
        if e.p.lamination == "stacked" and region then
            check(
                not c
                    and region.current == 0
                    and region.current_constraint ~= "zero"
                    and region.conductor ~= "stranded",
                "homogenised sheets cannot carry imposed or constrained transport currents"
            )
        end
        e.loss_sigma = e.p.conductivity * (region and region.fill or 1)
    end
    if coupled then
        network.prepare(m)
        if m.stats and m.A_im then
            network.recover(m)
        end
    else
        m.ac_network = nil
    end
end
-- The planar three-point rule integrates the consistent mass exactly.
local planar = { { 2 / 3, 1 / 6, 1 / 6 }, { 1 / 6, 2 / 3, 1 / 6 }, { 1 / 6, 1 / 6, 2 / 3 } }
--- SI quadrature shared by diffusion assembly and cycle-average observables.
-- @tparam table m Model.
-- @tparam table e Prepared triangle.
-- @treturn table Quadrature points with basis values and SI weighted measures.
function H.quadrature(m, e)
    if e.quadrature then
        return e.quadrature
    end
    local out = {}
    for _, p in ipairs(planar) do
        local x, y = 0, 0
        for i, id in ipairs(e.ids) do
            x = x + p[i] * m.nodes[id][1]
            y = y + p[i] * m.nodes[id][2]
        end
        out[#out + 1] = { n = p, r = x, z = y, weight = e.area / 3 }
    end
    return out
end
--- Refresh phasor norms and total section-current statistics.
-- Element display fields remain the phase-zero snapshot; peak is the largest
-- centroid complex vector norm, not the maximum instantaneous rotating field.
-- @tparam table m Solved model with both potential arrays and statistics.
function H.statistics(m)
    local peak, ir, ii = 0, 0, 0
    for _, e in ipairs(m.triangles) do
        local x, y = 0, 0
        for _, id in ipairs(e.ids) do
            x = x + m.nodes[id][1] / 3
            y = y + m.nodes[id][2] / 3
        end
        local f = H.field(m, e, x, y)
        peak = math.max(peak, math.sqrt(f.Bx[1] ^ 2 + f.Bx[2] ^ 2 + f.By[1] ^ 2 + f.By[2] ^ 2))
        for _, q in ipairs(H.quadrature(m, e)) do
            f = H.field(m, e, q.r, q.z)
            local area = q.weight / (e.quadrature and q.r or 1)
            ir, ii = ir + area * f.J[1], ii + area * f.J[2]
        end
    end
    m.stats.peak, m.stats.net_current, m.stats.net_current_imaginary = peak, ir, ii
end

--- Solve the complex diffusion system using the static real stiffness operator.
-- Real Dirichlet/Neumann/Robin data are peak, phase-zero boundary phasors.
-- @tparam table m Prepared model.
-- @tparam function assemble Shared real residual/tangent assembler.
-- @treturn table The solved model; failure raises before cache storage.
function H.solve(m, assemble)
    H.validate(m)
    if require("luafemm-equivalent").active(m) and not m.equivalent_tensors then
        return require("luafemm-equivalent").solve(m, assemble)
    end
    if m.options.problem_type == "axial" then
        return require("luafemm-axial").solve(m)
    end
    H.prepare(m)
    check(
        not m.gauge_nodes or next(m.gauge_nodes) == nil,
        "harmonic problems require anchored boundaries, not an automatic gauge"
    )
    if require("luafemm-fundamental").active(m) then
        m.A_im = {}
        for id in ipairs(m.nodes) do
            m.A_im[id] = 0
        end
        return require("luafemm").solve_real(m)
    end
    local start = os.clock()
    local R, K = assemble(m, m.A, true)
    local n = #m.free
    local br, bi, Ki = {}, {}, {}
    for i = 1, n do
        br[i], bi[i], Ki[i] = -R[i], 0, {}
    end
    local omega = 2 * math.pi * m.options.frequency
    for _, e in ipairs(m.triangles) do
        for _, q in ipairs(H.quadrature(m, e)) do
            for i, id in ipairs(e.ids) do
                local row = m.nodes[id].dof
                local si = m.nodes[id].dof_sign or 1
                if row then
                    local w = mu0 * q.weight * si * q.n[i]
                    br[row] = br[row] + w * (e.jre - e.current)
                    bi[row] = bi[row] + w * e.jim
                    for j, jd in ipairs(e.ids) do
                        local v = w * omega * e.sigma * q.n[j]
                        if e.ac_material then
                            local t = e.ac_material.imaginary
                            local ix, iy = q.br and q.br[i] or e.gy[i], q.bz and q.bz[i] or -e.gx[i]
                            local jx, jy = q.br and q.br[j] or e.gy[j], q.bz and q.bz[j] or -e.gx[j]
                            v = v
                                + q.weight
                                    * si
                                    * (ix * (t[1] * jx + t[2] * jy) + iy * (t[2] * jx + t[3] * jy))
                        end
                        local col = m.nodes[jd].dof
                        if col then
                            Ki[row][col] = (Ki[row][col] or 0) + v * (m.nodes[jd].dof_sign or 1)
                        else
                            bi[row] = bi[row] - v * m.A[jd]
                        end
                    end
                end
            end
        end
    end
    local assembly_seconds = os.clock() - start
    local ar, ai, it, res, telemetry
    if m.ac_network then
        ar, ai, it, res, telemetry = require("luafemm-ac-circuits").solve(m, K, Ki, br, bi)
    else
        ar, ai, it, res, telemetry = require("luafemm-linear").complex(m, K, Ki, br, bi)
    end
    m.A_im = {}
    for id, p in ipairs(m.nodes) do
        if p.dof then
            m.A[id] = (p.dof_sign or 1) * ar[p.dof]
            m.A_im[id] = (p.dof_sign or 1) * ai[p.dof]
        else
            m.A_im[id] = 0
        end
    end
    m.stats = {
        newton = 0,
        cg = 0,
        gmres = it,
        linear = telemetry.solves and telemetry
            or require("luafemm-linear").accumulate(nil, telemetry),
        assembly_seconds = assembly_seconds,
        residual = res,
        seconds = os.clock() - start,
        nodes = #m.nodes,
        triangles = #m.triangles,
    }
    require("luafemm").update_fields(m)
    require("luafemm-cache").store(m)
    return m
end
--- Return one real-valued view of a phasor; phase is omega*t in degrees.
-- @tparam number r Real phasor part.
-- @tparam number i Imaginary phasor part.
-- @tparam[opt=snapshot] string representation Scalar display convention.
-- @tparam[opt=0] number phase Temporal angle in degrees.
-- @treturn number Display value; zero-phasor argument is NaN.
function H.value(r, i, representation, phase)
    representation = representation or "snapshot"
    phase = phase or 0
    check(finite(phase), "phase must be finite")
    if representation == "real" then
        return r
    elseif representation == "imaginary" then
        return i
    elseif representation == "amplitude" then
        return math.sqrt(r * r + i * i)
    elseif representation == "rms" then
        return math.sqrt((r * r + i * i) / 2)
    elseif representation == "phase" then
        return r == 0 and i == 0 and 0 / 0 or math.deg(math.atan(i, r))
    elseif representation == "snapshot" then
        return r * math.cos(math.rad(phase)) - i * math.sin(math.rad(phase))
    end
    error("luafemm harmonic: unknown representation " .. tostring(representation))
end
--- Evaluate complex A, B, H, total J and induced J in a known triangle.
-- Returned values are pairs {real,imaginary}, in SI; J points out of the
-- planar section or in the positive azimuthal direction.
-- @tparam table m Solved model.
-- @tparam table e Known containing triangle.
-- @tparam number x SI horizontal coordinate or radius.
-- @tparam number y SI vertical or axial coordinate.
-- @treturn table A, Bx, By, Hx, Hy, J and Jeddy phasor pairs.
function H.field(m, e, x, y)
    if e.analytical_gap then
        local out = { A = {}, Bx = {}, By = {}, Hx = {}, Hy = {}, J = { 0, 0 }, Jeddy = { 0, 0 } }
        for k, A in ipairs({ m.A, m.A_im }) do
            out.Bx[k], out.By[k], out.A[k], out.Hx[k], out.Hy[k] =
                require("luafemm-airgap").field(m, A, x, y)
        end
        return out
    end
    if m.options.problem_type == "axial" then
        return require("luafemm-axial").field(m, e, x, y)
    end
    local p = m.nodes[e.ids[1]]
    local b = e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2])
    local c = e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
    local basis = { 1 - b - c, b, c }
    local r = { A = {}, Bx = {}, By = {}, Hx = {}, Hy = {} }
    for part, A in ipairs({ m.A, m.A_im }) do
        local a, bx, by = 0, 0, 0
        if m.options.problem_type == "axisymmetric" then
            bx, by, a = require("luafemm-axisymmetric").field(m, e, x, y, A)
        else
            for i, id in ipairs(e.ids) do
                a = a + basis[i] * A[id]
                bx = bx + e.gy[i] * A[id]
                by = by - e.gx[i] * A[id]
            end
        end
        local hx, hy = require("luafemm-static").element_law(m, e, bx, by, x, y)
        r.A[part], r.Bx[part], r.By[part], r.Hx[part], r.Hy[part] = a, bx, by, hx / mu0, hy / mu0
    end
    if e.ac_material then
        local apply = require("luafemm-losses").apply
        local xr, yr = apply(e.ac_material.imaginary, r.Bx[1], r.By[1])
        local xi, yi = apply(e.ac_material.imaginary, r.Bx[2], r.By[2])
        r.Hx[1], r.Hy[1] = r.Hx[1] - xi / mu0, r.Hy[1] - yi / mu0
        r.Hx[2], r.Hy[2] = r.Hx[2] + xr / mu0, r.Hy[2] + yr / mu0
    end
    local w = 2 * math.pi * m.options.frequency * e.sigma
    r.Jeddy = { w * r.A[2], -w * r.A[1] }
    r.J = { e.jre + r.Jeddy[1], e.jim + r.Jeddy[2] }
    if m.ac_network then
        local j, eddy = require("luafemm-ac-circuits").field(m, e, x, r.A)
        if j then
            r.J, r.Jeddy = j, eddy
        end
    end
    return r
end
--- Sample a solved harmonic field in model coordinates.
-- @tparam table m Solved harmonic model.
-- @tparam number x Model horizontal coordinate or radius.
-- @tparam number y Model vertical or axial coordinate.
-- @treturn table SI phasor pairs, using the first containing element at interfaces.
function H.sample(m, x, y)
    check(m.stats and m.A_im, "solve a harmonic model before sampling phasors")
    check(finite(x) and finite(y), "coordinates must be finite")
    x, y = x * m.options.unit, y * m.options.unit
    if require("luafemm-airgap").contains(m, x, y) then
        return H.field(m, { analytical_gap = true }, x, y)
    end
    if m.sliding and m.sliding.sector then
        local V = require("luafemm-sector-view")
        local e, xx, yy, angle, sign = V.locate(m, x, y)
        check(e, "sample outside physical sector domain")
        local f = H.field(m, e, xx, yy)
        for i = 1, 2 do
            f.Bx[i], f.By[i] = V.vector(f.Bx[i], f.By[i], angle, sign)
            f.Hx[i], f.Hy[i] = V.vector(f.Hx[i], f.Hy[i], angle, sign)
            for _, key in ipairs({ "A", "J", "Jeddy" }) do
                if f[key] then
                    f[key][i] = sign * f[key][i]
                end
            end
        end
        return f
    end
    for _, e in ipairs(m.triangles) do
        local p = m.nodes[e.ids[1]]
        local b = e.gx[2] * (x - p[1]) + e.gy[2] * (y - p[2])
        local c = e.gx[3] * (x - p[1]) + e.gy[3] * (y - p[2])
        if b >= -1e-10 and c >= -1e-10 and b + c <= 1 + 1e-10 then
            return H.field(m, e, x, y)
        end
    end
    error("luafemm: sample outside domain")
end
--- Create a lightweight static view for instantaneous field-line rendering.
-- The original phasors and computation signature are never changed.
-- @tparam table m Solved harmonic model.
-- @tparam[opt=0] number phase Temporal angle in degrees.
-- @treturn table Static rendering view sharing immutable geometry.
function H.snapshot(m, phase)
    check(
        m.options.problem_type ~= "axial",
        "axial magnetic field is out of plane; use current views"
    )
    local s = {}
    for k, v in pairs(m) do
        s[k] = v
    end
    s.options = {}
    for k, v in pairs(m.options) do
        s.options[k] = v
    end
    s.phasor_model, s.phasor_phase = m, phase or 0
    s.options.frequency = 0
    s.A, s.A_im, s.triangles = {}, nil, {}
    for i, a in ipairs(m.A) do
        s.A[i] = H.value(a, m.A_im[i], "snapshot", phase)
    end
    for i, e in ipairs(m.triangles) do
        local t = {}
        for k, v in pairs(e) do
            t[k] = v
        end
        local x, y = 0, 0
        for _, id in ipairs(e.ids) do
            x = x + m.nodes[id][1] / 3
            y = y + m.nodes[id][2] / 3
        end
        local f = H.field(m, e, x, y)
        t.bx, t.by =
            H.value(f.Bx[1], f.Bx[2], "snapshot", phase),
            H.value(f.By[1], f.By[2], "snapshot", phase)
        t.hx, t.hy =
            H.value(f.Hx[1], f.Hx[2], "snapshot", phase),
            H.value(f.Hy[1], f.Hy[2], "snapshot", phase)
        s.triangles[i] = t
    end
    return s
end
--- Integrate cycle-average energy and Ohmic losses on an element selection.
-- Conductive solids use |J_total|^2/(2 sigma); stranded homogenised copper
-- uses its fill-adjusted conductivity and prescribed current density.
-- @tparam table m Solved model with positive planar depth when applicable.
-- @tparam table chosen Boolean map keyed by triangle index.
-- @treturn table SI regional cycle-average quantities.
function H.integrals(m, chosen)
    if m.options.problem_type == "axial" then
        return require("luafemm-axial").integrals(m, chosen)
    end
    local axi = m.options.problem_type == "axisymmetric"
    local r = {
        area = 0,
        volume = 0,
        energy = 0,
        coenergy = 0,
        wire_transport_loss = 0,
        proximity_loss = 0,
        wire_energy = 0,
        joule_loss = 0,
        eddy_loss = 0,
        hysteresis_loss = 0,
        lamination_loss = 0,
        current_real = 0,
        current_imaginary = 0,
    }
    for i, e in ipairs(m.triangles) do
        if chosen[i] then
            r.area = r.area + e.area
            for _, q in ipairs(H.quadrature(m, e)) do
                local f = H.field(m, e, q.r, q.z)
                local volume = q.weight * (axi and 2 * math.pi or m.options.depth * m.options.unit)
                local area = q.weight / (axi and q.r or 1)
                r.volume = r.volume + volume
                r.current_real = r.current_real + area * f.J[1]
                r.current_imaginary = r.current_imaginary + area * f.J[2]
                local energy = 0
                for part = 1, 2 do
                    energy = energy + (f.Bx[part] * f.Hx[part] + f.By[part] * f.Hy[part]) / 4
                end
                local coenergy = energy
                if e.p.harmonic_law == "fundamental" then
                    local b = math.sqrt(f.Bx[1] ^ 2 + f.By[1] ^ 2)
                    local h, _, stored = require("luafemm-fundamental").scalar(e.p, b)
                    energy, coenergy = stored, b * h / 2 - stored
                end
                if e.analytical_gap then
                    energy, coenergy = 0, 0
                end
                local transport, proximity, local_energy =
                    require("luafemm-windings").density(m, e, f)
                r.wire_transport_loss = r.wire_transport_loss + volume * transport
                r.proximity_loss = r.proximity_loss + volume * proximity
                r.wire_energy = r.wire_energy + volume * local_energy
                r.energy = r.energy + volume * (energy + local_energy)
                r.coenergy = r.coenergy + volume * (coenergy + local_energy)
                r.joule_loss = r.joule_loss + volume * (transport + proximity)
                local intrinsic, sheet = require("luafemm-losses").density(m, e, f)
                r.hysteresis_loss = r.hysteresis_loss + volume * intrinsic
                r.lamination_loss = r.lamination_loss + volume * sheet
                local sigma = e.stranded and e.loss_sigma or e.sigma
                if sigma > 0 and not e.wire then
                    r.joule_loss = r.joule_loss + volume * (f.J[1] ^ 2 + f.J[2] ^ 2) / (2 * sigma)
                end
                if e.sigma > 0 then
                    r.eddy_loss = r.eddy_loss
                        + volume * (f.Jeddy[1] ^ 2 + f.Jeddy[2] ^ 2) / (2 * e.sigma)
                end
            end
        end
    end
    local gap_energy = require("luafemm-airgap").selected_energy(m, chosen)
    r.energy, r.coenergy = r.energy + gap_energy, r.coenergy + gap_energy
    r.magnetic_loss = r.hysteresis_loss + r.lamination_loss
    r.total_loss = r.joule_loss + r.magnetic_loss
    return r
end
--- Complex terminal quantities for harmonic circuits and optional branches.
-- Linkage is gauge independent only for closed planar windings. Terminal power
-- includes mutual exchange with other driven circuits; it is not a local loss.
-- @tparam table m Solved harmonic model with an unchanged signature.
-- @tparam string name Declared circuit name.
-- @tparam[opt] string winding Optional branch name; empty selects the whole port.
-- @treturn table SI complex terminal quantities; unavailable fields are absent.
function H.circuit(m, name, winding)
    check(m.stats and m.A_im, "solve before measuring a harmonic circuit")
    check(m.solution_signature == require("luafemm-cache").signature(m), "stale circuit solution")
    if m.ac_network then
        return require("luafemm-ac-circuits").result(m, name, winding)
    end
    if winding and winding ~= "" then
        check(
            m.options.problem_type == "axisymmetric" or m.options.depth,
            "set depth for circuit linkage"
        )
        require("luafemm-ac-circuits").prepare(m)
        require("luafemm-ac-circuits").recover(m)
        return require("luafemm-ac-circuits").result(m, name, winding)
    end
    local c = m.static_circuits and m.static_circuits[name]
    check(c, "unknown circuit " .. tostring(name))
    local axi = m.options.problem_type == "axisymmetric"
    check(axi or m.options.depth, "set depth for circuit linkage")
    local phase = math.rad(m.circuits[name].phase or 0)
    local ir, ii = c.current * math.cos(phase), c.current * math.sin(phase)
    local r = { current_real = ir, current_imaginary = ii, resistance = c.resistance }
    for _, g in pairs(c.windings) do
        if not axi and math.abs(g.turns) > 1e-10 then
            return require("luafemm-ac-circuits").decorate(r)
        end
    end
    local lr, li = 0, 0
    for _, e in ipairs(m.triangles) do
        local face = m.topology.faces[e.face]
        local region = face.region and m.regions[face.region]
        if region and region.circuit == name then
            for _, q in ipairs(H.quadrature(m, e)) do
                local ar, ai = 0, 0
                for i, id in ipairs(e.ids) do
                    ar = ar + q.n[i] * m.A[id]
                    ai = ai + q.n[i] * m.A_im[id]
                end
                local w = region.turns
                    * q.weight
                    * (axi and 2 * math.pi or m.options.depth * m.options.unit)
                    / c.areas[face.region]
                lr, li = lr + w * ar, li + w * ai
            end
        end
    end
    r.linkage_real, r.linkage_imaginary = lr, li
    -- Only a single driven port defines this self-response ratio. With several
    -- simultaneous excitations, report linkages and leave matrix assembly to
    -- the caller. Complex mutual entries follow from unit-current port solves.
    local d = ir * ir + ii * ii
    if m.operating and d > 0 then
        local isolated = true
        for other, circuit in pairs(m.circuits) do
            isolated = isolated and (other == name or circuit.current == 0)
        end
        if isolated then
            r.incremental_inductance_real = (lr * ir + li * ii) / d
            r.incremental_inductance_imaginary = (li * ir - lr * ii) / d
        end
    end
    if c.resistance then
        local omega = 2 * math.pi * m.options.frequency
        local vr, vi = c.resistance * ir - omega * li, c.resistance * ii + omega * lr
        r.voltage_real, r.voltage_imaginary = vr, vi
        r.power, r.reactive_power = (vr * ir + vi * ii) / 2, (vi * ir - vr * ii) / 2
        if d > 0 then
            r.impedance_real, r.impedance_imaginary =
                (vr * ir + vi * ii) / d, (vi * ir - vr * ii) / d
        end
    end
    return require("luafemm-ac-circuits").decorate(r)
end
return H
