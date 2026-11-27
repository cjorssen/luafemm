-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Fourier Dirichlet-to-Neumann operator for a source-free circular air gap.
-- P1 angular trace functions are projected by exact integrals. Each retained
-- Fourier mode solves the radial Laplace equation analytically in log(radius).
-- The Fourier construction is independent of FEMM's AGE stencil. Native AGE
-- choices dispatch to luafemm-age through the same assembly/sampling contract.
-- Bulk elements remain ordinary P1 triangles; the displayed band triangulation
-- carries no stiffness. All private operator/field caches are reconstructible.
local G = {}
local pi, mu0 = math.pi, 4 * math.pi * 1e-7
local function check(v, s)
    assert(v, "luafemm analytical gap: " .. s)
end
function G.active(m)
    return m.native_age or m.sliding and m.sliding.options.method ~= "p1"
end
local function angle(m, id)
    local o = m.sliding.options
    local p = m.nodes[id]
    local u = m.options.unit
    return math.atan(p[2] - o.center_y * u, p[1] - o.center_x * u) % (2 * pi)
end
-- Return coefficient rows: mean, then cosine/sine for each positive order.
local function trace(m, ids, modes, inner)
    local ring = {}
    local o = m.sliding.options
    local sectors = o.sectors or 1
    local anti = o.periodicity == "antiperiodic"
    local shift = m.sliding.sector and inner and math.rad(m.sliding.angle) or 0
    for sector = 0, sectors - 1 do
        for j, id in ipairs(ids) do
            -- The last sector endpoint is already identified with the first by
            -- strong periodic DOFs; retaining both would duplicate a trace hat.
            if sectors == 1 or j < #ids then
                ring[#ring + 1] = {
                    id,
                    (angle(m, id) + shift + 2 * pi * sector / sectors) % (2 * pi),
                    anti and (-1) ^ sector or 1,
                }
            end
        end
    end
    table.sort(ring, function(a, b)
        return a[2] < b[2]
    end)
    local rows = { {} }
    for n = 1, modes do
        rows[2 * n], rows[2 * n + 1] = {}, {}
    end
    for j, p in ipairs(ring) do
        local prev, nextp = ring[(j - 2) % #ring + 1], ring[j % #ring + 1]
        local t = p[2]
        local l = (t - prev[2]) % (2 * pi)
        local r = (nextp[2] - t) % (2 * pi)
        check(l > 1e-12 and r > 1e-12, "duplicate angular trace nodes")
        rows[1][p[1]] = (rows[1][p[1]] or 0) + p[3] * (l + r) / (4 * pi)
        for n = 1, modes do
            rows[2 * n][p[1]] = (rows[2 * n][p[1]] or 0)
                + p[3]
                    * ((math.cos(n * t) - math.cos(n * (t - l))) / l + (math.cos(n * t) - math.cos(
                        n * (t + r)
                    )) / r)
                    / (pi * n * n)
            rows[2 * n + 1][p[1]] = (rows[2 * n + 1][p[1]] or 0)
                + p[3]
                    * ((math.sin(n * t) - math.sin(n * (t - l))) / l + (math.sin(n * t) - math.sin(
                        n * (t + r)
                    )) / r)
                    / (pi * n * n)
        end
    end
    if anti then
        rows[1] = {}
    end
    if sectors > 1 then
        for n = 1, modes do
            if n % sectors ~= (anti and sectors / 2 or 0) then
                rows[2 * n], rows[2 * n + 1] = {}, {}
            end
        end
    end
    return rows
end
--- Rebuild the geometry-only operator; mark visual band triangles after restores.
function G.prepare(m)
    if not G.active(m) then
        return
    end
    for i, e in ipairs(m.triangles) do
        e.analytical_gap = m.sliding and i > m.sliding.bulk_count or nil
    end
    if m.airgap_operator then
        return m.airgap_operator
    end
    if m.native_age or m.sliding.options.method ~= "analytical" then
        local N = require("luafemm-age")
        m.native_age = m.native_age
            or { ages = N.records(m), edition = m.sliding.options.method, unit = m.options.unit }
        local op = N.operator(m.native_age.ages, m.native_age.edition)
        op.modes = 0
        for _, a in ipairs(m.native_age.ages) do
            op.modes = math.max(op.modes, math.ceil(a.header[9] / 2))
        end
        m.airgap_operator = op
        return op
    end
    local s, o = m.sliding, m.sliding.options
    local sectors = o.sectors or 1
    local count = (math.max(#s.inner, #s.outer) - (sectors > 1 and 1 or 0)) * sectors
    local modes = o.modes == 0 and math.ceil(count / 2) + 8 or o.modes
    check(
        modes >= math.ceil(count / 2) and modes <= 256,
        "modes must resolve both traces and be at most 256; use fewer segments or a coarser bulk mesh"
    )
    local ids = {}
    for _, ring in ipairs({ s.inner, s.outer }) do
        for _, id in ipairs(ring) do
            ids[#ids + 1] = id
        end
    end
    local op = {
        ids = ids,
        inner = trace(m, s.inner, modes, true),
        outer = trace(m, s.outer, modes, false),
        modes = modes,
        L = math.log(o.outer_radius / o.inner_radius),
        K = {},
    }
    for i = 1, #ids do
        op.K[i] = {}
        for j = 1, #ids do
            op.K[i][j] = 0
        end
    end
    local function term(a, b, scale)
        for i, id in ipairs(ids) do
            local v = scale / sectors * (a[id] or 0)
            if v ~= 0 then
                for j, jd in ipairs(ids) do
                    op.K[i][j] = op.K[i][j] + v * (b[jd] or 0)
                end
            end
        end
    end
    local d = {}
    for _, id in ipairs(ids) do
        d[id] = (op.inner[1][id] or 0) - (op.outer[1][id] or 0)
    end
    term(d, d, 2 * pi / op.L)
    for n = 1, modes do
        local z = math.exp(-n * op.L)
        local c = (1 + z * z) / (1 - z * z)
        local h = 2 * z / (1 - z * z)
        for _, k in ipairs({ 2 * n, 2 * n + 1 }) do
            local a, b = op.inner[k], op.outer[k]
            term(a, a, pi * n * c)
            term(b, b, pi * n * c)
            term(a, b, -pi * n * h)
            term(b, a, -pi * n * h)
        end
    end
    -- Enforce the exact constant null mode after floating-point accumulation.
    for i = 1, o.periodicity ~= "antiperiodic" and #ids or 0 do
        local sum = 0
        for j = 1, #ids do
            if j ~= i then
                sum = sum + op.K[i][j]
            end
        end
        op.K[i][i] = -sum
    end
    m.airgap_operator = op
    return op
end
--- Add the scaled weak operator to the residual and optional tangent matrix.
function G.assemble(m, A, R, K)
    local op = G.prepare(m)
    if not op then
        return
    end
    for i, id in ipairs(op.ids) do
        local p = m.nodes[id]
        if p.dof then
            local si = p.dof_sign or 1
            for j, jd in ipairs(op.ids) do
                local v = si * op.K[i][j]
                R[p.dof] = R[p.dof] + v * A[jd]
                local q = m.nodes[jd]
                if K and q.dof then
                    K[p.dof][q.dof] = (K[p.dof][q.dof] or 0) + v * (q.dof_sign or 1)
                end
            end
        end
    end
end
local function coefficients(m, A)
    local op = G.prepare(m)
    m.airgap_fields = m.airgap_fields or setmetatable({}, { __mode = "k" })
    if m.airgap_fields[A] then
        return op, m.airgap_fields[A]
    end
    local out = { {}, {} }
    for j, rows in ipairs({ op.inner, op.outer }) do
        for k, row in ipairs(rows) do
            local v = 0
            for id, w in pairs(row) do
                v = v + w * A[id]
            end
            out[j][k] = v
        end
    end
    m.airgap_fields[A] = out
    return op, out
end
--- Analytic Az and curl at SI coordinates; return Bx, By, Az, Hx, Hy.
function G.field(m, A, x, y)
    G.prepare(m)
    if m.native_age then
        return require("luafemm-age").field(m, A, x, y)
    end
    local op, a = coefficients(m, A)
    local o = m.sliding.options
    local u = m.options.unit
    x, y = x - o.center_x * u, y - o.center_y * u
    local r = math.sqrt(x * x + y * y)
    local t = math.log(r / (o.inner_radius * u))
    local theta = math.atan(y, x)
    local v = a[1][1] * (1 - t / op.L) + a[2][1] * t / op.L
    local dt, da = (a[2][1] - a[1][1]) / op.L, 0
    for n = 1, op.modes do
        local den = 1 - math.exp(-2 * n * op.L)
        local function sh(z)
            return math.exp(n * (z - op.L)) * (1 - math.exp(-2 * n * z)) / den
        end
        local function ch(z)
            return math.exp(n * (z - op.L)) * (1 + math.exp(-2 * n * z)) / den
        end
        local left, right = sh(op.L - t), sh(t)
        local dl, dr = -n * ch(op.L - t), n * ch(t)
        local c, s = math.cos(n * theta), math.sin(n * theta)
        local ac, bs =
            a[1][2 * n] * left + a[2][2 * n] * right,
            a[1][2 * n + 1] * left + a[2][2 * n + 1] * right
        v = v + ac * c + bs * s
        da = da + n * (-ac * s + bs * c)
        dt = dt
            + (a[1][2 * n] * dl + a[2][2 * n] * dr) * c
            + (a[1][2 * n + 1] * dl + a[2][2 * n + 1] * dr) * s
    end
    local br, bt = da / r, -dt / r
    local c, s = x / r, y / r
    local bx, by = br * c - bt * s, br * s + bt * c
    return bx, by, v, bx / mu0, by / mu0
end
function G.contains(m, x, y)
    if not G.active(m) then
        return false
    end
    G.prepare(m)
    if m.native_age then
        local u = m.native_age.unit
        for _, age in ipairs(m.native_age.ages) do
            local h = age.header
            local r = math.sqrt((x - h[7] * u) ^ 2 + (y - h[8] * u) ^ 2)
            if r >= h[4] * u and r <= h[5] * u then
                return true
            end
        end
        return false
    end
    local o = m.sliding.options
    local u = m.options.unit
    local r = math.sqrt((x - o.center_x * u) ^ 2 + (y - o.center_y * u) ^ 2)
    return r >= o.inner_radius * u and r <= o.outer_radius * u
end
--- Exact retained-mode energy (J) for a real potential vector.
function G.energy(m, A)
    local op = G.prepare(m)
    local v = 0
    for i, id in ipairs(op.ids) do
        for j, jd in ipairs(op.ids) do
            v = v + A[id] * op.K[i][j] * A[jd]
        end
    end
    return v * m.options.depth * m.options.unit / (2 * mu0)
end
--- Energy of a complete gap selection; modal energy is not element-additive.
function G.selected_energy(m, chosen)
    if not G.active(m) then
        return 0
    end
    local count, selected = 0, 0
    for i, e in ipairs(m.triangles) do
        if e.analytical_gap then
            count = count + 1
            if chosen[i] then
                selected = selected + 1
            end
        end
    end
    check(
        selected == 0 or selected == count,
        "energy selection must contain the complete analytical gap or none of it"
    )
    if selected == 0 then
        local complete = count == 0 and m.native_age
        for i in ipairs(m.triangles) do
            complete = complete and chosen[i]
        end
        if not complete then
            return 0
        end
    end
    local value = G.energy(m, m.A)
    if m.options.frequency > 0 then
        value = (value + G.energy(m, m.A_im)) / 2
    end
    return value
end
--- Integrate a traction callback on a circular contour inside the gap.
-- Callback receives SI x,y and outward normal line element nx,ny (metres).
function G.contour(m, callback)
    local op = G.prepare(m)
    if m.native_age then
        local u = m.native_age.unit
        for _, age in ipairs(m.native_age.ages) do
            local h = age.header
            local r = (h[4] + h[5]) * u / 2
            local n = 4 * h[9] + 4
            for i = 0, n - 1 do
                local t = math.rad(h[6]) * (i + 0.5) / n
                local c, s = math.cos(t), math.sin(t)
                local ds = r * math.rad(h[6]) / n
                callback(h[7] * u + r * c, h[8] * u + r * s, c * ds, s * ds)
            end
        end
        return
    end
    local o = m.sliding.options
    local u = m.options.unit
    local r = math.sqrt(o.inner_radius * o.outer_radius) * u
    local sectors = o.sectors or 1
    local count = 4 * op.modes + 4
    for i = 0, count - 1 do
        local t = math.rad(o.start_angle or 0) + 2 * pi * (i + 0.5) / (count * sectors)
        local c, s = math.cos(t), math.sin(t)
        local ds = r * 2 * pi / (count * sectors)
        callback(o.center_x * u + r * c, o.center_y * u + r * s, c * ds, s * ds)
    end
end
return G
