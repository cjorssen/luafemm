-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Native concentric FEMM air-gap elements, with explicit interpolation edition.
-- Independent construction: endpoint Hermite data -> four corner potentials ->
-- the rectangular 2/3 triangular, 1/3 bilinear energy blend. No expanded native
-- polynomial tables are embedded. ANS indices remain zero-based in raw records.
local N = {}
local pi = math.pi
local function check(v, s)
    assert(v, "luafemm AGE: " .. s)
end
local function shift_valid(value)
    return type(value) == "number" and value >= -1e-9 and value <= 1 + 1e-9
end
--- Four interpolation weights at t=1-shift, on the uniform stencil -1,0,1,2.
-- FEMM 4 uses central slopes; FEMM 5 also uses central second differences.
function N.weights(shift, edition)
    check(shift_valid(shift), "fractional shift must lie between zero and one")
    check(
        edition == "femm4" or edition == "femm5",
        "choose age interpolation=femm4 or femm5 explicitly"
    )
    local t = 1 - shift
    local h0, h1, d0, d1, s0, s1
    if edition == "femm4" then
        h0 = 2 * t ^ 3 - 3 * t ^ 2 + 1
        h1 = 1 - h0
        d0 = t ^ 3 - 2 * t ^ 2 + t
        d1 = t ^ 3 - t ^ 2
        s0, s1 = 0, 0
    else
        h0 = 1 - 10 * t ^ 3 + 15 * t ^ 4 - 6 * t ^ 5
        h1 = 1 - h0
        d0 = t - 6 * t ^ 3 + 8 * t ^ 4 - 3 * t ^ 5
        d1 = -4 * t ^ 3 + 7 * t ^ 4 - 3 * t ^ 5
        s0 = (t ^ 2 - 3 * t ^ 3 + 3 * t ^ 4 - t ^ 5) / 2
        s1 = (t ^ 3 - 2 * t ^ 4 + t ^ 5) / 2
    end
    return { -d0 / 2 + s0, h0 - d1 / 2 - 2 * s0 + s1, h1 + d0 / 2 + s0 - 2 * s1, d1 / 2 + s1 }
end
--- Dimensionless 10 by 10 symmetric stiffness (global solver scales by mu0).
-- The solver translates the corner origin by the outer fractional shift;
-- native postprocessing instead samples corners at the global angular grid.
function N.matrix(ci, co, aspect, edition)
    check(shift_valid(ci) and shift_valid(co), "fractional shifts must lie between zero and one")
    check(
        aspect > 0 and aspect < math.huge and 1 / aspect < math.huge,
        "invalid element aspect ratio"
    )
    if ci > co then
        ci, co = ci - co, 0
    else
        ci, co = 1 - co + ci, 1
    end
    local wi, wo = N.weights(ci, edition), N.weights(co, edition)
    local T = { {}, {}, {}, {} }
    for i = 1, 4 do
        for j = 1, 10 do
            T[i][j] = 0
        end
    end
    for j = 1, 4 do
        T[1][j], T[2][j + 1] = wi[j], wi[j]
        T[3][j + 5], T[4][j + 6] = wo[j], wo[j]
    end
    local k, l = aspect, 1 / aspect
    local d = 5 * (k + l) / 12
    local h, v, z = (-5 * k + l) / 12, (k - 5 * l) / 12, -(k + l) / 12
    local Q = { { d, h, v, z }, { h, d, z, v }, { v, z, d, h }, { z, v, h, d } }
    local out = {}
    for i = 1, 10 do
        out[i] = {}
        for j = 1, 10 do
            local sum = 0
            for a = 1, 4 do
                for b = 1, 4 do
                    sum = sum + T[a][i] * Q[a][b] * T[b][j]
                end
            end
            out[i][j] = sum
        end
    end
    return out
end
--- Native five-node stencils, including seam signs, for zero-based cell k.
function N.stencil(age, k)
    local n, q = age.header[9], age.quads
    local ids, w = {}, {}
    for side = 0, 1 do
        local entries = {
            { (k - 1) % n + 1, 1 },
            { k + 1, 1 },
            { k + 1, 3 },
            { k + 2, 3 },
            {
                k + 2 > n and 2 or k + 3,
                3,
            },
        }
        for j, p in ipairs(entries) do
            local col = p[2] + 4 * side
            ids[#ids + 1] = q[p[1]][col] + 1
            local sign = (age.header[1] == 1 and ((k == 0 and j == 1) or (k == n - 1 and j == 5)))
                    and -1
                or 1
            w[#w + 1] = sign * q[p[1]][col + 1]
        end
    end
    return ids, w
end
--- Assemble raw native records into an owned dense trace operator.
function N.operator(ages, edition)
    local cells = 0
    for _, age in ipairs(ages) do
        cells = cells + age.header[9]
        check(
            age.header[9] >= 3 and cells <= 2048,
            "native field import requires 3 to 2048 total gap cells"
        )
        local count = 360 / age.header[6]
        local rounded = math.floor(count + 0.5)
        check(
            math.abs(count - rounded) < 1e-7 and (age.header[1] == 0 or rounded % 2 == 0),
            "native AGE span must tile a full revolution with consistent periodicity"
        )
    end
    local ids, seen = {}, {}
    for _, age in ipairs(ages) do
        for _, q in ipairs(age.quads) do
            for j = 1, 7, 2 do
                local id = q[j] + 1
                if not seen[id] then
                    ids[#ids + 1] = id
                    seen[id] = #ids
                end
            end
        end
    end
    check(#ids <= 2048, "native trace operator exceeds 2048 nodes; use a coarser gap")
    local op = { ids = ids, K = {} }
    for i = 1, #ids do
        op.K[i] = {}
        for j = 1, #ids do
            op.K[i][j] = 0
        end
    end
    for _, age in ipairs(ages) do
        local h = age.header
        local aspect = 2 * (h[5] - h[4]) / (math.rad(h[6] / h[9]) * (h[5] + h[4]))
        local localK = N.matrix(h[10], h[11], aspect, edition)
        for k = 0, h[9] - 1 do
            local nn, w = N.stencil(age, k)
            for i = 1, 10 do
                for j = 1, 10 do
                    local a, b = seen[nn[i]], seen[nn[j]]
                    op.K[a][b] = op.K[a][b] + w[i] * w[j] * localK[i][j]
                end
            end
        end
    end
    return op
end
--- Native centre samples, in SI, before the native Fourier reconstruction.
function N.samples(age, A, edition, unit)
    local h = age.header
    local wi, wo = N.weights(h[10], edition), N.weights(h[11], edition)
    local r, dr = (h[4] + h[5]) * unit / 2, (h[5] - h[4]) * unit
    local dt = math.rad(h[6] / h[9])
    local out = {}
    for k = 0, h[9] - 1 do
        local nn, w = N.stencil(age, k)
        local a, b, c, d = 0, 0, 0, 0
        for j = 1, 4 do
            a = a + wi[j] * w[j] * A[nn[j]]
            b = b + wi[j] * w[j + 1] * A[nn[j + 1]]
            c = c + wo[j] * w[j + 5] * A[nn[j + 5]]
            d = d + wo[j] * w[j + 6] * A[nn[j + 6]]
        end
        out[k + 1] =
            { (b - a + d - c) / (2 * dt * r), (a + b - c - d) / (2 * dr), (a + b + c + d) / 4 }
    end
    return out
end
--- Generate representable full-circle AGE records from uniform physical rings.
-- Node coordinates already contain the rotor rotation: exported native angles
-- are therefore zero (the current geometry is the new reference position).
function N.records(m)
    local s, o = m.sliding, m.sliding.options
    local u = m.options.unit
    check(not s.sector, "native AGE generation currently requires a full circle")
    check(#s.inner == #s.outer, "native AGE traces require equal node counts; increase segments")
    local n = #s.inner
    local function ring(ids)
        local out = {}
        for _, id in ipairs(ids) do
            local p = m.nodes[id]
            local theta = math.atan(p[2] - o.center_y * u, p[1] - o.center_x * u) % (2 * pi)
            if 2 * pi - theta < 1e-10 then
                theta = 0
            end
            out[#out + 1] = { id = id, theta = theta }
        end
        table.sort(out, function(a, b)
            return a.theta < b.theta
        end)
        for j = 1, n do
            check(
                math.abs((out[j % n + 1].theta - out[j].theta) % (2 * pi) - 2 * pi / n) < 1e-7,
                "native AGE traces must be uniform; increase segments relative to mesh refinement"
            )
        end
        return out
    end
    local ri, ro = ring(s.inner), ring(s.outer)
    local age = {
        name = "luafemm sliding gap",
        header = {
            0,
            0,
            0,
            o.inner_radius,
            o.outer_radius,
            360,
            o.center_x,
            o.center_y,
            n,
            ri[1].theta * n / (2 * pi),
            ro[1].theta * n / (2 * pi),
        },
        quads = {},
        rings = {},
    }
    for j = 0, n do
        local a, b = j % n + 1, (j - 1) % n + 1
        age.quads[j + 1] = { ri[b].id - 1, 1, ri[a].id - 1, 1, ro[b].id - 1, 1, ro[a].id - 1, 1 }
        if j < n then
            age.rings[j + 1] = { ri[a].id - 1, ro[a].id - 1 }
        end
    end
    return { age }
end
--- Fourier reconstruction of native midgap samples, extended radially by Laplace.
-- At midgap this is the FEMM postprocessor's DFT; elsewhere it is a harmonic
-- continuation, not the bulk P1 interpolation. Samples own no physical state.
local function coefficients(m, A)
    m.airgap_fields = m.airgap_fields or setmetatable({}, { __mode = "k" })
    if m.airgap_fields[A] then
        return m.airgap_fields[A]
    end
    local result = {}
    for index, age in ipairs(m.native_age.ages) do
        local h = age.header
        local samples = N.samples(age, A, m.native_age.edition, m.native_age.unit)
        local n = h[9]
        local dt = math.rad(h[6] / n)
        local count = h[1] == 0 and math.floor(n / 2) + 1 or math.floor((n + 1) / 2)
        local rows = {}
        for j = 0, count - 1 do
            local order = h[1] == 0 and math.floor(360 / h[6] + 0.5) * j
                or math.floor(180 / h[6] + 0.5) * (2 * j + 1)
            local rc, rs, tc, ts, mean = 0, 0, 0, 0, 0
            local div = (order == 0 or h[1] == 0 and n % 2 == 0 and j == count - 1) and n or n / 2
            for k, s in ipairs(samples) do
                local c, t = math.cos((k - 0.5) * dt * order), math.sin((k - 0.5) * dt * order)
                rc, rs, tc, ts =
                    rc + s[1] * c / div,
                    rs + s[1] * t / div,
                    tc + s[2] * c / div,
                    ts + s[2] * t / div
                mean = mean + s[3] / n
            end
            rows[#rows + 1] = { order, rc, rs, tc, ts, mean }
        end
        result[index] = rows
    end
    m.airgap_fields[A] = result
    return result
end
function N.field(m, A, x, y)
    local data = m.native_age
    local u = data.unit
    local all = coefficients(m, A)
    for index, age in ipairs(data.ages) do
        local h = age.header
        local xx, yy = x - h[7] * u, y - h[8] * u
        local r = math.sqrt(xx * xx + yy * yy)
        if r >= h[4] * u * (1 - 1e-12) and r <= h[5] * u * (1 + 1e-12) then
            local R = (h[4] + h[5]) * u / 2
            local theta = math.atan(yy, xx)
            local v, br, bt = 0, 0, 0
            for _, row in ipairs(all[index]) do
                local order, rc, rs, tc, ts, mean = table.unpack(row)
                if order == 0 then
                    v = v + mean - R * tc * math.log(r / R)
                    bt = bt + R / r * tc
                else
                    local z = order * math.log(r / R)
                    check(
                        math.abs(z) < 600,
                        "radial harmonic continuation exceeds finite range; use a thinner gap"
                    )
                    local ch, sh =
                        (math.exp(z) + math.exp(-z)) / 2, (math.exp(z) - math.exp(-z)) / 2
                    local ac, as = -R / order * (rs * ch + tc * sh), R / order * (rc * ch - ts * sh)
                    local dc, ds = -R / r * (rs * sh + tc * ch), R / r * (rc * sh - ts * ch)
                    local c, s = math.cos(order * theta), math.sin(order * theta)
                    v = v + ac * c + as * s
                    br = br + order / r * (-ac * s + as * c)
                    bt = bt - dc * c - ds * s
                end
            end
            local c, s = xx / r, yy / r
            local bx, by = br * c - bt * s, br * s + bt * c
            return bx, by, v, bx / (4 * pi * 1e-7), by / (4 * pi * 1e-7)
        end
    end
    error("luafemm AGE: sample outside gap")
end
--- Replace the display strip by circular AGE boundaries in a FEM document.
-- The native file has no interpolation-edition tag; the caller must retain
-- that choice externally and open the answer with the matching FEMM release.
function N.document(m, d, t, labels)
    local o = m.sliding.options
    local bid = #d.boundaries + 1
    d.boundaries[bid] = {
        bdryname = m.native_age.ages[1].name,
        bdrytype = 6,
        innerangle = 0,
        outerangle = 0,
        a_0 = 0,
        a_1 = 0,
        a_2 = 0,
        phi = 0,
        c0 = 0,
        c1 = 0,
    }
    local function ring(p)
        local r = math.sqrt((p[1] - o.center_x) ^ 2 + (p[2] - o.center_y) ^ 2)
        if math.abs(r - o.inner_radius) < o.outer_radius * 1e-9 then
            return 1
        end
        if math.abs(r - o.outer_radius) < o.outer_radius * 1e-9 then
            return 2
        end
    end
    local segments = {}
    for _, s in ipairs(d.segments) do
        local p, q = d.points[s[1] + 1], d.points[s[2] + 1]
        if ring(p) and ring(p) == ring(q) then
            local a = math.atan(p[2] - o.center_y, p[1] - o.center_x)
            local b = math.atan(q[2] - o.center_y, q[1] - o.center_x)
            local delta = (b - a) % (2 * pi)
            local first, last = s[1], s[2]
            if delta > pi then
                first, last = last, first
                delta = 2 * pi - delta
            end
            d.arcs[#d.arcs + 1] =
                { first, last, math.deg(delta), 360 / #m.sliding.inner, bid, 0, 0 }
        else
            segments[#segments + 1] = s
        end
    end
    d.segments = segments
    local removed = {}
    for _, f in ipairs(t.faces) do
        local x, y = f.seed[1] - o.center_x, f.seed[2] - o.center_y
        local r = math.sqrt(x * x + y * y)
        if not f.hole and r > o.inner_radius and r < o.outer_radius then
            removed[labels[f.id] + 1] = true
            d.holes[#d.holes + 1] = { f.seed[1], f.seed[2], 0 }
        end
    end
    local map, new = {}, {}
    for i, l in ipairs(d.labels) do
        if not removed[i] then
            new[#new + 1] = l
            map[i] = #new
        end
    end
    for id, label in pairs(labels) do
        labels[id] = map[label + 1] and map[label + 1] - 1 or nil
    end
    d.labels = new
    d.header.comment = "luafemm native "
        .. m.native_age.edition
        .. " AGE; use the matching interpolation basis"
end
--- Add a sampled native gap to an immutable contour view. Bulk triangles retain
-- their native coordinates; separate vertex IDs prevent cross-operator joins.
function N.view(m)
    local v = { options = m.options, stats = m.stats, nodes = {}, A = {}, triangles = {} }
    for i, p in ipairs(m.nodes) do
        v.nodes[i] = p
        v.A[i] = m.A[i]
    end
    for _, e in ipairs(m.triangles) do
        if not e.analytical_gap then
            v.triangles[#v.triangles + 1] = e
        end
    end
    local function triangle(a, b, c)
        local p, q, r = v.nodes[a], v.nodes[b], v.nodes[c]
        local twice = (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])
        local aa, bb, cc = v.A[a], v.A[b], v.A[c]
        v.triangles[#v.triangles + 1] = {
            ids = { a, b, c },
            bx = (aa * (r[1] - q[1]) + bb * (p[1] - r[1]) + cc * (q[1] - p[1])) / twice,
            by = -(aa * (q[2] - r[2]) + bb * (r[2] - p[2]) + cc * (p[2] - q[2])) / twice,
        }
    end
    for _, age in ipairs(m.native_age.ages) do
        local h, u = age.header, m.native_age.unit
        local n = math.max(16, 2 * h[9]) * (m.options.contour_subdivisions or 1)
        local nr = 2 * (m.options.contour_subdivisions or 1)
        local rows = {}
        for i = 0, n do
            rows[i] = {}
            for j = 0, nr do
                local t = math.rad(h[6]) * i / n
                local r = (h[4] + (h[5] - h[4]) * j / nr) * u
                local x, y = h[7] * u + r * math.cos(t), h[8] * u + r * math.sin(t)
                local id = #v.nodes + 1
                v.nodes[id] = { x, y }
                local _, _, a = N.field(m, m.A, x, y)
                v.A[id] = a
                rows[i][j] = id
            end
        end
        for i = 0, n - 1 do
            for j = 0, nr - 1 do
                local a, b, c, d = rows[i][j], rows[i + 1][j], rows[i + 1][j + 1], rows[i][j + 1]
                triangle(a, d, c)
                triangle(a, c, b)
            end
        end
    end
    return v
end
return N
