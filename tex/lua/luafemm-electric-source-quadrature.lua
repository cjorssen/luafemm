-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Adaptive smooth contrast loads, with positive Gauss--Duffy quadrature.
-- Geometry has already certified that no singular source touches the cell.
-- One rule and four child rules estimate local vector-integral error relative
-- to an absolute-integrand scale, avoiding division by a cancelling component.
-- Length factors give SI gradients times model-coordinate basis gradients;
-- the axisymmetric weight includes the complete revolution. Budget/depth
-- exhaustion fails explicitly instead of returning an unqualified partial load.
local Q = {}
local gauss = {
    { 0.06943184420297371, 0.17392742256872693 },
    { 0.33000947820757187, 0.32607257743127307 },
    { 0.66999052179242813, 0.32607257743127307 },
    { 0.93056815579702629, 0.17392742256872693 },
}
local function mid(a, b)
    return { (a[1] + b[1]) / 2, (a[2] + b[2]) / 2 }
end
function Q.integrator(m, base)
    local h, count = m.hybrid, 0
    local function rule(a, b, c, delta)
        count = count + 1
        assert(
            count <= h.quadrature_cells,
            "luafemm electric source correction: source quadrature cell budget exceeded"
        )
        local twice = (b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])
        local ix, iy, scale = 0, 0, 0
        for _, u in ipairs(gauss) do
            for _, v in ipairs(gauss) do
                local l1, l2, l3 = 1 - u[1], u[1] * (1 - v[1]), u[1] * v[1]
                local x, y = l1 * a[1] + l2 * b[1] + l3 * c[1], l1 * a[2] + l2 * b[2] + l3 * c[2]
                local s = assert(base(x, y))
                local measure = twice
                    * u[2]
                    * v[2]
                    * u[1]
                    * m.unit
                    * (m.fe.axisymmetric and 2 * math.pi * x * m.unit or 1)
                local ex = delta.xx * s.Ex + delta.xy * s.Ey
                local ey = delta.xy * s.Ex + delta.yy * s.Ey
                ix, iy = ix + measure * ex, iy + measure * ey
                scale = scale + math.abs(measure) * math.sqrt(ex * ex + ey * ey)
            end
        end
        return { ix, iy, scale }
    end
    local function integrate(a, b, c, delta, depth, coarse)
        coarse = coarse or rule(a, b, c, delta)
        local ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
        local cells = { { a, ab, ca }, { ab, b, bc }, { ca, bc, c }, { ab, bc, ca } }
        local values, total = {}, { 0, 0, 0 }
        for j, cell in ipairs(cells) do
            values[j] = rule(cell[1], cell[2], cell[3], delta)
            for k = 1, 3 do
                total[k] = total[k] + values[j][k]
            end
        end
        local error = math.sqrt((total[1] - coarse[1]) ^ 2 + (total[2] - coarse[2]) ^ 2)
        if error <= h.quadrature_tolerance * math.max(total[3], 1e-300) then
            return total[1], total[2]
        end
        assert(
            depth < h.quadrature_depth,
            "luafemm electric source correction: source quadrature tolerance not reached within maximum depth"
        )
        local ix, iy = 0, 0
        for j, cell in ipairs(cells) do
            local x, y = integrate(cell[1], cell[2], cell[3], delta, depth + 1, values[j])
            ix, iy = ix + x, iy + y
        end
        return ix, iy
    end
    return function(t)
        local n = m.fe.nodes
        return integrate(
            n[t.ids[1]],
            n[t.ids[2]],
            n[t.ids[3]],
            { xx = t.tensor.xx - h.epsilon_r, yy = t.tensor.yy - h.epsilon_r, xy = t.tensor.xy },
            0
        )
    end, function()
        return count
    end
end
return Q
