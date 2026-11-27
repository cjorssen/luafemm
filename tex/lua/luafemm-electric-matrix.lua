-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Maxwell capacitance coefficients on one owned, unchanged finite-element mesh.
-- Q_i = offset_i + sum_j C_ij V_j. Unselected voltage electrodes retain their
-- operating voltages in the offset; charge-controlled conductors are eliminated
-- with their prescribed charges. Unit columns have zero independent sources.
-- @module luafemm-electric-matrix
local M = {}
local E = require("luafemm-electrostatic")
local F = require("luafemm-electric-fe")
local copy = require("luafemm-physical").copy
local function check(v, message)
    assert(v, "luafemm electric matrix: " .. message)
end
function M.calculate(model, o)
    o = o or {}
    check(model.fe, "finite-element model required")
    local m = E.snapshot(E.freeze(model))
    local names, selected = {}, {}
    local ports = o.electrodes
    if type(ports) == "string" then
        local list = {}
        for name in ports:gmatch("[^,]+") do
            list[#list + 1] = name:match("^%s*(.-)%s*$")
        end
        ports = #list > 0 and list or nil
    end
    for _, name in ipairs(ports or m.fe.electrodes) do
        local p = m.fe.electrodes[name]
        check(p, "unknown electrode: " .. tostring(name))
        if ports or p.condition == "voltage" then
            check(p.condition == "voltage", "matrix ports must be voltage-controlled")
            check(not selected[name], "duplicate matrix port")
            names[#names + 1], selected[name] = name, true
        end
    end
    check(#names > 0, "at least one voltage-controlled port is required")
    local normalization = o.normalization or "total"
    if normalization == "" then
        normalization = "total"
    end
    check(normalization == "total" or normalization == "per length", "unknown normalization")
    check(
        normalization ~= "per length" or not m.fe.axisymmetric,
        "per length is not an axisymmetric normalization"
    )
    local scale = normalization == "total" and m.fe.depth or 1
    local function solve(column)
        local trial = copy(m)
        local f = trial.fe
        trial.frozen, f.solution, f.imported_solution = false, nil, nil
        if trial.hybrid then
            trial.hybrid.solution = nil
        end
        for id, owner in pairs(f.owners) do
            if not f.floating[id] then
                if column then
                    f.fixed[id] = owner == column and 1 or 0
                elseif selected[owner] then
                    f.fixed[id] = 0
                end
            end
        end
        if column then
            if trial.hybrid then
                trial.sources = {}
            end
            f.point_loads = {}
            for _, t in ipairs(f.triangles) do
                t.density = 0
            end
            for _, e in ipairs(f.edges) do
                e.loads, e.displacement = 0, 0
            end
            for _, name in ipairs(f.electrodes) do
                if f.electrodes[name].condition == "charge" then
                    f.electrodes[name].charge = 0
                end
            end
        end
        F.solve(trial)
        return trial
    end
    local baseline = solve(nil)
    local result = {
        electrodes = names,
        coefficients = {},
        offset = {},
        normalization = normalization,
        units = normalization == "total" and "F" or "F/m",
        model = m,
    }
    for i, name in ipairs(names) do
        result.coefficients[i], result.offset[i] = {}, baseline.fe.charges[name] * scale
    end
    for j, name in ipairs(names) do
        local column = solve(name)
        for i, row in ipairs(names) do
            result.coefficients[i][j] = column.fe.charges[row] * scale
        end
    end
    return result
end
function M.value(matrix, row, column, quantity)
    local function index(name)
        for i, electrode in ipairs(matrix.electrodes) do
            if name == electrode then
                return i
            end
        end
        error("luafemm electric matrix: unknown port: " .. tostring(name))
    end
    local i = index(row)
    if quantity == "offset charge" then
        return matrix.offset[i]
    end
    check(quantity == nil or quantity == "capacitance", "unknown matrix quantity")
    return matrix.coefficients[i][index(column)]
end
return M
