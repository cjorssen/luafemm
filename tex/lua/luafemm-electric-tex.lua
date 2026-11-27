-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Typed TeX adapter for electric models, views and owned profiles.
-- IDs are carried by TeX-local macros, not by a mutable global current model.
-- Source data have picture lifetime; named snapshots/profiles have document
-- lifetime and are replaced atomically. Returned snapshots are always copies.
-- @module luafemm-electric-tex
local T = {}
local E = require("luafemm-electrostatic")
local R = require("luafemm-electric-view")
local P = require("luafemm-path")
local C = require("luafemm-physical")
local models, names, profiles, seeds, magnetic = {}, {}, {}, {}, {}
local matrices, forces = {}, {}
local function check(v, s)
    assert(v, "luafemm electric: " .. s)
end
local function model(id)
    return assert(models[tonumber(id)], "luafemm: no electric picture")
end
local function inverse(f, x, y)
    x, y = x - f[5], y - f[6]
    local det = f[1] * f[4] - f[2] * f[3]
    return (f[4] * x - f[3] * y) / det, (-f[2] * x + f[1] * y) / det
end
-- Decimal coordinate parsing and transformed PGF arithmetic use finite TeX
-- precision. The declared model-unit tolerance snaps only near the rectangular
-- FE boundary, never an interior material interface. Zero disables snapping.
local function snap(m, points)
    if not m.fe then
        return points
    end
    local tolerance = m.fe.boundary_tolerance
    for _, p in ipairs(points) do
        for j = 1, 2 do
            for _, bound in ipairs({ m.fe.bounds[j], m.fe.bounds[j + 2] }) do
                if math.abs(p[j] - bound) <= tolerance then
                    p[j] = bound
                end
            end
        end
    end
    return points
end
local function input_file(path)
    local resolved
    if resolvers and resolvers.findfile then
        resolved = resolvers.findfile(path)
    elseif kpse then
        resolved = kpse.find_file(path, "tex")
    end
    if resolved and resolved ~= "" then
        return resolved
    end
    local candidate = require("luafemm-tex").output_file(path)
    local stream = io.open(candidate, "rb")
    if stream then
        stream:close()
        return candidate
    end
    return path
end
function T.new(o, transform, basis)
    local m
    if o.problem_type == "axisymmetric" then
        if o.xmin_explicit ~= "1" then
            o.xmin = "0"
        end
        if o.depth_explicit ~= "1" then
            o.depth = nil
        end
    end
    if o.native_file and o.native_file ~= "" then
        check(
            o.snapshot == ""
                and o.file == ""
                and o.native_unit_explicit ~= "1"
                and o.classic_explicit ~= "1"
                and o.file_override ~= "1",
            "native file owns its physical model"
        )
        m = require("luafemm-tex").native_read(o.native_file)
        check(
            m.descriptor and m.descriptor.physics == "electrostatics",
            "use femm/problem for a magnetic native file"
        )
    elseif o.file and o.file ~= "" then
        check(o.snapshot == "", "file and snapshot are mutually exclusive")
        check(
            o.classic_explicit ~= "1" and o.file_override ~= "1",
            "imported file owns its physical model"
        )
        m = require("luafemm-electric-exchange").read(input_file(o.file), {
            length_unit = o.length_unit,
            edition = o.file_edition,
            mesh_size = o.import_mesh_size ~= "" and o.import_mesh_size or nil,
            curve_tolerance = o.import_curve_tolerance,
        })
    elseif o.snapshot ~= "" then
        check(names[o.snapshot], "unknown electric snapshot: " .. o.snapshot)
        check(o.classic_explicit ~= "1", "snapshot owns its geometry and sources")
        m = E.snapshot(names[o.snapshot])
        if m.fe and m.fe.axisymmetric and o.xmin_explicit ~= "1" then
            o.xmin = "0"
        end
        -- Restored physical units and permittivity cannot be redeclared.
        check(o.length_unit == m.descriptor.length_unit, "snapshot length unit mismatch")
        check(
            math.abs(E.number(o.epsilon_r, "relative permittivity") * E.epsilon0 / m.epsilon - 1)
                < 1e-12,
            "snapshot permittivity mismatch"
        )
        m.bounds = E.new({
            length_unit = o.length_unit,
            epsilon_r = o.epsilon_r,
            xmin = o.xmin,
            xmax = o.xmax,
            ymin = o.ymin,
            ymax = o.ymax,
        }).bounds
    else
        m = E.new(o)
    end
    m.frame = P.frame(transform, basis)
    models[#models + 1] = m
    magnetic[#models] = require("luafemm").current
    require("luafemm").current = nil
    tex.sprint("\\expandafter\\def\\csname luafemm@electricid\\endcsname{" .. #models .. "}")
end
function T.finish(id, name)
    local m = model(id)
    if name ~= "" then
        names[name] = E.snapshot(m)
    end
    require("luafemm").current = magnetic[tonumber(id)]
end
function T.snapshot(name)
    check(names[name], "unknown electric snapshot: " .. name)
    return E.snapshot(names[name])
end
--- Register the transformed centre, including noncentral node anchors.
-- PGF's final shape frame already includes positioning and transform shape.
local function source_center(id, transform, cx, cy)
    local m = model(id)
    local f = P.frame(transform, { "1pt", "0pt", "0pt", "1pt" })
    cx, cy = tonumber(cx:match("([%d%.%+%-eE]+)pt")), tonumber(cy:match("([%d%.%+%-eE]+)pt"))
    local x, y = inverse(m.frame, f[1] * cx + f[3] * cy + f[5], f[2] * cx + f[4] * cy + f[6])
    return m, x, y
end
function T.charge(id, transform, cx, cy, q, r)
    local m, x, y = source_center(id, transform, cx, cy)
    E.add(m, x, y, q, r)
end
--- Source-correction nodes use typed strengths and a separate physical guard.
-- Only the regular axis is snapped, using the existing declaration tolerance.
function T.source(id, transform, cx, cy, o)
    local m, x, y = source_center(id, transform, cx, cy)
    if m.fe and m.fe.axisymmetric and math.abs(x) <= m.fe.boundary_tolerance then
        x = 0
    end
    require("luafemm-electric-hybrid").add(m, x, y, o)
end
--- Capture closed dielectric regions or one connected physical boundary path.
function T.physical(id, kind, path, o)
    local m = model(id)
    check(m.fe, "physical electric paths require model=finite elements")
    local tolerance = E.positive(o.tolerance, "curve tolerance")
    local F = require("luafemm-electric-fe")
    if kind == "region" then
        local contours = P.read_contours(path, m.frame, tolerance)
        for _, c in ipairs(contours) do
            snap(m, c)
        end
        F.region(m, contours, o)
    else
        local points, closed = P.read(path, m.frame, tolerance, true)
        F.path(m, snap(m, points), closed, kind == "charge layer" and "charge" or kind, o)
    end
end
local function list(text, pairs)
    if not text or text == "" then
        return nil
    end
    local result = {}
    for entry in (text .. ","):gmatch("(.-),") do
        check(entry:match("%S"), "empty numeric list entry")
        if pairs then
            local x, y = entry:match("^%s*([^/]+)/([^/]+)%s*$")
            check(x, "seed syntax is x/y, x/y, ...")
            result[#result + 1] = { E.number(x, "seed x"), E.number(y, "seed y") }
        else
            result[#result + 1] = E.number(entry, "potential level")
        end
    end
    return result
end
local function viewoptions(o)
    for k, v in pairs(o) do
        if v == "" then
            o[k] = nil
        end
    end
    o.seeds = list(o.seeds, true)
    o.levels = list(o.levels, false)
    return o
end
local function draw(m, paths, style)
    local f = m.frame
    for _, path in ipairs(paths) do
        local tokens = { "\\draw[" .. (path.style or style) .. "] " }
        for i, p in ipairs(path.points) do
            tokens[#tokens + 1] = string.format(
                "%s(%.7fpt,%.7fpt)",
                i == 1 and "" or "--",
                f[1] * p[1] + f[3] * p[2] + f[5],
                f[2] * p[1] + f[4] * p[2] + f[6]
            )
        end
        tokens[#tokens + 1] = path.closed and "--cycle;" or ";"
        tex.sprint(table.concat(tokens))
    end
end
function T.draw(id, kind, o)
    local m = model(id)
    o = viewoptions(o)
    if kind == "mesh" then
        check(m.fe, "electric mesh requires model=finite elements")
        draw(m, require("luafemm-electric-fe").mesh_paths(m), "femm electric mesh")
    elseif kind == "geometry" and m.fe then
        draw(m, require("luafemm-electric-fe").geometry(m))
    elseif kind == "geometry" then
        check(m.classic, "electric geometry requires a symmetry model")
        draw(m, require("luafemm-electric-classics").geometry(m, o.plate_span))
    elseif kind == "lines" then
        o.seeds = o.seeds or seeds[tonumber(id)]
        local lines = R.lines(m, o)
        for _, line in ipairs(lines) do
            if line.start == "steps" or line.finish == "steps" then
                texio.write_nl(
                    "term and log",
                    "luafemm electric warning: streamline reached max steps"
                )
                break
            end
        end
        draw(m, lines, "femm electric field")
    elseif kind == "contours" then
        draw(m, R.contours(m, o), "femm equipotential")
    elseif kind == "vectors" then
        draw(m, R.vectors(m, o), "femm electric vector")
    else
        error("luafemm: unknown electric view")
    end
end
--- Capture owns evaluated geometry; later camera changes never call the provider.
function T.scene(id, options, view)
    local m = options.model and options.model ~= "" and T.snapshot(options.model) or model(id)
    view = viewoptions(view)
    -- Match ordinary electric drawing: explicit seed coordinates override a
    -- seed path captured in the current picture. Named snapshots own physics,
    -- not the originating picture's transient seed-path registry.
    if options.mode ~= "spatial" and (not options.model or options.model == "") then
        view.seeds = view.seeds or seeds[tonumber(id)]
    end
    local scene = require("luafemm-electric-scene").capture(m, options, view)
    require("luafemm-scene3d-tex").store(options.name, scene)
end
function T.capture(id, name, path, o, kind)
    local m = model(id)
    local tolerance = E.positive(o.tolerance, "curve tolerance")
    local points = snap(m, P.read(path, m.frame, tolerance, true))
    if kind == "seeds" then
        -- Arclength sampling is geometrical: it must not freeze the source list.
        local empty = E.new({ length_unit = m.descriptor.length_unit })
        local profile = R.profile(empty, points, { samples = o.samples })
        local result = {}
        for _, v in ipairs(profile.data) do
            result[#result + 1] = { v.x, v.y }
        end
        seeds[tonumber(id)] = result
    else
        check(name:match("%S"), "electric profile needs a name")
        profiles[name] = R.profile(m, points, o)
    end
end
function T.profile(name)
    check(profiles[name], "unknown electric profile: " .. name)
    return C.copy(profiles[name])
end
function T.coordinates(name, component, abscissa)
    local p = assert(profiles[name], "luafemm: unknown electric profile: " .. name)
    local components = { Er = "Ex", Ez = "Ey", Dr = "Dx", Dz = "Dy" }
    local axes = { r = "x", z = "y" }
    if components[component] or axes[abscissa] then
        check(
            p.model.fe and p.model.fe.axisymmetric,
            "radial/axial aliases require an axisymmetric FE profile"
        )
        component, abscissa = components[component] or component, axes[abscissa] or abscissa
    end
    check(
        ({ V = true, Ex = true, Ey = true, E = true, Dx = true, Dy = true, D = true })[component],
        "unknown electric component"
    )
    check(({ s = true, x = true, y = true, t = true })[abscissa], "unknown electric abscissa")
    local result = {}
    for _, f in ipairs(p.data) do
        result[#result + 1] = string.format(
            "(%.12g,%s)",
            f[abscissa],
            f[component] == f[component] and string.format("%.12g", f[component]) or "nan"
        )
    end
    return table.concat(result, " ")
end
function T.value(id, name, quantity, electrode)
    local m = name ~= "" and assert(names[name], "luafemm: unknown electric snapshot") or model(id)
    local value
    if
        quantity == "source interaction energy"
        or quantity == "source interaction energy per length"
    then
        check(m.hybrid, "source interaction energy requires model=source correction")
        check(
            not electrode or electrode == "",
            "source interaction energy is a whole-model quantity"
        )
        check(
            quantity ~= "source interaction energy per length" or not m.fe.axisymmetric,
            "per length is not an axisymmetric normalization"
        )
        value = require("luafemm-electric-interaction").calculate(m)
        if quantity == "source interaction energy" then
            value = value * m.fe.depth
        end
    else
        value = E.properties(m, electrode)[quantity]
    end
    check(value ~= nil, "quantity is not available with this model/normalization: " .. quantity)
    check(value == value and math.abs(value) < math.huge, "integral property overflow")
    tex.sprint(string.format("%.12g", value))
end
function T.matrix(id, o)
    check(o.name and o.name:match("%S"), "matrix needs a name")
    local m = o.model ~= "" and assert(names[o.model], "luafemm: unknown electric snapshot")
        or model(id)
    matrices[o.name] = require("luafemm-electric-matrix").calculate(m, o)
end
function T.matrix_snapshot(name)
    return C.copy(assert(matrices[name], "luafemm: unknown electric matrix"))
end
function T.matrix_value(name, row, column, quantity)
    local value = require("luafemm-electric-matrix").value(
        assert(matrices[name], "luafemm: unknown electric matrix"),
        row,
        column,
        quantity
    )
    check(value == value and math.abs(value) < math.huge, "matrix value overflow")
    tex.sprint(string.format("%.12g", value))
end
function T.matrix_coordinates(name, column)
    local matrix = assert(matrices[name], "luafemm: unknown electric matrix")
    local values = {}
    for i, row in ipairs(matrix.electrodes) do
        values[#values + 1] = string.format(
            "(%d,%.12g)",
            i,
            require("luafemm-electric-matrix").value(matrix, row, column)
        )
    end
    return table.concat(values, " ")
end
function T.force(id, name, path, o)
    check(name and name:match("%S"), "force contour needs a name")
    local m = model(id)
    local points, closed = P.read(path, m.frame, E.positive(o.tolerance, "curve tolerance"), true)
    check(closed, "force contour must be closed")
    -- The same physical-box tolerance used by declarations preserves a regular
    -- axis closure after PGF fixed-point rotations and inverse frame conversion.
    forces[name] = require("luafemm-electric-force").calculate(m, snap(m, points), o)
end
function T.force_snapshot(name)
    return C.copy(assert(forces[name], "luafemm: unknown electric force"))
end
function T.force_value(name, quantity)
    local value = assert(forces[name], "luafemm: unknown electric force")[quantity]
    check(value ~= nil, "force quantity is unavailable for this geometry: " .. quantity)
    check(value == value and math.abs(value) < math.huge, "force overflow")
    tex.sprint(string.format("%.12g", value))
end
--- Export an owned view of a named or picture-local model. Never replace the
-- imported source; native representation checks run before the atomic write.
function T.export(id, o)
    local m = o.model ~= "" and assert(names[o.model], "luafemm: unknown electric snapshot")
        or model(id)
    require("luafemm-electric-exchange").write(
        m,
        require("luafemm-tex").output_file(o.file),
        o.format,
        o.edition
    )
end
--- Resolve a model for native persistence without forcing a problem solve.
function T.native_model(id, name)
    return name ~= "" and assert(names[name], "luafemm: unknown electric snapshot") or model(id)
end
return T
