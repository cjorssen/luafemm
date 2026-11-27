-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- LuaTeX boundary: validated key values, path capture and generated TeX.
-- Numerical modules do not access TeX globals. This module deliberately does.
-- @module luafemm-tex
local T = {}
local M = require("luafemm")
local abs, max = math.abs, math.max
local function check(value, message)
    assert(value, "luafemm: " .. message)
end

--- Resolve native documents through the active format, then the output folder.
function T.native_read(file)
    local ext = file:lower():match("%.([^.]+)$")
    check(ext == "femx" or ext == "ansx", "native file needs a .femx or .ansx extension")
    local resolved = resolvers and resolvers.findfile and resolvers.findfile(file)
        or kpse and kpse.find_file(file, "tex")
    if not resolved or resolved == "" then
        local candidate = T.output_file(file)
        local stream = io.open(candidate, "rb")
        if stream then
            stream:close()
            resolved = candidate
        else
            resolved = file
        end
    end
    return require("luafemm-extended").read(resolved, ext == "ansx" and "answer" or "problem")
end
function T.native_save(id, o)
    local m
    if o.physics == "electrostatics" then
        m = require("luafemm-electric-tex").native_model(id, o.model)
    else
        check(o.physics == "magnetism", "unknown persistence physics")
        local biot = require("luafemm-biot-tex")
        if o.biot_id and o.biot_id ~= "" or o.model ~= "" and biot.has_model(o.model) then
            check(
                not (o.model ~= "" and T.models and T.models[o.model] and biot.has_model(o.model)),
                "ambiguous magnetic native model name"
            )
            m = require("luafemm-biot-tex").native_model(o.biot_id, o.model)
        else
            m = o.model ~= "" and T.models and T.models[o.model] or M.current
        end
        check(m, "no magnetic model to save")
    end
    check(o.file:match("%S"), "native save requires a filename")
    local ext = o.kind == "answer" and ".ansx" or ".femx"
    check(o.file:sub(-#ext):lower() == ext, "native filename extension disagrees with kind")
    require("luafemm-extended").write(m, T.output_file(o.file), o.kind)
end

--- Parse a finite numeric key without evaluating user text as Lua code.
-- @tparam string value Numeric literal supplied through luaescapestring.
-- @tparam string key Name included in diagnostics.
-- @param default Value returned for an empty string, if supplied.
function T.number(value, key, default)
    if value == "" and default ~= nil then
        return default
    end
    local number = tonumber(value)
    check(
        number and number == number and abs(number) < math.huge,
        "expected a finite number for " .. key .. ", got " .. tostring(value)
    )
    return number
end

--- Apply one key-based exterior condition to the current model.
function T.boundary(side, options)
    for key, value in pairs(options) do
        if key ~= "type" then
            options[key] = T.number(value, "boundary " .. key)
        end
    end
    M.boundary(M.current, side, options)
end

--- Initialise the current picture from validated strings and the PGF transform.
function T.new_problem(
    options,
    transform,
    basis,
    cache,
    import_file,
    explicit,
    name,
    operating,
    sliding,
    import_format,
    import_bias,
    age_interpolation,
    time,
    native_file
)
    for key, value in pairs(options) do
        if value == "" then
            options[key] = nil
        elseif
            key ~= "mesher"
            and key ~= "interpolation"
            and key ~= "field_model"
            and key ~= "problem_type"
            and key ~= "material_edition"
            and key ~= "material_file"
            and key ~= "linear_solver"
            and key ~= "preconditioner"
            and key ~= "linear_ordering"
            and not (key == "linear_tolerance" and value == "auto")
        then
            options[key] = T.number(value, key)
        end
    end
    options.log = function(message)
        texio.write_nl("luafemm: " .. message)
    end
    -- Follow TeX's output directory for relative cache names. The standalone
    -- Lua API deliberately retains ordinary working-directory semantics.
    if
        cache
        and cache.file
        and status
        and status.output_directory
        and status.output_directory ~= ""
        and not cache.file:match("^[/\\]")
        and not cache.file:match("^%a:")
    then
        cache.file = status.output_directory .. "/" .. cache.file
    end
    options.cache = cache
    -- Only explicitly supplied numerical keys override a retained template.
    local numerical, retained = {}, {}
    for key in (explicit or ""):gmatch("[^,]+") do
        if require("luafemm-linear").option_defaults[key] ~= nil then
            numerical[key] = options[key]
        else
            retained[#retained + 1] = key
        end
    end
    local template_explicit = table.concat(retained, ",")
    local position = sliding and sliding.name ~= ""
    local temporal = time and (time.previous ~= "" or time.state ~= "")
    if time and not temporal then
        check(
            T.number(time.dt, "time step") == 0
                and time.motion == "fixed"
                and time.source_scale == ""
                and time.boundary_scale == ""
                and T.number(time.inertia, "inertia") == 0
                and T.number(time.load, "load torque") == 0
                and T.number(time.damping, "damping") == 0,
            "time stepping keys require previous state"
        )
    end
    if time and time.import_time ~= "" then
        check(import_file and import_file ~= "", "import time requires an ANS import")
    end
    local offset_x = sliding and T.number(sliding.offset_x or "0", "rotor offset x") or 0
    local offset_y = sliding and T.number(sliding.offset_y or "0", "rotor offset y") or 0
    check(
        position or temporal or (offset_x == 0 and offset_y == 0),
        "rotor offsets require a sliding band"
    )
    if not position and not temporal then
        check(
            not sliding or T.number(sliding.angle, "rotor angle") == 0,
            "rotor angle requires a sliding band"
        )
    end
    if not operating or operating.name == "" then
        check(
            position or temporal or not T.operating_currents or not next(T.operating_currents),
            "circuit currents require an operating point or sliding band"
        )
        check(
            position or not T.operating_phases or not next(T.operating_phases),
            "circuit phases require an operating point or sliding band"
        )
        check(
            not operating or operating.mode == "incremental",
            "permeability requires an operating point"
        )
    end
    if native_file and native_file ~= "" then
        check(
            not temporal
                and not position
                and (not operating or operating.name == "")
                and (not import_file or import_file == ""),
            "native file conflicts with another model selector"
        )
        check(
            (explicit or ""):gsub(",", "") == "",
            "native file owns physical and numerical settings"
        )
        M.current = T.native_read(native_file)
        check(M.current.options, "use femm/electrostatic for an electric native file")
        check(not M.current.biot, "use femm/biot for a prescribed-current native file")
    elseif temporal then
        check(
            not position
                and (not operating or operating.name == "")
                and (not import_file or import_file == ""),
            "time states cannot be combined with other problem selectors"
        )
        check(template_explicit == "", "a time state owns geometry and physical problem settings")
        check(not next(T.operating_phases or {}), "time states do not have excitation phases")
        local S = require("luafemm-time-state")
        local name = time.previous ~= "" and time.previous or time.state
        local state = T.time_states and T.time_states[name]
        check(state, "unknown time state " .. name)
        check(time.previous == "" or time.state == "", "choose state or previous state")
        if time.previous ~= "" then
            local mode = time.motion
            local o = {
                dt = T.number(time.dt, "time step"),
                source_scale = time.source_scale ~= ""
                        and T.number(time.source_scale, "source scale")
                    or nil,
                boundary_scale = time.boundary_scale ~= ""
                        and T.number(time.boundary_scale, "boundary scale")
                    or nil,
                motion = mode,
                currents = T.operating_currents,
                cache = cache,
                numerical = numerical,
            }
            if mode == "prescribed" then
                o.angle = T.number(sliding.angle, "rotor angle")
                o.offset_x = offset_x
                o.offset_y = offset_y
            else
                check(
                    T.number(sliding.angle, "rotor angle") == 0 and offset_x == 0 and offset_y == 0,
                    "position keys require motion=prescribed in a time step"
                )
            end
            if mode == "dynamic" then
                o.inertia = T.number(time.inertia, "inertia")
                o.load = T.number(time.load, "load torque")
                o.damping = T.number(time.damping, "damping")
            else
                check(
                    T.number(time.inertia, "inertia") == 0
                        and T.number(time.load, "load torque") == 0
                        and T.number(time.damping, "damping") == 0,
                    "mechanical coefficients require motion=dynamic"
                )
            end
            M.current = require("luafemm-time").new(state, o)
        else
            check(
                T.number(time.dt, "time step") == 0
                    and time.motion == "fixed"
                    and not next(T.operating_currents or {})
                    and not next(numerical)
                    and time.source_scale == ""
                    and time.boundary_scale == ""
                    and T.number(time.inertia, "inertia") == 0
                    and T.number(time.load, "load torque") == 0
                    and T.number(time.damping, "damping") == 0
                    and T.number(sliding.angle, "rotor angle") == 0
                    and offset_x == 0
                    and offset_y == 0,
                "viewing a state cannot advance or change its excitation"
            )
            M.current = S.restore(state)
        end
        M.current.options.log = options.log
    elseif position then
        check(
            not operating or operating.name == "",
            "cannot combine sliding band and operating point"
        )
        check(not import_file or import_file == "", "cannot combine sliding band and FEM import")
        check(template_explicit == "", "a sliding band owns its problem settings")
        local template = T.sliding_bands and T.sliding_bands[sliding.name]
        check(template, "unknown sliding band " .. sliding.name)
        M.current = require("luafemm-sliding").new(template, {
            numerical = numerical,
            angle = T.number(sliding.angle, "rotor angle"),
            offset_x = offset_x ~= 0 and offset_x or nil,
            offset_y = offset_y ~= 0 and offset_y or nil,
            phases = T.operating_phases,
            currents = T.operating_currents,
            cache = cache,
        })
    elseif operating and operating.name ~= "" then
        check(not import_file or import_file == "", "cannot combine operating point and FEM import")
        for key in template_explicit:gmatch("[^,]+") do
            check(
                key == "frequency",
                "an operating point owns its geometry and problem settings; only frequency and numerical policy may be overridden"
            )
        end
        local bias = T.operating_points and T.operating_points[operating.name]
        check(bias, "unknown operating point " .. operating.name)
        M.current = require("luafemm-operating").new(bias, {
            numerical = numerical,
            mode = operating.mode,
            currents = T.operating_currents,
            phases = T.operating_phases,
            frequency = options.frequency,
            cache = cache,
        })
    elseif import_file and import_file ~= "" then
        check(options.field_model ~= "ideal", "FEM import cannot be combined with ideal mode")
        local overrides = { cache = cache, log = options.log }
        for key in (explicit or ""):gmatch("[^,]+") do
            check(
                key ~= "xmin" and key ~= "xmax" and key ~= "ymin" and key ~= "ymax",
                "imported domain bounds come from the FEM file"
            )
            overrides[key] = options[key]
        end
        local file
        if resolvers and resolvers.findfile then
            file = resolvers.findfile(import_file)
        elseif kpse then
            file = kpse.find_file(import_file, "tex")
        end
        if file == "" then
            file = nil
        end
        if not file and status and status.output_directory and status.output_directory ~= "" then
            local candidate = status.output_directory .. "/" .. import_file
            local stream = io.open(candidate, "rb")
            if stream then
                stream:close()
                file = candidate
            end
        end
        local format = import_format or "auto"
        check(format == "auto" or format == "fem" or format == "ans", "unknown import format")
        if format == "auto" then
            format = import_file:lower():match("%.([^.]+)$")
            check(
                format == "fem" or format == "ans",
                "set import format for files without a .fem or .ans extension"
            )
        end
        check(
            not import_bias or import_bias == "" or format == "fem",
            "import bias file requires format=fem"
        )
        if age_interpolation and age_interpolation ~= "" then
            check(format == "ans", "age interpolation requires an ANS snapshot")
            overrides.age_interpolation = age_interpolation
        end
        if time and time.import_time ~= "" then
            check(format == "ans", "import time requires an ANS snapshot")
            overrides.time = T.number(time.import_time, "import time")
        end
        M.current = (format == "ans" and M.import_ans or M.import_fem)(
            file or import_file,
            overrides,
            import_bias
        )
    else
        if time then
            check(
                T.number(time.dt, "time step") == 0
                    and time.motion == "fixed"
                    and time.source_scale == ""
                    and time.boundary_scale == "",
                "time stepping keys require previous state"
            )
        end
        M.current = M.new(options)
    end
    if M.current.options.material_edition == "femm4" then
        T.register_material_choices("femm4")
    end
    T.models = T.models or {}
    T.model_count = (T.model_count or 0) + 1
    local id = tostring(T.model_count)
    T.models[id] = M.current
    if name and name ~= "" then
        check(not name:match("^%d+$"), "numeric model names are reserved")
        check(not T.models[name], "model name already registered: " .. name)
        T.models[name] = M.current
    end
    tex.sprint("\\expandafter\\def\\csname luafemm@picturemodel\\endcsname{" .. id .. "}")
    M.current.frame = require("luafemm-path").frame(transform, basis)
    -- Shape dimensions use the original coordinate basis, independently of
    -- text metrics and later local x/y vector changes.
    M.current.component_basis = require("luafemm-path").frame("{1}{0}{0}{1}{0pt}{0pt}", basis)
end

-- Time-state names have document lifetime. Reusing a name atomically replaces
-- its owned snapshot; already derived models and profiles remain independent.
function T.time_state(name, model, o)
    check(type(name) == "string" and name:match("%S"), "time state needs a name")
    local m = model ~= "" and T.models and T.models[model] or M.current
    check(m, "unknown state model")
    local options = { initial = o.initial }
    if o.time ~= "" then
        options.time = T.number(o.time, "initial time")
    end
    if o.speed ~= "" then
        options.speed = T.number(o.speed, "initial speed")
    end
    if o.band ~= "" then
        options.template = T.sliding_bands and T.sliding_bands[o.band]
        check(options.template, "unknown state motion template")
    end
    local S = require("luafemm-time-state")
    local state = S.capture(m, options)
    if o.file ~= "" then
        S.write(state, T.time_file(o.file))
    end
    T.time_states = T.time_states or {}
    T.time_states[name] = state
end
function T.time_file(file)
    if
        status
        and status.output_directory
        and status.output_directory ~= ""
        and not file:match("^[/\\]")
        and not file:match("^%a:")
    then
        return status.output_directory .. "/" .. file
    end
    return file
end
function T.restore_time_state(name, file)
    check(name:match("%S") and file:match("%S"), "restore state needs name and file")
    local state = require("luafemm-time-state").read(T.time_file(file))
    T.time_states = T.time_states or {}
    T.time_states[name] = state
end

--- Mesh and capture a named moving-band template after all physical paths.
-- A picture-level capture is bound to its originating model identifier.
function T.sliding_band(name, model, options)
    check(type(name) == "string" and name:match("%S"), "sliding band needs a name")
    T.sliding_bands = T.sliding_bands or {}
    check(not T.sliding_bands[name], "duplicate sliding band " .. name)
    local m = model and model ~= "" and T.models and T.models[model] or nil
    if not model or model == "" then
        m = M.current
    end
    check(m, "unknown model for sliding band")
    for k, v in pairs(options) do
        if k ~= "method" and k ~= "periodicity" then
            options[k] = T.number(v, k)
        end
    end
    T.sliding_bands[name] = require("luafemm-sliding").capture(m, options)
end

--- Capture a named, owned bias after completing the picture's solve.
-- The optional model id binds end-picture hooks to their originating picture.
function T.operating_point(name, model)
    check(type(name) == "string" and name:match("%S"), "operating point needs a name")
    T.operating_points = T.operating_points or {}
    check(not T.operating_points[name], "duplicate operating point " .. name)
    local m = M.current
    if model and model ~= "" then
        m = T.models and T.models[model]
    end
    check(m, "unknown model for operating point")
    M.solve(m)
    T.operating_points[name] = require("luafemm-operating").capture(m)
end

--- Print the current cache outcome; this command never starts a calculation.
function T.cache_status()
    check(M.current, "no current model")
    tex.sprint(M.current.cache_info.status)
end

--- Resolve one PGF material style; omitted overrides inherit library properties.
function T.material(name, library, mur, bh, remanence, coercivity, extra)
    local function optional(value, key)
        if value ~= "" then
            return T.number(value, key)
        end
    end
    local points = M.parse_pairs(bh)
    local properties = {
        library = library ~= "" and library or nil,
        mur = optional(mur, "mu r"),
        bh = next(points) and points or nil,
        remanence = optional(remanence, "remanence"),
        coercivity = optional(coercivity, "coercivity"),
    }
    for key, value in pairs(extra or {}) do
        if key == "lamination" or key == "wire_model" or key == "harmonic_law" then
            properties[key] = value ~= "inherit" and value or nil
        else
            properties[key] = optional(value, key)
        end
    end
    M.material(M.current, name, properties)
end

--- Print a solved-model statistic as plain numeric TeX input.
function T.statistic(key, mesh)
    local m = M.current
    check(m, "no current model")
    local source
    if mesh then
        source = m.mesh_stats
    else
        source = m.stats
    end
    local value = source and source[key]
    if not mesh and value == nil then
        value = require("luafemm-linear").statistic(m, key)
    end
    check(type(value) == "number", "statistic unavailable: " .. key)
    tex.sprint(string.format(mesh and "%.6g" or "%.8g", value))
end

--- Print the field magnitude in tesla at a point in model units.
function T.field_value(x, y)
    if M.current.options.problem_type == "axial" then
        local f = require("luafemm-harmonic").sample(M.current, T.number(x, "x"), T.number(y, "y"))
        tex.sprint(string.format("%.3f", math.abs(f.Bz[1])))
        return
    end
    local bx, by = M.sample(M.current, T.number(x, "x"), T.number(y, "y"))
    tex.sprint(string.format("%.3f", math.sqrt(bx * bx + by * by)))
end

--- Capture a closed material path and apply the PGF coordinate-precision convention.
-- @tparam table m Unmeshed model with a frame.
-- @tparam string s Evaluated soft-path tokens.
-- @tparam string material Model material name.
-- @tparam number current Current density in A/m^2.
-- @tparam number tolerance Polygonalisation tolerance in model units.
-- @tparam number angle Magnetisation direction in degrees.
-- @tparam number mesh_size Local spacing, zero for inheritance.
-- @tparam[opt] string turns Component winding ampere-turns (single contour only).
-- @tparam[opt="nonzero"] string fill_rule Evaluated native PGF fill rule.
function T.capture_path(
    m,
    s,
    material,
    current,
    tolerance,
    angle,
    mesh_size,
    turns,
    fill_rule,
    ideal_domain,
    static_options
)
    check(m and m.frame, "set femm/problem on the tikzpicture first")
    local contours = require("luafemm-path").read_contours(s, m.frame, tolerance)
    -- Round below the useful precision of PGF fixed-point coordinates.
    -- Snap domain endpoints to suppress TeX dimension-rounding artifacts.
    for _, points in ipairs(contours) do
        for _, p in ipairs(points) do
            for axis = 1, 2 do
                p[axis] = math.floor(p[axis] * 10000 + 0.5) / 10000
                local lo = axis == 1 and m.options.xmin or m.options.ymin
                local hi = axis == 1 and m.options.xmax or m.options.ymax
                if abs(p[axis] - lo) < 0.001 then
                    p[axis] = lo
                end
                if abs(p[axis] - hi) < 0.001 then
                    p[axis] = hi
                end
            end
        end
    end
    if turns and turns ~= "" then
        check(#contours == 1, "component winding needs one contour")
        current = require("luafemm-components").current(
            contours[1],
            m.options.unit,
            T.number(turns, "ampere turns")
        )
    end
    check(
        ideal_domain == nil or ideal_domain == "true" or ideal_domain == "false",
        "invalid ideal domain flag"
    )
    M.region_contours(m, material, contours, current, angle, mesh_size, fill_rule)
    require("luafemm-static").region(m.regions[#m.regions], static_options)
    m.regions[#m.regions].ideal_domain = ideal_domain == "true"
    -- The evaluated soft path is passive PGF token data, never executed on reload.
    -- Retain it before flattening so a later geometry compilation can reproduce
    -- the same primitives independently of the drawn approximation.
    m.regions[#m.regions].source = {
        soft_path = s,
        frame = { table.unpack(m.frame) },
        tolerance = tolerance,
        turns = turns,
    }
end

-- Saved shape IDs keep anchor geometry independent of the current model.
-- These records hold no solver/model reference and never enter disk caches.
local components = {}

--- Create immutable node geometry; PGF stores the returned numeric ID.
-- @tparam string kind Component type.
-- @tparam table options Escaped numeric strings from component keys.
-- @tparam[opt] table machine Structured machine options for rotor/stator nodes.
function T.component_new(kind, options, machine)
    check(M.current and M.current.frame, "set femm/problem before declaring components")
    check(not M.current.nodes, "declare components before meshing")
    if kind == "rotor" or kind == "stator" then
        check(M.current.options.field_model == "open", "machine nodes require field model=open")
        options = machine or options
    else
        for key, value in pairs(options) do
            options[key] = T.number(value, key)
        end
    end
    if kind ~= "rotor" and kind ~= "stator" then
        options.ideal = M.current.options.field_model == "ideal" and 1 or 0
    end
    check(
        M.current.options.problem_type == "planar",
        "component nodes are planar; draw meridional regions for axisymmetry"
    )
    local c = require("luafemm-components").make(kind, options)
    c.basis = M.current.component_basis
    components[#components + 1] = c
    tex.sprint(tostring(#components))
end

--- Install indexed shape anchors before PGF uses an anchor for node placement.
-- Definitions select the saved ID at lookup time; they do not capture a model.
-- @tparam integer id Saved component geometry ID.
function T.component_anchors(id)
    local c = assert(components[id])
    if c.machine then
        for anchor in pairs(c.anchors) do
            if anchor:match("^slot ") or anchor:match("^gap ") then
                tex.sprint(
                    "\\expandafter\\gdef\\csname pgf@anchor@femm component@"
                        .. anchor
                        .. "\\endcsname{\\csname luafemm@componentanchor\\endcsname{"
                        .. anchor
                        .. "}}"
                )
            end
        end
    end
end

local function component_point(c, p)
    local f = c.basis
    return f[1] * p[1] + f[3] * p[2], f[2] * p[1] + f[4] * p[2]
end

--- Emit one saved anchor in the shape's local PGF coordinate system.
function T.component_anchor(id, name)
    local c = assert(components[id])
    local p = c.anchors[name]
    check(p, "anchor '" .. name .. "' is unavailable on " .. c.kind)
    local x, y = component_point(c, p)
    tex.sprint(string.format("\\pgfpoint{%.9fpt}{%.9fpt}", x, y))
end

--- Intersect a connection ray with the component envelope, not its iron.
function T.component_border(id, x, y)
    local c = assert(components[id])
    local f = c.basis
    local function dimension(s)
        return assert(tonumber(s:match("([%d%.%+%-eE]+)pt")))
    end
    x, y = dimension(x), dimension(y)
    local det = f[1] * f[4] - f[2] * f[3]
    local u, v = (f[4] * x - f[3] * y) / det, (-f[2] * x + f[1] * y) / det
    local scale = max(abs(u) / c.half_width, abs(v) / c.half_height)
    if scale == 0 then
        tex.sprint("\\pgfpointorigin")
    else
        local px, py = component_point(c, { u / scale, v / scale })
        tex.sprint(string.format("\\pgfpoint{%.9fpt}{%.9fpt}", px, py))
    end
end

--- Draw/capture each part once, at the final node transform, before its text.
-- A component-local magnetization direction is mapped into the model frame.
-- Coil ampere-turns are converted after capture, using the rounded polygon.
function T.component_paths(id, transform, name)
    local c = assert(components[id])
    local n = require("luafemm-path").frame(transform, { "1pt", "0pt", "0pt", "1pt" })
    local f = M.current.frame
    local det = f[1] * f[4] - f[2] * f[3]
    local function mapped(p, vector)
        local x, y = component_point(c, p)
        x, y = n[1] * x + n[3] * y, n[2] * x + n[4] * y
        if not vector then
            x, y = x + n[5] - f[5], y + n[6] - f[6]
        end
        return { (f[4] * x - f[3] * y) / det, (-f[2] * x + f[1] * y) / det }
    end
    local function angle_value(degrees)
        local theta = math.rad(degrees)
        local v = mapped({ math.cos(theta), math.sin(theta) }, true)
        check(v[1] ^ 2 + v[2] ^ 2 > 0, "singular component transform")
        return math.deg(math.atan(v[2], v[1]))
    end
    if name and name ~= "" then
        M.current.component_cycles = M.current.component_cycles or {}
        local cycles = {}
        for cycle, points in pairs(c.cycles or {}) do
            cycles[cycle] = {}
            for _, p in ipairs(points) do
                cycles[cycle][#cycles[cycle] + 1] = mapped(p)
            end
        end
        M.current.component_cycles[name] = cycles
    end

    if c.machine then
        M.current.machine_nodes = M.current.machine_nodes or {}
        if name and name ~= "" then
            M.current.machine_nodes[name] = {
                center = mapped({ 0, 0 }),
                east = mapped({ c.airgap_radius, 0 }),
                north = mapped({ 0, c.airgap_radius }),
                kind = c.kind,
                tolerance = c.options.curve_tolerance,
                nominal_radius = c.airgap_radius,
            }
        end
    end
    for _, r in ipairs(c.regions) do
        local paths = {}
        for _, contour in ipairs(r.contours or { r.points }) do
            local points = {}
            for i, p in ipairs(contour) do
                local x, y = component_point(c, p)
                points[i] = string.format("(%.9fpt,%.9fpt)", x, y)
            end
            paths[#paths + 1] = table.concat(points, "--") .. "--cycle"
        end
        local path = table.concat(paths, " ")
        if c.machine then
            tex.sprint(
                string.format(
                    "\\csname luafemm@machinepart\\endcsname{%s}{%s}{%.17g}{%.17g}{%s}{%.17g}{%s}{%s}",
                    r.role,
                    r.turns and string.format("%.17g", r.turns) or "",
                    r.mesh_size,
                    angle_value(c.options.magnetization_angle),
                    r.material,
                    r.current_density or 0,
                    path,
                    r.style or ""
                )
            )
        else
            tex.sprint(
                string.format(
                    "\\csname luafemm@componentpart\\endcsname{%s}{%s}{%.17g}{%.17g}{%s}",
                    r.role,
                    r.turns and string.format("%.17g", r.turns) or "",
                    r.mesh_size,
                    angle_value(r.angle or c.options.magnetization_angle),
                    path
                )
            )
        end
    end
end

local function canvas(m, x, y)
    local f = m.frame
    return string.format("(%.7fpt,%.7fpt)", f[1] * x + f[3] * y + f[5], f[2] * x + f[4] * y + f[6])
end

local maps = {}
function T.visual(o, vectors)
    local V = require("luafemm-visual")
    local m = M.current
    local v = V.prepare(m, o)
    if vectors then
        for _, p in ipairs(V.vectors(v, o)) do
            tex.sprint(
                "\\draw[femm vector] "
                    .. canvas(m, p[1], p[2])
                    .. "--"
                    .. canvas(m, p[3], p[4])
                    .. ";"
            )
        end
    else
        check(o.colormap:match("^[%w%-]+$"), "colormap needs a simple PGFPlots identifier")
        local cells, lo, hi = V.cells(v, o)
        for i, q in ipairs(cells) do
            local points = {}
            for _, id in ipairs(m.triangles[i].ids) do
                local p = m.nodes[id]
                points[#points + 1] = canvas(m, p[1] / m.options.unit, p[2] / m.options.unit)
            end
            tex.sprint(
                string.format(
                    "\\pgfplotscolormapaccess[0:1]{%.12g}{%s}\\def\\luafemmcellcolor{\\csname pgfutil@definecolor\\endcsname{luafemmcell}{rgb}}\\expandafter\\luafemmcellcolor\\expandafter{\\pgfmathresult}\\path[fill=luafemmcell,femm scalar cell] %s--cycle;",
                    math.max(0, math.min(1, (q - lo) / (hi - lo))),
                    o.colormap,
                    table.concat(points, "--")
                )
            )
        end
        maps[o.name] = { minimum = lo, maximum = hi, colormap = o.colormap, unit = v.unit }
    end
end
function T.colorbar(name)
    local r = assert(maps[name], "luafemm visual: unknown colour-map name")
    token.set_macro("luafemm@mapmin", string.format("%.17g", r.minimum))
    token.set_macro("luafemm@mapmax", string.format("%.17g", r.maximum))
    token.set_macro("luafemm@mapcolor", r.colormap)
    token.set_macro("luafemm@mapunit", r.unit)
end

--- Print an exact section flux for PGF's generic number formatter.
function T.flux_number(name, snap_tolerance)
    local p = require("luafemm-profiles").get(name)
    if snap_tolerance == "auto" then
        snap_tolerance = p.model.options.problem_type == "axisymmetric" and 0 or 0.003
    elseif type(snap_tolerance) == "string" then
        snap_tolerance = T.number(snap_tolerance, "snap tolerance")
    end
    check(p.offset == 0, "flux requires offset=0; declare the displaced section explicitly")
    tex.sprint(string.format("%.12g", M.flux(p.model, p.points, snap_tolerance)))
end

--- Draw a component mean contour and optionally register its PGFPlots profile.
function T.mean_picture(m, options)
    local points = require("luafemm-ideal").mean(m, options.component, options.cycle, options.kind)
    if options.profile ~= "" then
        require("luafemm-profiles").register(
            m,
            options.profile,
            points,
            { samples = options.samples }
        )
    end
    local out = { "\\draw[femm mean path,femm current mean] " }
    for i, p in ipairs(points) do
        out[#out + 1] = (i > 1 and "--" or "") .. canvas(m, p[1], p[2])
    end
    out[#out + 1] = ";"
    tex.sprint(table.concat(out))
end

--- Solve if needed and emit contours in the original picture frame.
-- @tparam table m Model with frame.
-- @tparam number count Number of contour levels.
function T.tikz_field_picture(m, count, phase)
    if not m.stats then
        M.solve(m)
    end
    for _, path in ipairs(M.contours(m, count, nil, phase)) do
        local out = { "\\draw[femm field,femm current field] " }
        for i, p in ipairs(path) do
            out[#out + 1] = (i > 1 and "--" or "") .. canvas(m, p.x, p.y)
        end
        out[#out + 1] = ";"
        tex.sprint(table.concat(out))
    end
end

--- Render phase-specific physical eddy-current paths or normal glyphs.
-- @tparam table m Model with captured root frame.
-- @tparam table options Validated current-view keys.
function T.current_picture(m, options)
    local view = require("luafemm-currents").geometry(m, options)
    for _, path in ipairs(view.paths) do
        local out = {
            string.format(
                "\\draw[femm eddy line,femm current eddy,opacity=%.8g] ",
                path.opacity or 1
            ),
        }
        for i, p in ipairs(path) do
            out[#out + 1] = (i > 1 and "--" or "") .. canvas(m, p.x, p.y)
        end
        out[#out + 1] = ";"
        tex.sprint(table.concat(out))
    end
    for _, p in ipairs(view.symbols) do
        local centre = canvas(m, p.x, p.y)
        tex.sprint(
            string.format(
                "\\scope[shift={%s},femm eddy symbol,femm current eddy,opacity=%.8g]",
                centre,
                p.opacity
            )
        )
        tex.sprint("\\draw (0,0) circle[radius=\\pgfkeysvalueof{/femm/currents/symbol radius}];")
        if p.dot then
            tex.sprint(
                "\\fill[fill=.] (0,0) circle[radius=.28*\\pgfkeysvalueof{/femm/currents/symbol radius}];"
            )
        else
            tex.sprint(
                "\\draw (-.55*\\pgfkeysvalueof{/femm/currents/symbol radius},-.55*\\pgfkeysvalueof{/femm/currents/symbol radius})--(.55*\\pgfkeysvalueof{/femm/currents/symbol radius},.55*\\pgfkeysvalueof{/femm/currents/symbol radius}) (-.55*\\pgfkeysvalueof{/femm/currents/symbol radius},.55*\\pgfkeysvalueof{/femm/currents/symbol radius})--(.55*\\pgfkeysvalueof{/femm/currents/symbol radius},-.55*\\pgfkeysvalueof{/femm/currents/symbol radius});"
            )
        end
        tex.sprint("\\endscope")
    end
end

--- Repeat oriented field contours using a declared periodic rigid transform.
-- Antiperiodic copies reverse arrows; no material geometry or integral is copied.
function T.tikz_repeated_field(m, count, copies, pair_index)
    check(
        copies >= 1 and copies % 1 == 0 and copies <= 1000,
        "copies must be an integer in [1,1000]"
    )
    check(pair_index >= 1 and pair_index % 1 == 0, "invalid pair index")
    check(count >= 1 and count % 1 == 0, "line count must be a positive integer")
    M.solve(m)
    local pair = m.pairs and m.pairs[pair_index]
    check(pair, "unknown periodic pair")
    local levels
    if pair.sign == -1 then
        local peak = 0
        for _, a in ipairs(m.A) do
            peak = math.max(peak, math.abs(a))
        end
        if peak < 1e-20 then
            return
        end
        levels = {}
        for i = 1, count do
            levels[i] = -peak + 2 * peak * (i - 0.5) / count
        end
    end
    local paths = M.contours(m, count, levels)
    local c, s = math.cos(math.rad(pair.angle)), math.sin(math.rad(pair.angle))
    for k = 0, copies - 1 do
        for _, path in ipairs(paths) do
            local out = { "\\draw[femm field,femm current field] " }
            for j = 1, #path do
                local p = path[pair.sign == -1 and k % 2 == 1 and #path - j + 1 or j]
                local x, y = p.x, p.y
                for step = 1, k do
                    local a, b = x - pair.cx, y - pair.cy
                    x, y = pair.cx + c * a - s * b + pair.x, pair.cy + s * a + c * b + pair.y
                end
                out[#out + 1] = (j > 1 and "--" or "") .. canvas(m, x, y)
            end
            out[#out + 1] = ";"
            tex.sprint(table.concat(out))
        end
    end
end

--- Mesh if needed and emit each edge exactly once in the picture frame.
-- @tparam table m Model with frame.
function T.tikz_mesh_picture(m)
    if not m.nodes then
        M.mesh(m)
    end
    local seen = {}
    local out = { "\\draw[femm mesh,femm current mesh] " }
    for _, e in ipairs(m.triangles) do
        for i = 1, 3 do
            local a, b = e.ids[i], e.ids[i % 3 + 1]
            local k = math.min(a, b) .. ":" .. max(a, b)
            if not seen[k] then
                seen[k] = true
                local p, q = m.nodes[a], m.nodes[b]
                local u = m.options.unit
                out[#out + 1] = canvas(m, p[1] / u, p[2] / u)
                    .. "--"
                    .. canvas(m, q[1] / u, q[2] / u)
                    .. " "
            end
        end
    end
    out[#out + 1] = ";"
    tex.sprint(table.concat(out))
end

--- Export a named or current model; relative paths follow TeX's output directory.
function T.export(format, file, name, biasfile)
    local m = name ~= "" and T.models and T.models[name] or M.current
    check(m and (name == "" or T.models[name]), "unknown model " .. name)
    if
        status
        and status.output_directory
        and status.output_directory ~= ""
        and not file:match("^[/\\]")
        and not file:match("^%a:")
    then
        file = status.output_directory .. "/" .. file
    end
    if format == "fem" then
        M.export_fem(m, file, biasfile)
    elseif format == "ans" then
        M.export_ans(m, file, biasfile)
    else
        error("luafemm: unknown export format")
    end
end

--- Draw the compiled interfaces without registering physical geometry twice.
function T.geometry_picture(m)
    check(m and m.frame, "create a model before drawing its geometry")
    local t = m.topology or require("luafemm-topology").native(m)
    if m.sliding and m.sliding.sector then
        for _, s in ipairs(require("luafemm-sector-view").geometry(m, t)) do
            tex.sprint(
                "\\draw[femm geometry,femm current geometry,femm group/"
                    .. tostring(s.group or 0)
                    .. "/.try] "
                    .. canvas(m, s[1][1], s[1][2])
                    .. "--"
                    .. canvas(m, s[2][1], s[2][2])
                    .. ";"
            )
        end
        return
    end
    for _, s in ipairs(t.segments) do
        local p, q = t.points[s[1]], t.points[s[2]]
        tex.sprint(
            "\\draw[femm geometry,femm current geometry,femm group/"
                .. tostring(s.group or 0)
                .. "/.try] "
                .. canvas(m, p[1], p[2])
                .. "--"
                .. canvas(m, q[1], q[2])
                .. ";"
        )
    end
end

-- Identifiers contain only ASCII letters, digits and hyphens: safe PGF keys.
--- Install safe ASCII identifiers as PGF choice handlers at library load time.
-- The generic loader calls this with @ a letter; no document data are emitted.
function T.register_material_choices(edition)
    for _, r in ipairs(require("luafemm-materials").list(edition)) do
        tex.sprint(
            "\\pgfkeys{/femm/region/material/"
                .. r.id
                .. "/.code={\\csname luafemm@selectlibrary\\endcsname{"
                .. r.id
                .. "}}}"
        )
        tex.sprint(
            "\\pgfkeys{/tikz/femm/materials/library/"
                .. r.id
                .. "/.code={\\expandafter\\def\\csname luafemm@library\\endcsname{"
                .. r.id
                .. "}}}"
        )
    end
end

--- Register and draw an exact scan between two concentric circular components.
local function machine_gap(m, o)
    local machines = m.machine_nodes or {}
    local r, s = machines[o.rotor], machines[o.stator]
    check(
        r and s and r.kind == "rotor" and s.kind == "stator",
        "gap profile needs named rotor and stator"
    )
    local function radius(v)
        local x, y = v.east[1] - v.center[1], v.east[2] - v.center[2]
        local u, w = v.north[1] - v.center[1], v.north[2] - v.center[2]
        local rr = math.sqrt(x * x + y * y)
        check(
            math.abs(rr - math.sqrt(u * u + w * w)) < 1e-5 * rr
                and math.abs(x * u + y * w) < 1e-5 * rr * rr,
            "gap profile requires circular components (no shear or unequal scaling)"
        )
        return rr
    end
    local ri, ro = radius(r), radius(s)
    check(
        math.abs(r.center[1] - s.center[1]) + math.abs(r.center[2] - s.center[2]) < 0.003,
        "gap components must be concentric"
    )
    check(ri < ro, "rotor must fit strictly inside stator")
    local tol = math.max(r.tolerance * ri / r.nominal_radius, s.tolerance * ro / s.nominal_radius)
    check(ro - ri > 4 * tol, "reduce curve tolerance to resolve the air gap")
    return r.center, ri, ro, tol
end
--- Register a gap scan in the current model and emit its model-frame annotation.
-- @tparam table o Node names, position, profile name and angular sampling options.
function T.gap_profile(o)
    local m = M.current
    local center, ri, ro, tol = machine_gap(m, o)
    check(o.position > 0 and o.position < 1, "gap profile must lie strictly inside the gap")
    local rad = ri + (ro - ri) * o.position
    local P = require("luafemm-profiles")
    check(rad > ri + tol and rad < ro - tol, "profile too close to polygonal air-gap walls")
    local p = P.circle(m, o.name, center[1], center[2], rad, o)
    local points = {}
    for _, q in ipairs(p.points) do
        points[#points + 1] = canvas(m, q[1], q[2])
    end
    tex.sprint("\\draw[femm gap profile] " .. table.concat(points, "--") .. ";")
end

--- Add an air-only refinement band with margins from polygonal pole faces.
-- @tparam table o Rotor/stator names and a positive model-frame mesh_size.
function T.air_gap(o)
    local m = M.current
    check(not m.nodes, "declare air-gap refinement before solving")
    check(o.mesh_size > 0, "air-gap mesh size must be positive")
    local center, ri, ro, tol = machine_gap(m, o)
    local paths = {}
    for _, r in ipairs({ ri + 2 * tol, ro - 2 * tol }) do
        local n = math.max(36, math.ceil(math.pi / math.acos(1 - tol / r)))
        local points = {}
        for i = 0, n - 1 do
            local a = 2 * math.pi * i / n
            points[#points + 1] =
                canvas(m, center[1] + r * math.cos(a), center[2] + r * math.sin(a))
        end
        paths[#paths + 1] = table.concat(points, "--") .. "--cycle"
    end
    tex.sprint(
        string.format(
            "\\path[femm air gap,even odd rule,femm/region={material=air,current density=0,mesh size=%.17g}] %s;",
            o.mesh_size,
            table.concat(paths, " ")
        )
    )
end
--- Serialize a harmonic spectrum or scalar for generic TeX/PGFPlots.
-- @tparam table o Profile, component, order, normalization and scalar quantity.
-- @tparam boolean plot Emit coordinate pairs rather than a scalar.
function T.harmonic(o, plot)
    local h = require("luafemm-profiles").harmonics(
        o.profile,
        o.component,
        o.maximum_order,
        o.fundamental
    )
    if plot then
        local out = {}
        for i = 1, o.maximum_order do
            local value = h[i].amplitude
            if o.normalize == "true" then
                check(h[o.fundamental].amplitude > 1e-30, "cannot normalize a zero fundamental")
                value = value / h[o.fundamental].amplitude
            end
            out[#out + 1] = string.format("(%d,%.12g)", i, value)
        end
        tex.sprint(table.concat(out, " "))
    else
        local value
        if o.quantity == "thd" or o.quantity == "mean" then
            value = h[o.quantity]
        else
            check(h[o.order], "harmonic order exceeds maximum order")
            value = h[o.order][o.quantity]
        end
        check(type(value) == "number" and value == value, "undefined harmonic statistic")
        tex.sprint(string.format("%.12g", value))
    end
end
--- Capture a computational domain or a tagged path in the physical frame.
function T.static_path(kind, text, o)
    local m = M.current
    check(m and m.frame, "set femm/problem before paths")
    local path = require("luafemm-path")
    local F = require("luafemm-frontiers")
    local function snap(p)
        for _, v in ipairs(p) do
            for axis = 1, 2 do
                v[axis] = math.floor(v[axis] * 10000 + 0.5) / 10000
                local lo = axis == 1 and m.options.xmin or m.options.ymin
                local hi = axis == 1 and m.options.xmax or m.options.ymax
                if math.abs(v[axis] - lo) < 0.001 then
                    v[axis] = lo
                end
                if math.abs(v[axis] - hi) < 0.001 then
                    v[axis] = hi
                end
            end
        end
    end
    if kind == "domain" then
        local contours = path.read_contours(text, m.frame, o.tolerance)
        for _, p in ipairs(contours) do
            snap(p)
        end
        F.domain(m, contours, o.fill_rule)
    else
        local p, closed = path.read(text, m.frame, o.tolerance, true)
        if closed then
            table.remove(p)
        end
        snap(p)
        F.path(m, o.name, p, closed, o)
    end
end

--- Dispatch queued declarations after the problem begin hook.
function T.static_declare(kind, o)
    check(M.current, "set femm/problem first")
    if kind == "pair" then
        o.sign = o.type == "antiperiodic" and -1 or 1
        require("luafemm-frontiers").pair(M.current, o.first, o.second, o)
    else
        require("luafemm-static").circuit(M.current, o.name, o)
    end
end

--- Resolve relative generated files against the engine output directory.
function T.output_file(file)
    if
        file ~= ""
        and status
        and status.output_directory
        and status.output_directory ~= ""
        and not file:match("^[/\\]")
        and not file:match("^%a:")
    then
        return status.output_directory .. "/" .. file
    end
    return file
end
--- Copy a passive answer; resolve inputs through the format resolver, then the
-- engine output directory used by an earlier export in the same document.
function T.copy_answer(input, output, rings)
    local file
    if resolvers and resolvers.findfile then
        file = resolvers.findfile(input)
    elseif kpse then
        file = kpse.find_file(input, "tex")
    end
    if not file or file == "" then
        local candidate = T.output_file(input)
        local stream = io.open(candidate, "rb")
        if stream then
            stream:close()
            file = candidate
        else
            file = input
        end
    end
    return require("luafemm-ans").copy(file, T.output_file(output), rings)
end
return T
