-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Persistent mesh and solution reuse, independent of TeX and drawing state.
-- The private, versioned .lfc format is passive data, never a Lua program.
-- Exact canonical descriptors, rather than a short hash, decide cache hits.
-- An Adler-32 checksum detects accidental file damage, not hostile tampering.
-- @module luafemm-cache
local C = {}
-- Bump this revision for incompatible data changes OR numerical changes
-- within a development series whose public version string is unchanged.
local format_version = 6
local max_bytes = 256 * 1024 * 1024
local modes = { off = true, auto = true, refresh = true, frozen = true }
local function check(value, message)
    assert(value, "luafemm cache: " .. message)
end
local function finite(value)
    return type(value) == "number" and value == value and math.abs(value) < math.huge
end

-- Length-prefixed strings can contain arbitrary material names. Sorted keys
-- and 17 significant digits make descriptors deterministic and lossless for
-- the binary64 values used by LuaTeX, regardless of table insertion order.
local function encode(value)
    local output = {}
    local function emit(v, depth)
        check(depth <= 32, "data nesting limit exceeded")
        local kind = type(v)
        if kind == "string" then
            output[#output + 1] = "s" .. #v .. ":" .. v
        elseif kind == "number" then
            check(finite(v), "non-finite value")
            output[#output + 1] = "n" .. string.format("%.17g", v == 0 and 0 or v) .. ";"
        elseif kind == "boolean" then
            output[#output + 1] = v and "t" or "f"
        elseif kind == "table" then
            local keys = {}
            for key in pairs(v) do
                check(type(key) == "string" or finite(key), "unsupported table key")
                keys[#keys + 1] = key
            end
            table.sort(keys, function(a, b)
                if type(a) ~= type(b) then
                    return type(a) < type(b)
                end
                return a < b
            end)
            output[#output + 1] = "m" .. #keys .. ":"
            for _, key in ipairs(keys) do
                emit(key, depth + 1)
                emit(v[key], depth + 1)
            end
        else
            error("luafemm cache: unsupported value type " .. kind)
        end
    end
    emit(value, 0)
    return table.concat(output)
end

local function decode(text)
    local at, budget = 1, 5000000
    local function count()
        local stop = text:find(":", at, true)
        check(stop and stop - at < 12, "invalid length")
        local digits = text:sub(at, stop - 1)
        check(digits:match("^%d+$"), "invalid length")
        at = stop + 1
        return tonumber(digits)
    end
    local function read(depth)
        budget = budget - 1
        check(depth <= 32 and budget >= 0, "data complexity limit exceeded")
        local tag = text:sub(at, at)
        at = at + 1
        if tag == "s" then
            local length = count()
            check(length <= #text - at + 1, "truncated string")
            local value = text:sub(at, at + length - 1)
            at = at + length
            return value
        elseif tag == "n" then
            local stop = text:find(";", at, true)
            check(stop and stop - at <= 32, "invalid number")
            local value = tonumber(text:sub(at, stop - 1))
            check(finite(value), "invalid number")
            at = stop + 1
            return value
        elseif tag == "t" or tag == "f" then
            return tag == "t"
        elseif tag == "m" then
            local length = count()
            check(length <= budget / 2, "table size limit exceeded")
            local value = {}
            for _ = 1, length do
                local key = read(depth + 1)
                check(type(key) == "string" or finite(key), "invalid table key")
                check(value[key] == nil, "duplicate table key")
                value[key] = read(depth + 1)
            end
            return value
        end
        error("luafemm cache: invalid data tag")
    end
    local value = read(0)
    check(at == #text + 1, "trailing data")
    return value
end

local function checksum(text)
    local a, b = 1, 0
    for i = 1, #text do
        a = (a + text:byte(i)) % 65521
        b = (b + a) % 65521
    end
    return string.format("%08x", b * 65536 + a)
end

--- Validate copied constructor options; caching is opt-in.
-- Lua callbacks cannot be fingerprinted faithfully (closures may change).
-- @tparam table m Newly constructed model.
function C.configure(m)
    local o = m.options
    local c = o.cache or { mode = "off" }
    check(type(c) == "table", "cache must be a table")
    c.mode = c.mode or "auto"
    check(modes[c.mode], "unknown mode " .. tostring(c.mode))
    if c.mode ~= "off" then
        check(
            type(c.file) == "string" and c.file:match("%S") and not c.file:find("%c"),
            "file must be a nonempty filename stem without control characters"
        )
        check(
            not o.source and not o.boundary,
            "source/boundary callbacks require cache mode off; closures cannot be fingerprinted"
        )
    end
    o.cache = c
    m.cache_info =
        { status = c.mode == "off" and "off" or "pending", file = c.file and c.file .. ".lfc" }
end

local function descriptors(m)
    local o = m.options
    local mesh = { regions = {} }
    for _, key in ipairs({ "unit", "h", "xmin", "xmax", "ymin", "ymax", "mesher", "field_model" }) do
        mesh[key] = o[key]
    end
    -- Native ideal topology is derived from the region geometry below. Its
    -- face currents/materials are physical data, not reasons to remesh.
    mesh.topology = o.field_model ~= "ideal" and m.topology or nil
    if m.import_document then
        mesh.import_document = {}
        for key, value in pairs(m.import_document) do
            if key ~= "filename" then
                mesh.import_document[key] = value
            end
        end
    end
    mesh.curve_tolerance = o.curve_tolerance
    mesh.mesh_tolerance = o.mesh_tolerance or 0
    mesh.min_angle, mesh.max_area, mesh.max_steiner = o.min_angle, o.max_area, o.max_steiner
    local physical = {
        materials = {},
        regions = {},
        tolerance = o.tolerance,
        max_newton = o.max_newton,
        boundary = m.boundaries,
        ideal_sources = m.ideal_sources,
    }
    local function material(name)
        local p = m.materials[name]
        physical.materials[name] =
            { mur = p.mur, bh = p.bh, hc = p.hc, interpolation = p.interpolation }
    end
    material("air")
    for i, r in ipairs(m.regions) do
        mesh.regions[i] = {
            contours = r.contours,
            fill_rule = r.fill_rule,
            ideal_domain = r.ideal_domain,
            mesh_size = r.mesh_size,
            curve_tolerance = r.source and r.source.tolerance,
        }
        physical.regions[i] = { material = r.material, current = r.current, mx = r.mx, my = r.my }
        material(r.material)
    end
    return encode(mesh), encode(physical)
end

--- Exact passive signature used to reject stale solution exports.
function C.signature(m)
    local mesh, physical = descriptors(m)
    return mesh .. physical .. encode({ depth = m.options.depth })
end
C.encode = encode

local function status(m, value, reason)
    local info = m.cache_info
    if info.status ~= value or info.reason ~= reason then
        info.status, info.reason = value, reason
        if m.options.log then
            m.options.log(
                "cache "
                    .. value
                    .. ": "
                    .. (info.file or "")
                    .. (reason and " (" .. reason .. ")" or "")
            )
        end
    end
end

local function array(value, limit)
    check(type(value) == "table" and #value <= limit, "invalid array")
    local size = #value
    for key in pairs(value) do
        check(finite(key) and key % 1 == 0 and key >= 1 and key <= size, "invalid array index")
    end
    for i = 1, size do
        check(value[i] ~= nil, "sparse array")
    end
    return size
end

-- Validate before touching the live model. The checksum catches ordinary
-- corruption; bounds and topology checks also reject malformed data records.
local function validate(record, m, same_physics)
    local mesh = record.mesh
    check(type(mesh) == "table", "missing mesh")
    local n = array(mesh.nodes, 150000)
    local nt = array(mesh.elements, 300000)
    check(n >= 3 and nt >= 1 and type(mesh.boundary) == "table", "incomplete mesh")
    local o, area, edges, used = m.options, 0, {}, {}
    local x0, x1 = o.xmin * o.unit, o.xmax * o.unit
    local y0, y1 = o.ymin * o.unit, o.ymax * o.unit
    local epsilon = math.max(x1 - x0, y1 - y0) * 1e-9
    for id, p in ipairs(mesh.nodes) do
        check(array(p, 2) == 2 and finite(p[1]) and finite(p[2]), "invalid node")
        check(
            p[1] >= x0 - epsilon
                and p[1] <= x1 + epsilon
                and p[2] >= y0 - epsilon
                and p[2] <= y1 + epsilon,
            "node outside domain"
        )
        if mesh.boundary[id] and not m.topology then
            check(
                math.min(
                    math.abs(p[1] - x0),
                    math.abs(p[1] - x1),
                    math.abs(p[2] - y0),
                    math.abs(p[2] - y1)
                ) <= epsilon,
                "invalid exterior node"
            )
        end
    end
    for id, value in pairs(mesh.boundary) do
        check(
            finite(id) and id % 1 == 0 and id >= 1 and id <= n and value == true,
            "invalid boundary index"
        )
    end
    for _, ids in ipairs(mesh.elements) do
        check(array(ids, 3) == 3, "invalid element")
        for j, id in ipairs(ids) do
            check(finite(id) and id % 1 == 0 and id >= 1 and id <= n, "invalid node index")
            used[id] = true
            local other = ids[j % 3 + 1]
            check(finite(other), "invalid node index")
            local key = math.min(id, other) .. ":" .. math.max(id, other)
            local e = edges[key] or { count = 0, a = id, b = other }
            e.count = e.count + 1
            check(e.count <= 2, "non-manifold mesh")
            edges[key] = e
        end
        local p, q, r = mesh.nodes[ids[1]], mesh.nodes[ids[2]], mesh.nodes[ids[3]]
        local twice = (q[1] - p[1]) * (r[2] - p[2]) - (r[1] - p[1]) * (q[2] - p[2])
        check(finite(twice) and twice > 0, "nonpositive element area")
        area = area + twice / 2
    end
    local ne, exterior = 0, {}
    for _, e in pairs(edges) do
        ne = ne + 1
        if e.count == 1 then
            check(mesh.boundary[e.a] and mesh.boundary[e.b], "unprotected mesh boundary")
            exterior[e.a], exterior[e.b] = true, true
        end
    end
    for id = 1, n do
        check(used[id] and (not mesh.boundary[id] or exterior[id]), "invalid mesh vertex")
    end
    local expected = (x1 - x0) * (y1 - y0)
    if m.topology then
        expected = 0
        for _, f in ipairs(m.topology.faces) do
            if not f.hole then
                expected = expected + f.area * o.unit ^ 2
            end
        end
        for _, e in ipairs(mesh.elements) do
            local p, q, r = mesh.nodes[e[1]], mesh.nodes[e[2]], mesh.nodes[e[3]]
            local f = require("luafemm-topology").locate(
                m.topology,
                (p[1] + q[1] + r[1]) / (3 * o.unit),
                (p[2] + q[2] + r[2]) / (3 * o.unit)
            )
            check(f and not f.hole, "cached triangle outside domain")
        end
        check(
            n - ne + nt == require("luafemm-topology").characteristic(m.topology),
            "invalid domain topology"
        )
    else
        check(n - ne + nt == 1, "invalid mesh topology")
    end
    check(math.abs(area / expected - 1) < 1e-8, "invalid mesh coverage")
    if o.mesher == "delaunay" then
        check(
            type(mesh.stats) == "table" and type(mesh.stats.boundary_nodes) == "table",
            "missing mesh statistics"
        )
        check(
            encode(mesh.stats.boundary_nodes) == encode(mesh.boundary),
            "inconsistent boundary tags"
        )
        array(mesh.stats.segments, 300000)
        for _, segment in ipairs(mesh.stats.segments) do
            check(
                type(segment) == "table" and finite(segment[1]) and finite(segment[2]),
                "invalid constraint"
            )
            local key = math.min(segment[1], segment[2]) .. ":" .. math.max(segment[1], segment[2])
            check(edges[key], "missing constraint")
        end
        for _, key in ipairs({
            "seconds",
            "flips",
            "min_angle",
            "max_area",
            "mesh_tolerance",
            "exact_orientation",
            "exact_incircle",
            "exact_diametral",
            "steiner_points",
            "segment_splits",
            "circumcenters",
        }) do
            check(finite(mesh.stats[key]), "invalid mesh statistic " .. key)
        end
    end
    if record.solution then
        local solution = record.solution
        check(
            type(solution) == "table"
                and array(solution.A, n) == n
                and type(solution.stats) == "table",
            "invalid solution"
        )
        local fixed = same_physics
                and o.field_model ~= "ideal"
                and require("luafemm-boundary").prepare(m, mesh.nodes, mesh.elements, mesh.boundary)
            or {}
        if same_physics and o.field_model == "ideal" then
            require("luafemm-ideal").validate_potential(m, mesh.nodes, mesh.elements, solution.A)
        end
        for id, value in ipairs(solution.A) do
            check(finite(value) and (fixed[id] == nil or value == fixed[id]), "invalid potential")
        end
        for _, key in ipairs({ "newton", "cg", "residual", "seconds", "nodes", "triangles" }) do
            check(
                finite(solution.stats[key]) and solution.stats[key] >= 0,
                "invalid solution statistic"
            )
        end
        check(
            solution.stats.nodes == n and solution.stats.triangles == nt,
            "solution size mismatch"
        )
    end
end

local function read_file(file)
    local stream, message = io.open(file, "rb")
    check(stream, "cannot read " .. file .. ": " .. tostring(message))
    local text = stream:read(max_bytes + 1)
    stream:close()
    check(text and #text <= max_bytes, "empty or oversized file")
    local size, sum, start = text:match("^LUAFEMM%-CACHE (%d+) ([0-9a-f]+)\n()")
    check(start and tonumber(size) == #text - start + 1, "invalid header or truncated file")
    local payload = text:sub(start)
    check(checksum(payload) == sum, "checksum mismatch")
    return decode(payload)
end

--- Restore compatible data on demand; never import graphical state.
-- In auto mode a stale solution may still supply a reusable mesh. Frozen
-- mode requires the requested level and raises before mutating the model.
-- @tparam table m Current model with all regions declared.
-- @tparam boolean solution Whether a converged solution is requested.
-- @treturn boolean Whether the requested level was restored.
function C.restore(m, solution)
    local c = m.options.cache
    if c.mode == "off" then
        return false
    end
    if c.mode == "refresh" then
        status(m, "refresh")
        return false
    end
    local mesh_key, solution_key = descriptors(m)
    local record, reason = m._cache_record, m._cache_reason
    if record == nil then
        local ok, value = pcall(function()
            local r = read_file(m.cache_info.file)
            check(type(r) == "table" and r.format == format_version, "incompatible cache format")
            check(r.version == require("luafemm").version, "incompatible package version")
            check(r.mesh_key == mesh_key, "geometry or mesh settings changed")
            validate(r, m, r.solution_key == solution_key)
            return r
        end)
        record, reason = ok and value or false, not ok and tostring(value) or nil
        if reason then
            reason = reason:match("luafemm cache: (.*)") or reason
        end
        m._cache_record, m._cache_reason = record, reason
    end
    local match = record and record.mesh_key == mesh_key
    local solved = match and record.solution and record.solution_key == solution_key
    if c.mode == "frozen" then
        check(
            match and (not solution or solved),
            "frozen cache unavailable: "
                .. m.cache_info.file
                .. " ("
                .. (reason or "solution missing or physical/solver settings changed")
                .. ")"
        )
    end
    if not match then
        status(m, "miss", reason)
        return false
    end
    local femm = require("luafemm")
    if not m.nodes then
        local mesh = record.mesh
        -- prepare_mesh owns/mutates node tables; the saved record remains passive.
        local nodes, elements = {}, {}
        for i, p in ipairs(mesh.nodes) do
            nodes[i] = { p[1], p[2] }
        end
        for i, e in ipairs(mesh.elements) do
            elements[i] = { e[1], e[2], e[3] }
        end
        m.mesh_stats, m.xs, m.ys = mesh.stats, mesh.xs, mesh.ys
        femm.prepare_mesh(m, nodes, elements, mesh.boundary)
    end
    if solution and solved then
        m.A, m.stats = record.solution.A, record.solution.stats
        femm.update_fields(m)
        status(m, "solution")
        return true
    end
    status(m, "mesh", solution and "solution missing or physical/solver settings changed" or nil)
    return not solution
end

local serial = 0
local function write_file(file, record)
    local payload = encode(record)
    check(#payload < max_bytes - 100, "cache exceeds size limit")
    local text = "LUAFEMM-CACHE " .. #payload .. " " .. checksum(payload) .. "\n" .. payload
    -- Write beside the target, then rename on the same filesystem. A failed
    -- write or rename leaves the previous cache intact. Concurrent writers
    -- should use distinct stems; no global temporary directory is needed.
    serial = serial + 1
    local nonce = tostring({}):gsub("[^%w]", "")
    local temporary = file .. ".tmp-" .. os.time() .. "-" .. nonce .. "-" .. serial
    local stream, message = io.open(temporary, "wb")
    check(stream, "cannot write " .. temporary .. ": " .. tostring(message))
    local written, write_error = stream:write(text)
    local closed, close_error = stream:close()
    if not written or not closed then
        os.remove(temporary)
        error("luafemm cache: write failed: " .. tostring(write_error or close_error))
    end
    local renamed, rename_error = os.rename(temporary, file)
    if not renamed then
        os.remove(temporary)
        error("luafemm cache: cannot replace " .. file .. ": " .. tostring(rename_error))
    end
end

--- Atomically save the prepared mesh and, if converged, the solution.
-- Frames, plots, logs, closures and prepared material references are omitted.
-- @tparam table m Meshed or solved model.
function C.store(m)
    local c = m.options.cache
    if c.mode == "off" or c.mode == "frozen" then
        return
    end
    local mesh_key, solution_key = descriptors(m)
    local nodes, elements, boundary = {}, {}, {}
    for i, p in ipairs(m.nodes) do
        nodes[i] = { p[1], p[2] }
        if m.exterior_nodes[i] then
            boundary[i] = true
        end
    end
    for i, e in ipairs(m.triangles) do
        elements[i] = e.ids
    end
    local record = {
        format = format_version,
        version = require("luafemm").version,
        mesh_key = mesh_key,
        solution_key = solution_key,
        mesh = {
            nodes = nodes,
            elements = elements,
            boundary = boundary,
            stats = m.mesh_stats,
            xs = m.xs,
            ys = m.ys,
        },
        solution = m.stats and { A = m.A, stats = m.stats } or nil,
    }
    write_file(m.cache_info.file, record)
    m._cache_record, m._cache_reason = record, nil
end

return C
