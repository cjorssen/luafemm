-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Passive magnetic FEMM 4.0 reader/writer and planar DC model adapter.
-- Never evaluates imported Lua or TeX. File indices are converted at the
-- adapter boundary; passive documents retain original units and metadata.
-- @module luafemm-fem
local F = {}
local topology = require("luafemm-topology")
F.units = {
    inches = 0.0254,
    millimeters = 0.001,
    centimeters = 0.01,
    meters = 1,
    mils = 0.0000254,
    microns = 0.000001,
}
local sections = {
    { "pointprops", "point", "points_properties" },
    { "bdryprops", "bdry", "boundaries" },
    { "blockprops", "block", "materials" },
    { "circuitprops", "circuit", "circuits" },
    { "numpoints", 4, "points" },
    { "numsegments", 6, "segments" },
    { "numarcsegments", 7, "arcs" },
    { "numholes", 3, "holes" },
    { "numblocklabels", 9, "labels" },
}
local function check(v, s)
    assert(v, "luafemm FEM: " .. s)
end
local function finite(n)
    return type(n) == "number" and n == n and math.abs(n) < math.huge
end
local function number(s)
    local n = tonumber(s)
    check(finite(n), "invalid number " .. tostring(s))
    return n
end
local function copy(t)
    if type(t) ~= "table" then
        return t
    end
    local r = {}
    for k, v in pairs(t) do
        r[k] = copy(v)
    end
    return r
end
F.copy = copy
local function value(s)
    s = s:match("^%s*(.-)%s*$")
    if s:sub(1, 1) == '"' then
        check(s:sub(-1) == '"', "unterminated quoted string")
        return s:sub(2, -2)
    end
    return tonumber(s) or s
end
local function row(s)
    local out = {}
    local pre, quoted = s:match('^(.-)%s*"(.*)"%s*$')
    for v in (pre or s):gmatch("%S+") do
        out[#out + 1] = number(v)
    end
    if quoted then
        out[#out + 1] = quoted
    end
    return out
end

--- Parse a passive document. Errors include the source name and physical line.
-- Unknown scalar metadata are retained; physical support is checked separately.
-- @tparam string text Complete file bytes (at most 64 MiB).
-- @tparam[opt] string name Source name for diagnostics.
-- @treturn table Passive document in file coordinates.
function F.parse(text, name)
    check(type(text) == "string" and #text <= 64 * 1024 * 1024, "input exceeds 64 MiB")
    text = text:gsub("^\239\187\191", "")
    local lines = {}
    for line in (text .. "\n"):gmatch("(.-)\n") do
        lines[#lines + 1] = line:gsub("\r$", "")
    end
    local at = 0
    local d = { header = {}, filename = name }
    local function nextline()
        repeat
            at = at + 1
        until at > #lines or lines[at]:match("%S")
        check(at <= #lines, "unexpected end of file")
        return lines[at]
    end
    local function count(v)
        local n = number(v)
        check(n >= 0 and n % 1 == 0 and n <= 300000, "invalid record count")
        return n
    end
    local function parse()
        local seen = {}
        while at < #lines do
            while at < #lines and not lines[at + 1]:match("%S") do
                at = at + 1
            end
            if at == #lines then
                break
            end
            local line = nextline()
            local key, v = line:match("^%s*%[([^%]]+)%]%s*=%s*(.-)%s*$")
            check(key, "expected a section header (ANS input is not supported)")
            key = key:lower()
            check(not seen[key], "duplicate section " .. key)
            seen[key] = true
            local section
            for _, s in ipairs(sections) do
                if s[1] == key then
                    section = s
                    break
                end
            end
            if not section then
                d.header[key] = value(v)
            else
                local list = {}
                d[section[3]] = list
                for _ = 1, count(v) do
                    if type(section[2]) == "number" then
                        local r = row(nextline())
                        check(#r >= section[2], "short " .. key .. " record")
                        list[#list + 1] = r
                    else
                        check(
                            nextline():lower():match("^%s*<begin" .. section[2] .. ">%s*$"),
                            "missing property start"
                        )
                        local p = {}
                        while true do
                            local l = nextline()
                            if l:lower():match("^%s*<end" .. section[2] .. ">%s*$") then
                                break
                            end
                            local k, w = l:match("^%s*<([^>]+)>%s*=%s*(.-)%s*$")
                            check(k, "invalid property")
                            k = k:lower()
                            check(
                                p[k] == nil and (k ~= "bhpoints" or p.bh == nil),
                                "duplicate property " .. k
                            )
                            if k == "bhpoints" then
                                p.bh = {}
                                for _ = 1, count(w) do
                                    local r = row(nextline())
                                    check(#r == 2, "invalid B-H pair")
                                    p.bh[#p.bh + 1] = r
                                end
                            else
                                p[k] = value(w)
                            end
                        end
                        list[#list + 1] = p
                    end
                end
            end
        end
        check(tonumber(d.header.format) == 4, "only magnetic format 4.0 is supported")
        for _, s in ipairs(sections) do
            check(d[s[3]], "missing " .. s[1])
        end
        check(F.units[tostring(d.header.lengthunits):lower()], "unknown length unit")
        for _, s in ipairs(d.segments) do
            for i = 1, 2 do
                check(s[i] % 1 == 0 and s[i] >= 0 and s[i] < #d.points, "invalid segment endpoint")
            end
            check(s[1] ~= s[2], "zero-length segment")
        end
        for _, s in ipairs(d.arcs) do
            for i = 1, 2 do
                check(s[i] % 1 == 0 and s[i] >= 0 and s[i] < #d.points, "invalid arc endpoint")
            end
            check(s[1] ~= s[2] and s[3] > 0 and s[3] <= 180 and s[4] > 0, "invalid circular arc")
        end
        return d
    end
    local ok, result = pcall(parse)
    if not ok then
        error((name or "<string>") .. ":" .. at .. ": " .. tostring(result), 0)
    end
    return result
end

--- Read a file without executing any of its contents.
-- @tparam string path Input filename in the current filesystem namespace.
-- @treturn table Passive document with its source filename.
function F.read(path)
    local f, e = io.open(path, "rb")
    check(f, "cannot open " .. path .. ": " .. tostring(e))
    local s = f:read(64 * 1024 * 1024 + 1)
    f:close()
    return F.parse(s, path)
end
local function scalar(v)
    if type(v) == "number" then
        check(finite(v), "non-finite output")
        return string.format("%.17g", v)
    end
    check(type(v) == "string" and not v:find('[\r\n"]'), "invalid quoted text")
    return '"' .. v .. '"'
end
local function keys(t)
    local r = {}
    for k in pairs(t) do
        r[#r + 1] = k
    end
    table.sort(r)
    return r
end

--- Serialize a passive document, retaining unknown scalar metadata.
-- @tparam table d Passive document from read/parse or document.
-- @treturn string Canonical LF-separated magnetic problem.
function F.serialize(d)
    local out = { "[Format] = 4.0" }
    for _, k in ipairs(keys(d.header)) do
        if k ~= "format" then
            local v = d.header[k]
            if k == "lengthunits" or k == "problemtype" or k == "coordinates" then
                out[#out + 1] = "[" .. k .. "] = " .. v
            else
                out[#out + 1] = "[" .. k .. "] = " .. scalar(v)
            end
        end
    end
    for _, s in ipairs(sections) do
        local list = d[s[3]] or {}
        out[#out + 1] = "[" .. s[1] .. "] = " .. #list
        for _, p in ipairs(list) do
            if type(s[2]) == "number" then
                local r = {}
                for _, v in ipairs(p) do
                    r[#r + 1] = scalar(v)
                end
                out[#out + 1] = table.concat(r, "\t")
            else
                out[#out + 1] = "<Begin" .. s[2] .. ">"
                for _, k in ipairs(keys(p)) do
                    if k == "bh" then
                        out[#out + 1] = "<BHPoints> = " .. #p.bh
                        for _, v in ipairs(p.bh) do
                            out[#out + 1] = scalar(v[1]) .. "\t" .. scalar(v[2])
                        end
                    else
                        out[#out + 1] = "<" .. k .. "> = " .. scalar(p[k])
                    end
                end
                out[#out + 1] = "<End" .. s[2] .. ">"
            end
        end
    end
    return table.concat(out, "\n") .. "\n"
end

--- Atomically replace an output file after all data have been validated.
-- @tparam string path Destination filename; its parent must exist.
-- @tparam string text Already validated complete output bytes.
-- @tparam[opt] string input Protected source filename.
function F.write_text(path, text, input)
    check(
        type(path) == "string" and path:match("%S") and not path:find("%c"),
        "invalid output filename"
    )
    check(not input or path ~= input, "refusing to overwrite imported source")
    if input then
        local ok, lfs = pcall(require, "lfs")
        if ok then
            local a, b = lfs.attributes(input), lfs.attributes(path)
            check(
                not (a and b and a.ino and a.ino ~= 0 and a.dev == b.dev and a.ino == b.ino),
                "refusing to overwrite imported source through an alias"
            )
        end
    end
    local tmp = path .. ".tmp-" .. tostring({}):gsub("[^%w]", "")
    local f, e = io.open(tmp, "wb")
    check(f, "cannot write " .. tmp .. ": " .. tostring(e))
    local ok, err = f:write(text)
    local closed, ce = f:close()
    if not ok or not closed then
        os.remove(tmp)
        error(tostring(err or ce))
    end
    local done, why = os.rename(tmp, path)
    if not done then
        os.remove(tmp)
        error(tostring(why))
    end
end

local function reference(list, id, what)
    check(
        finite(id) and id % 1 == 0 and id >= 0 and id <= #list,
        "invalid " .. what .. " reference"
    )
    return id ~= 0 and list[id] or nil
end
local function zero(p, fields, what)
    for _, k in ipairs(fields) do
        check((p[k] or 0) == 0, "unsupported " .. what .. " (" .. k .. ")")
    end
end

local known_header = {}
for k in
    (
        "format frequency precision minangle depth lengthunits problemtype coordinates acsolver "
        .. "prevsoln prevtype comment dosmartmesh forcemaxmesh extzo extro extri"
    ):gmatch("%S+")
do
    known_header[k] = true
end
local known_material = {}
for k in
    (
        "blockname mu_x mu_y h_c h_cangle j_re j_im sigma d_lam phi_h phi_hx phi_hy "
        .. "lamtype lamfill nstrands wired bh"
    ):gmatch("%S+")
do
    known_material[k] = true
end
local known_boundary = {}
for k in
    ("bdryname bdrytype a_0 a_1 a_2 phi c0 c1 c0i c1i mu_ssd sigma_ssd innerangle outerangle"):gmatch(
        "%S+"
    )
do
    known_boundary[k] = true
end
local function known(record, allowed, kind)
    for k in pairs(record) do
        check(allowed[k], "unsupported " .. kind .. " field " .. k)
    end
end

--- Convert a supported planar DC document into an owned unsolved Lua model.
-- Options override mesh controls; source coordinates are normalized directly.
-- @tparam table d Passive magnetic problem; copied into the model.
-- @tparam[opt] table options Explicit model-unit, mesh and solver overrides.
-- @treturn table Unmeshed model with tagged polygonal topology.
function F.model(d, options)
    check(
        not options or options.field_model ~= "ideal",
        "FEM import cannot be combined with ideal mode"
    )
    local M = require("luafemm")
    local o = copy(options or {})
    local h = d.header
    known(h, known_header, "problem")
    for _, p in ipairs(d.materials) do
        known(p, known_material, "material")
    end
    for _, p in ipairs(d.boundaries) do
        known(p, known_boundary, "boundary")
    end
    for _, p in ipairs(d.circuits) do
        known(
            p,
            { circuitname = true, totalamps_re = true, totalamps_im = true, circuittype = true },
            "circuit"
        )
    end
    local kind = tostring(h.problemtype):lower()
    check(kind == "planar" or kind == "axisymmetric", "unknown problem type")
    check(
        not o.problem_type or o.problem_type == kind,
        "imported problem type cannot be overridden"
    )
    o.problem_type = kind
    check((h.frequency or 0) == 0, "AC problems are unsupported")
    check(
        (h.prevsoln or "") == "" and (h.prevtype or 0) == 0,
        "previous-solution modes are unsupported"
    )
    check(#d.points_properties == 0, "point properties are unsupported")
    check(
        kind == "axisymmetric" or (finite(h.depth or 1) and (h.depth or 1) > 0),
        "invalid planar depth"
    )
    local unit = F.units[tostring(h.lengthunits):lower()]
    local modelunit = o.unit or 0.001
    local scale = unit / modelunit
    local raw = {}
    local xmin, ymin, xmax, ymax = math.huge, math.huge, -math.huge, -math.huge
    local function point(x, y)
        local p = { x * scale, y * scale }
        xmin = math.min(xmin, p[1])
        xmax = math.max(xmax, p[1])
        ymin = math.min(ymin, p[2])
        ymax = math.max(ymax, p[2])
        return p
    end
    local points = {}
    for i, p in ipairs(d.points) do
        check(p[3] == 0, "point conditions are unsupported")
        points[i] = point(p[1], p[2])
    end
    local periodic_counts = {}
    local function condition(id)
        local b = reference(d.boundaries, id, "boundary")
        if not b then
            return nil
        end
        zero(b, { "c0i", "c1i" }, "complex boundary")
        if b.bdrytype == 0 then
            check(
                tostring(h.coordinates or "cartesian"):lower() == "cartesian"
                    or ((b.a_1 or 0) == 0 and (b.a_2 or 0) == 0),
                "polar-varying potential is unsupported"
            )
            local phase = math.cos(math.rad(b.phi or 0))
            return {
                type = "dirichlet",
                value = (b.a_0 or 0) * phase,
                dx = (b.a_1 or 0) * phase / unit,
                dy = (b.a_2 or 0) * phase / unit,
                coefficient = 0,
            }
        end
        if b.bdrytype == 4 or b.bdrytype == 5 then
            periodic_counts[id] = (periodic_counts[id] or 0) + 1
            return {
                type = "periodic",
                path = "fem-pair-" .. id .. "-" .. periodic_counts[id],
                fem_group = id,
                sign = b.bdrytype == 4 and 1 or -1,
                value = 0,
                dx = 0,
                dy = 0,
                coefficient = 0,
            }
        end
        check(b.bdrytype == 2, "unsupported boundary type " .. tostring(b.bdrytype))
        check((b.c0 or 0) >= 0, "negative Robin coefficient is unsupported")
        return { type = "robin", value = b.c1 or 0, dx = 0, dy = 0, coefficient = b.c0 or 0 }
    end
    for _, s in ipairs(d.segments) do
        raw[#raw + 1] = {
            points[s[1] + 1],
            points[s[2] + 1],
            condition = condition(s[4]),
            hidden = s[5],
            group = s[6],
        }
    end
    local tolerance = o.curve_tolerance or 0.03
    check(finite(tolerance) and tolerance > 0, "invalid curve tolerance")
    for _, s in ipairs(d.arcs) do
        local boundary = reference(d.boundaries, s[5], "boundary")
        check(
            not boundary or (boundary.bdrytype ~= 4 and boundary.bdrytype ~= 5),
            "periodic arc import is unsupported; use paired straight segments"
        )
        local a, b = points[s[1] + 1], points[s[2] + 1]
        local theta = math.rad(s[3])
        local dx, dy = b[1] - a[1], b[2] - a[2]
        local chord = math.sqrt(dx * dx + dy * dy)
        local cx, cy =
            (a[1] + b[1]) / 2 - dy / (2 * math.tan(theta / 2)),
            (a[2] + b[2]) / 2 + dx / (2 * math.tan(theta / 2))
        local radius = chord / (2 * math.sin(theta / 2))
        local start = math.atan(a[2] - cy, a[1] - cx)
        local step = math.min(math.rad(s[4]), 2 * math.acos(math.max(-1, 1 - tolerance / radius)))
        check(step > 0, "arc tolerance is below floating-point resolution")
        local n = math.ceil(theta / step)
        check(n <= 40000, "arc subdivision limit exceeded")
        local prev = a
        for j = 1, n do
            local ang = start + theta * j / n
            local p = j == n and b
                or point(
                    (cx + radius * math.cos(ang)) / scale,
                    (cy + radius * math.sin(ang)) / scale
                )
            raw[#raw + 1] = { prev, p, condition = condition(s[5]), hidden = s[6], group = s[7] }
            prev = p
        end
    end
    check(xmax > xmin and ymax > ymin, "empty or degenerate domain")
    o.unit = modelunit
    o.xmin, o.xmax, o.ymin, o.ymax = xmin, xmax, ymin, ymax
    o.h = o.h or math.max(xmax - xmin, ymax - ymin) / 12
    o.mesher = "delaunay"
    if kind == "planar" then
        o.depth = o.depth or (h.depth or 1) * scale
    else
        o.depth = nil
    end
    o.tolerance = o.tolerance or h.precision
    local m = M.new(o)
    m.import_document = copy(d)
    for i, p in ipairs(d.materials) do
        zero(p, { "j_im", "phi_h", "phi_hx", "phi_hy" }, "material physics")
        check(
            (p.lamtype or 0) >= 0 and (p.lamtype or 0) <= 2,
            "unsupported imported lamination type"
        )
        check(
            (p.lamtype or 0) ~= 0 or (p.lamfill or 1) == 1,
            "nonlinear fill mixture import is unsupported"
        )
        local name = "fem:" .. i
        M.material(m, name, {
            mur = p.mu_x or 1,
            -- FEMM's directional sheet models use the permeability along
            -- the sheets as the base in both principal directions.
            mur_x = (p.lamtype == 2 and p.mu_y or p.mu_x) or 1,
            mur_y = (p.lamtype == 1 and p.mu_x or p.mu_y) or 1,
            conductivity = (p.sigma or 0) * 1e6,
            lamination = p.lamtype == 1 and "parallel x"
                or p.lamtype == 2 and "parallel y"
                or "none",
            lamination_fill = p.lamfill or 1,
            current = (p.j_re or 0) * 1e6,
            coercivity = p.h_c or 0,
            bh = p.bh and #p.bh > 0 and p.bh or nil,
            interpolation = o.interpolation or "femm",
            fem_name = p.blockname,
        })
    end
    local t = topology.arrange(
        raw,
        o.mesh_tolerance and o.mesh_tolerance > 0 and o.mesh_tolerance
            or math.max(xmax - xmin, ymax - ymin) * 1e-12
    )
    for i, c in ipairs(d.circuits) do
        check(
            c.circuittype == 1,
            "parallel FEMM circuit import is unsupported; export resolved series windings"
        )
        zero(c, { "totalamps_im" }, "complex circuit")
        require("luafemm-static").circuit(
            m,
            c.circuitname or ("fem-circuit-" .. i),
            { current = c.totalamps_re or 0 }
        )
    end
    local default
    for i, l in ipairs(d.labels) do
        reference(d.materials, l[3], "material")
        reference(d.circuits, l[5], "circuit")
        check(l[3] > 0, "missing material label")
        check((l[9] or 0) % 2 == 0, "external axisymmetric region is unsupported")
        require("luafemm-static-exchange").parse_expression(l[10], scale)
        local f = topology.locate(t, l[1] * scale, l[2] * scale, true)
        check(f, "label outside bounded domain")
        check(not f.label, "duplicate/conflicting block labels")
        f.label = i
        if math.floor((l[9] or 0) / 2) % 2 == 1 then
            check(not default, "multiple default labels")
            default = i
        end
    end
    for _, p in ipairs(d.holes) do
        local f = topology.locate(t, p[1] * scale, p[2] * scale, true)
        check(f and not f.label and not f.hole, "invalid or conflicting hole seed")
        f.hole = true
    end
    for _, f in ipairs(t.faces) do
        if not f.hole then
            local id = f.label or default
            check(id, "unlabeled region without default material")
            local l = d.labels[id]
            local current = m.materials["fem:" .. l[3]].current
            local circuit = reference(d.circuits, l[5], "circuit")
            if circuit then
                check(circuit.circuittype == 1, "parallel circuits are unsupported")
                zero(circuit, { "totalamps_im" }, "complex circuit")
                check(f.label, "default labels with circuits are unsupported")
                current = (circuit.totalamps_re or 0) * (l[8] or 1) / (f.area * modelunit ^ 2)
            end
            f.material, f.current, f.angle = "fem:" .. l[3], current, l[6]
            f.mx, f.my = math.cos(math.rad(f.angle)), math.sin(math.rad(f.angle))
            f.mesh_size = l[4] > 0 and l[4] * scale or 0
            f.group = l[7]
            M.region_contours(
                m,
                f.material,
                f.contours,
                f.current,
                f.angle,
                f.mesh_size,
                "even odd"
            )
            local region = m.regions[#m.regions]
            region.face, f.region = f.id, #m.regions
            local law, offset = require("luafemm-static-exchange").parse_expression(l[10], scale)
            if law then
                region.angle = offset
                region.mx, region.my = math.cos(math.rad(offset)), math.sin(math.rad(offset))
                f.angle, f.mx, f.my = region.angle, region.mx, region.my
                require("luafemm-static").region(region, law)
            end
            if circuit then
                region.current = 0
                require("luafemm-static").region(
                    region,
                    { circuit = circuit.circuitname or ("fem-circuit-" .. l[5]), turns = l[8] or 1 }
                )
                if law then
                    require("luafemm-static").region(region, law)
                end
            end
        end
    end
    m.topology = t
    require("luafemm-static").excite(m, t)
    local periodic_groups = {}
    for _, segment in ipairs(t.segments) do
        local c = segment.condition
        if c and c.fem_group then
            local group = periodic_groups[c.fem_group] or {}
            group[#group + 1] = segment
            periodic_groups[c.fem_group] = group
        end
    end
    check(
        kind ~= "axisymmetric" or not next(periodic_groups),
        "axisymmetric periodic input is unsupported"
    )
    if next(periodic_groups) then
        require("luafemm-static-exchange").import_pairs(m, t, periodic_groups)
    end
    return m
end

--- Compile the exact computational polygonal problem for FEM/ANS output.
-- Region currents become effective material records. No mesh is generated.
-- @tparam table m Native or imported model with positive planar depth.
-- @return Passive FEM document, computational topology, face-to-label map.
function F.document(m)
    if m.sliding then
        require("luafemm-sliding").validate(m)
    end
    check(
        m.options.problem_type ~= "axial",
        "axial H diffusion has no FEM vector-potential representation"
    )
    check(not m.operating, "operating-point FEM export is not implemented")
    check(m.options.frequency == 0, "harmonic FEM export is not implemented")
    check(m.options.field_model ~= "ideal", "ideal floating boundaries cannot be exported to FEMM")
    check(not m.options.source and not m.options.boundary, "callbacks cannot be exported exactly")
    for _, p in pairs(m.materials) do
        check(
            not require("luafemm-losses").active(p),
            "loss material metadata cannot be exported by the static FEM codec"
        )
    end
    local o = m.options
    local axi = o.problem_type == "axisymmetric"
    if axi then
        require("luafemm-axisymmetric").validate(m)
    end
    check(axi or (finite(o.depth) and o.depth > 0), "set a positive problem depth before export")
    local unitname
    for name, u in pairs(F.units) do
        if math.abs(u / o.unit - 1) < 1e-12 then
            unitname = name
        end
    end
    check(unitname, "model length unit cannot be represented in FEMM")
    require("luafemm-frontiers").compile(m)
    local t = require("luafemm-static-exchange").pairs(m, m.topology or topology.native(m))
    local d = {
        header = {
            format = 4,
            frequency = 0,
            precision = o.tolerance,
            minangle = math.max(1, o.min_angle),
            depth = o.depth or 1,
            lengthunits = unitname,
            problemtype = o.problem_type,
            coordinates = "cartesian",
            acsolver = 0,
            prevtype = 0,
            prevsoln = "",
            comment = "luafemm static computational geometry; voltage/parallel windings resolved at the operating point",
        },
        points_properties = {},
        boundaries = {},
        materials = {},
        circuits = {},
        points = {},
        segments = {},
        arcs = {},
        holes = {},
        labels = {},
    }
    for i, p in ipairs(t.points) do
        d.points[i] = { p[1], p[2], 0, 0 }
    end
    local periodic_ids = {}
    for _, s in ipairs(t.segments) do
        local c = s.condition
        local id = 0
        if c then
            check(
                c.type == "dirichlet" or (c.dx == 0 and c.dy == 0),
                "varying mixed data cannot be exported exactly"
            )
            id = c.fem_pair and periodic_ids[c.fem_pair] or #d.boundaries + 1
            if c.fem_pair then
                periodic_ids[c.fem_pair] = id
            end
            local b = {
                bdryname = "boundary " .. id,
                bdrytype = c.type == "periodic" and (c.sign == 1 and 4 or 5)
                    or c.type == "dirichlet" and 0
                    or 2,
                a_0 = 0,
                a_1 = 0,
                a_2 = 0,
                phi = 0,
                c0 = 0,
                c1 = 0,
                c0i = 0,
                c1i = 0,
                mu_ssd = 0,
                sigma_ssd = 0,
            }
            if c.type == "dirichlet" then
                b.a_0, b.a_1, b.a_2 = c.value, c.dx * o.unit, c.dy * o.unit
            else
                b.c0, b.c1 = c.coefficient, c.value
            end
            d.boundaries[id] = b
        end
        d.segments[#d.segments + 1] = { s[1] - 1, s[2] - 1, -1, id, s.hidden or 0, s.group or 0 }
    end
    local circuit_ids = {}
    local names = {}
    for name in pairs(m.static_circuits or {}) do
        names[#names + 1] = name
    end
    table.sort(names)
    for _, name in ipairs(names) do
        local c = m.static_circuits[name]
        local windings = {}
        for w in pairs(c.windings) do
            windings[#windings + 1] = w
        end
        table.sort(windings)
        circuit_ids[name] = {}
        for _, w in ipairs(windings) do
            local id = #d.circuits + 1
            circuit_ids[name][w] = id
            local label = #windings == 1 and name or name .. ":" .. w
            d.circuits[id] = {
                circuitname = label,
                totalamps_re = c.windings[w].current,
                totalamps_im = 0,
                circuittype = 1,
            }
        end
    end
    local labels = {}
    for _, f in ipairs(t.faces) do
        if f.hole then
            d.holes[#d.holes + 1] = { f.seed[1], f.seed[2], 0 }
        else
            local p = m.materials[f.material]
            local id = #d.materials + 1
            local mux, muy = require("luafemm-static-exchange").permeability(p)
            local r = f.region and m.regions[f.region]
            local circuit = r and r.circuit and circuit_ids[r.circuit][r.winding] or 0
            local turns = circuit ~= 0
                    and r.turns * f.area * o.unit ^ 2 / m.static_circuits[r.circuit].areas[f.region]
                or 1
            local sigma = p.conductivity
            if circuit ~= 0 then
                local length = axi
                        and 2 * math.pi * m.static_circuits[r.circuit].moments[f.region] / m.static_circuits[r.circuit].areas[f.region] / o.unit
                    or o.depth
                sigma = sigma * r.fill * length / (r.length > 0 and r.length or length)
            end
            d.materials[id] = {
                blockname = (p.fem_name or f.material) .. " [" .. id .. "]",
                mu_x = mux,
                mu_y = muy,
                h_c = p.hc,
                h_cangle = 0,
                j_re = circuit == 0 and f.current / 1e6 or 0,
                j_im = 0,
                sigma = sigma / 1e6,
                d_lam = 0,
                phi_h = 0,
                phi_hx = 0,
                phi_hy = 0,
                lamtype = 0,
                lamfill = 1,
                nstrands = 0,
                wired = 0,
                bh = copy(p.bh or {}),
            }
            d.labels[#d.labels + 1] = {
                f.seed[1],
                f.seed[2],
                id,
                f.mesh_size or 0,
                circuit,
                f.angle or 0,
                f.group or 0,
                turns,
                0,
                require("luafemm-static-exchange").expression(r),
            }
            labels[f.id] = #d.labels - 1
        end
    end
    return d, t, labels
end

--- Export the computational problem; this never starts a magnetic solve.
-- @tparam table m Model with its physical declarations complete.
-- @tparam string path Destination including .fem.
-- @treturn table Exported passive document.
function F.write(m, path)
    check(not m.operating, "operating-point FEM export is not implemented")
    check(m.options.frequency == 0, "harmonic FEM export is not implemented")
    local d = F.document(m)
    F.write_text(path, F.serialize(d), m.import_document and m.import_document.filename)
    return d
end
return F
