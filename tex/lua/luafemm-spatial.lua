-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Conservative broad-phase box queries for normalised meshing coordinates.
-- This index never makes a geometric decision: callers retain their distance
-- tests and exact-sign predicates. Results have ascending integer IDs, preserving
-- exhaustive-search tie breaking. Large boxes live in an overflow set instead
-- of allocating an unbounded number of grid entries. No third-party code copied.
-- @module luafemm-spatial
local S = {}
local floor, min, max = math.floor, math.min, math.max

--- Create a mutable uniform-grid index.
-- @tparam number size Positive cell width in the caller's coordinate system.
-- @return Index with insert(id,xmin,ymin,xmax,ymax), remove(id), query(box).
function S.new(size)
    assert(size > 0 and size < math.huge, "luafemm spatial: invalid cell width")
    local cells, records, overflow = {}, {}, {}
    local index = { queries = 0, candidates = 0 }
    local function bounds(x0, y0, x1, y1)
        local magnitude = max(1, math.abs(x0), math.abs(y0), math.abs(x1), math.abs(y1))
        assert(
            x0 == x0 and y0 == y0 and x1 == x1 and y1 == y1 and magnitude < math.huge,
            "luafemm spatial: non-finite search bounds"
        )
        assert(x0 <= x1 and y0 <= y1, "luafemm spatial: reversed search bounds")
        -- Outward padding covers roundoff in normalised box/disk construction.
        -- It is deliberately independent of mesh tolerance and never merges points.
        local pad = 320 * 2 ^ -52 * magnitude
        assert(
            math.abs(x0 - pad) < math.huge
                and math.abs(y0 - pad) < math.huge
                and math.abs(x1 + pad) < math.huge
                and math.abs(y1 + pad) < math.huge,
            "luafemm spatial: search bounds overflow; review geometry scale"
        )
        return { x0 - pad, y0 - pad, x1 + pad, y1 + pad }
    end
    local function range(b)
        local x0, y0, x1, y1 =
            floor(b[1] / size), floor(b[2] / size), floor(b[3] / size), floor(b[4] / size)
        local n = (x1 - x0 + 1) * (y1 - y0 + 1)
        -- A cell quotient can overflow even though all supplied bounds are
        -- finite. Never enter a numeric for-loop with infinite endpoints.
        if
            n ~= n
            or n == math.huge
            or math.abs(x0) == math.huge
            or math.abs(y0) == math.huge
            or math.abs(x1) == math.huge
            or math.abs(y1) == math.huge
        then
            return 0, 0, 0, 0, math.huge
        end
        return x0, y0, x1, y1, n
    end
    function index.remove(id)
        local r = records[id]
        if not r then
            return
        end
        for _, k in ipairs(r.cells) do
            local bucket = cells[k]
            bucket[id] = nil
            if not next(bucket) then
                cells[k] = nil
            end
        end
        records[id], overflow[id] = nil, nil
    end
    function index.insert(id, x0, y0, x1, y1)
        index.remove(id)
        local b = bounds(x0, y0, x1, y1)
        local r = { box = b, cells = {} }
        records[id] = r
        local a, c, d, e, n = range(b)
        if n > 64 then
            overflow[id] = true
            return
        end
        for x = a, d do
            for y = c, e do
                local k = x .. ":" .. y
                cells[k] = cells[k] or {}
                cells[k][id] = true
                r.cells[#r.cells + 1] = k
            end
        end
    end
    function index.query(x0, y0, x1, y1)
        index.queries = index.queries + 1
        local b = bounds(x0, y0, x1, y1)
        local found, ids = {}, {}
        local function consider(id)
            if found[id] then
                return
            end
            found[id] = true
            local q = records[id].box
            if q[1] <= b[3] and q[3] >= b[1] and q[2] <= b[4] and q[4] >= b[2] then
                ids[#ids + 1] = id
            end
        end
        local a, c, d, e, n = range(b)
        if n > 4096 then
            for id in pairs(records) do
                consider(id)
            end
        else
            for x = a, d do
                for y = c, e do
                    for id in pairs(cells[x .. ":" .. y] or {}) do
                        consider(id)
                    end
                end
            end
            for id in pairs(overflow) do
                consider(id)
            end
        end
        table.sort(ids)
        index.candidates = index.candidates + #ids
        return ids
    end
    return index
end

--- Insert the bounding box of a segment, optionally expanded by a distance.
function S.segment(index, id, a, b, pad)
    pad = pad or 0
    index.insert(
        id,
        min(a[1], b[1]) - pad,
        min(a[2], b[2]) - pad,
        max(a[1], b[1]) + pad,
        max(a[2], b[2]) + pad
    )
end

--- Insert a box containing a segment's open diametral disk.
function S.disk(index, id, a, b)
    local dx, dy = b[1] - a[1], b[2] - a[2]
    local x, y = a[1] + dx / 2, a[2] + dy / 2
    local r = math.sqrt(dx * dx + dy * dy) / 2
    index.insert(id, x - r, y - r, x + r, y + r)
end
return S
