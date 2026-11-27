-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Passive material-library parser shared by generation and explicit legacy input.
-- The format is data, never Lua/TeX. Bounds and ASCII identifiers protect the
-- generic PGF choice generator; duplicate display names retain distinct IDs.
local C = {}
local function check(v, s)
    assert(v, "luafemm catalogue: " .. s)
end
local function number(s)
    local v = tonumber(s)
    check(v and v == v and math.abs(v) < math.huge, "invalid numeric value")
    return v
end
function C.parse(text)
    check(type(text) == "string" and #text <= 16 * 1024 * 1024, "input exceeds 16 MiB")
    local records, folders, used = {}, {}, {}
    local current, remaining, folder
    local function quoted(s)
        return assert(s:match('^"(.*)"$'), "luafemm catalogue: expected quoted name")
    end
    local function slug(s)
        return s:lower()
            :gsub("%%", " percent ")
            :gsub("[^a-z0-9]+", "-")
            :gsub("^%-", "")
            :gsub("%-$", "")
    end
    for line in (text .. "\n"):gmatch("(.-)\n") do
        line = line:gsub("\r", ""):match("^%s*(.-)%s*$")
        if line == "<BeginBlock>" then
            check(not current and not folder, "nested block")
            current = { fields = {}, bh = {}, path = table.concat(folders, " / ") }
            remaining = 0
        elseif line == "<EndBlock>" then
            check(current and remaining == 0 and current.name, "incomplete material")
            local id = slug(current.name)
            if used[id] then
                id = id .. "-" .. slug(folders[#folders] or "root")
            end
            check(id ~= "" and not used[id], "duplicate identifier " .. id)
            used[id] = true
            current.id = id
            records[#records + 1] = current
            check(#records <= 10000, "too many materials")
            current = nil
        elseif line == "<BeginFolder>" then
            check(not current and not folder, "nested folder header")
            folder = true
        elseif line == "<EndFolder>" then
            check(not current and #folders > 0 and not folder, "unmatched folder")
            table.remove(folders)
        elseif current and remaining > 0 then
            local b, h = line:match("^(%S+)%s+(%S+)$")
            check(b and h, "invalid B-H point")
            current.bh[#current.bh + 1] = { number(b), number(h) }
            remaining = remaining - 1
        elseif line ~= "" then
            local key, value = line:match("^<([^>]+)>%s*=%s*(.*)$")
            check(key, "invalid record")
            if key == "FolderName" then
                check(folder, "unexpected folder name")
                folders[#folders + 1] = quoted(value)
                folder = false
            elseif key == "FolderURL" or key == "FolderVendor" then
                check(not current, "unexpected folder metadata")
            elseif current and key == "BlockName" then
                current.name = quoted(value)
            elseif current then
                check(
                    key:match("^[A-Za-z_][A-Za-z0-9_]*$") and current.fields[key] == nil,
                    "invalid/duplicate field"
                )
                current.fields[key] = number(value)
                if key == "BHPoints" then
                    remaining = current.fields[key]
                    check(
                        remaining >= 0 and remaining % 1 == 0 and remaining <= 100000,
                        "invalid B-H count"
                    )
                end
            else
                check(false, "unexpected property " .. key)
            end
        end
    end
    check(not current and not folder and #folders == 0 and #records > 0, "incomplete database")
    return records
end
function C.read(path)
    check(
        type(path) == "string" and path:match("%S") and not path:find("%c"),
        "invalid catalogue filename"
    )
    local f = assert(io.open(path, "rb"), "luafemm catalogue: cannot open " .. path)
    local text = f:read(16 * 1024 * 1024 + 1)
    f:close()
    return C.parse(text), text
end
--- Portable byte fingerprint for provenance, not a cryptographic signature.
function C.fingerprint(text)
    local a, b = 1, 0
    for i = 1, #text do
        a = (a + text:byte(i)) % 65521
        b = (b + a) % 65521
    end
    return string.format("%08x:%d", b * 65536 + a, #text)
end
return C
