-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Bind all package modules to one installation for the process lifetime.
-- Selection never changes package.path. Missing siblings are fatal: falling
-- through to another TEXMF tree could silently combine incompatible solvers.
-- @module luafemm-loader
local L = {}
local names = {
    "luafemm-ac-circuits",
    "luafemm-ans",
    "luafemm-axial",
    "luafemm-axisymmetric",
    "luafemm-bh",
    "luafemm-boundary",
    "luafemm-cache",
    "luafemm-complex",
    "luafemm-components",
    "luafemm-currents",
    "luafemm-fem",
    "luafemm-frontiers",
    "luafemm-geometry",
    "luafemm-harmonic",
    "luafemm-ideal",
    "luafemm-integrals",
    "luafemm-losses",
    "luafemm-machines",
    "luafemm-materials-data",
    "luafemm-materials",
    "luafemm-mesh",
    "luafemm-operating",
    "luafemm-path",
    "luafemm-post",
    "luafemm-predicates",
    "luafemm-profiles",
    "luafemm-refine",
    "luafemm-shapes",
    "luafemm-sliding",
    "luafemm-static-exchange",
    "luafemm-static",
    "luafemm-tex",
    "luafemm-topology",
    "luafemm",
}
--- Install sibling-only preload functions before loading the public module.
-- @tparam string directory Directory containing the selected Lua modules.
-- @treturn string Selected directory, for provenance diagnostics.
function L.select(directory)
    local selected = package.loaded["luafemm-installation"]
    assert(not selected or selected == directory, "luafemm: conflicting installations")
    if selected then
        return selected
    end
    for _, name in ipairs(names) do
        -- Existing process-lifetime modules cannot safely be replaced.
        assert(
            not package.loaded[name],
            "luafemm: load luafemm before its private modules: " .. name
        )
        local path = directory .. "/" .. name .. ".lua"
        local file = assert(io.open(path, "rb"), "luafemm: incomplete installation: " .. path)
        file:close()
        package.preload[name] = function()
            return assert(loadfile(path))(name)
        end
    end
    package.loaded["luafemm-installation"] = directory
    return directory
end
return L
