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
    "luafemm-adaptive",
    "luafemm-visual",
    "luafemm-exterior",
    "luafemm-equivalent",
    "luafemm-extended",
    "luafemm-physical",
    "luafemm-electrostatic",
    "luafemm-electric-classics",
    "luafemm-electric-fe",
    "luafemm-electric-quadrature",
    "luafemm-electric-hybrid-boundary",
    "luafemm-electric-interaction",
    "luafemm-electric-hybrid",
    "luafemm-electric-source-quadrature",
    "luafemm-electric-hybrid-view",
    "luafemm-fee",
    "luafemm-electric-native",
    "luafemm-electric-exchange",
    "luafemm-electric-metric",
    "luafemm-electric-matrix",
    "luafemm-electric-constraints",
    "luafemm-electric-periodic",
    "luafemm-electric-force",

    "luafemm-dielectrics",
    "luafemm-dielectrics-data",
    "luafemm-electric-fe-view",
    "luafemm-electric-view",
    "luafemm-electric-tex",
    "luafemm-time",
    "luafemm-time-state",
    "luafemm-time-exchange",
    "luafemm-ac-circuits",
    "luafemm-network",
    "luafemm-ee",
    "luafemm-flux-drive",
    "luafemm-airgap",
    "luafemm-age",
    "luafemm-sector-gap",
    "luafemm-sector-view",
    "luafemm-studies",
    "luafemm-response-exchange",
    "luafemm-ans",
    "luafemm-answer-reader",
    "luafemm-harmonic-exchange",
    "luafemm-axial",
    "luafemm-axisymmetric",
    "luafemm-bh",
    "luafemm-boundary",
    "luafemm-cache",
    "luafemm-complex",
    "luafemm-linear",
    "luafemm-factor",
    "luafemm-cg",
    "luafemm-components",
    "luafemm-currents",
    "luafemm-fem",
    "luafemm-frontiers",
    "luafemm-fundamental",
    "luafemm-windings",
    "luafemm-geometry",
    "luafemm-harmonic",
    "luafemm-ideal",
    "luafemm-integrals",
    "luafemm-losses",
    "luafemm-machines",
    "luafemm-materials-data",
    "luafemm-catalogue",
    "luafemm-lamination",
    "luafemm-materials",
    "luafemm-mesh",
    "luafemm-spatial",
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
