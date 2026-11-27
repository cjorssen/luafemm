-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
--- Bounded FEMM 5 real planar transient answer interchange.
-- DeltaT and A_previous describe time history, never a permeability bias.
-- Native answers have no absolute clock or Lua motion template; the caller
-- supplies time explicitly and uses .lft for a complete moving-model restart.
local X = {}
local function check(v, s)
    assert(v, "luafemm transient exchange: " .. s)
end
function X.validate(m)
    require("luafemm-time").eligible(m)
    require("luafemm-time").validate(m)
    check(
        m.options.problem_type == "planar" and not m.sliding,
        "native transient answers require fixed planar geometry; use a time-state restart"
    )
    for _, r in ipairs(m.regions) do
        check(
            r.current_constraint ~= "zero",
            "isolated conductor export needs native circuit voltage records; use a time-state restart"
        )
    end
end
--- Import an ordinary time answer through the validated static mesh reader,
-- then attach the distinct owned time derivative and rebuild instantaneous J.
function X.model(a, options)
    options = options or {}
    local time = options.time
    check(
        type(time) == "number" and time == time and time < math.huge and time >= a.dt,
        "supply import time >= DeltaT in seconds; the answer has no absolute clock"
    )
    check(#a.ages == 0, "native transient AGE snapshots are unsupported")
    for _, s in ipairs(a.states) do
        check(s[1] == 1, "transient voltage-gradient conductor snapshots are unsupported")
    end
    local F = require("luafemm-fem")
    local b = F.copy(a)
    b.dt, b.time_previous = nil, nil
    b.document.header.deltat = nil
    b.document.header.prevsoln = ""
    for _, p in ipairs(b.nodes) do
        p[5] = nil
    end
    local overrides = F.copy(options)
    overrides.time = nil
    local M = require("luafemm")
    local m = require("luafemm-answer-reader").model(b, overrides)
    require("luafemm-time").eligible(m)
    local previous = F.copy(m)
    previous.A = F.copy(a.time_previous)
    previous.import_answer = nil
    M.update_fields(previous)
    local old = require("luafemm-time-state").capture(previous, { time = time - a.dt })
    local d = old.data
    for _, pair in ipairs(a.pairs) do
        local x, y = d.A[pair[1] + 1], d.A[pair[2] + 1]
        check(
            math.abs(x - (pair[3] == 0 and y or -y))
                <= 1e-7 * math.max(1e-14, math.abs(x), math.abs(y)),
            "previous time potential violates periodic pairing"
        )
    end
    m.transient = {
        dt = a.dt,
        time = time,
        previous = F.copy(a.time_previous),
        linkage = d.linkage,
        energy = d.energy,
        source_scale = 1,
        boundary_scale = 1,
        motion = "fixed",
    }
    m.time_sources = d.sources
    m.time_mesh = { nodes = d.nodes, elements = d.elements }
    -- Preserve time provenance in the snapshot signature, not as a DC bias.
    m.import_answer.time_previous = F.copy(a.time_previous)
    m.import_answer.dt = a.dt
    M.update_fields(m)
    m.transient_contract = require("luafemm-cache").signature(m)
    return m
end
return X
