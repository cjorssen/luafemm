-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Optional independent ANS postprocessor fixtures; invoked by post-interop.py.
dofile("tests/bootstrap.lua")
local f, post = require("luafemm"), require("luafemm-post")
local make = dofile("tests/u_model.lua")
local out = assert(io.open("build/interop/post.tsv", "w"))
for _, h in ipairs({ 4, 2 }) do
    local m = make(h, 75, 3000)
    m.options.depth = 20
    m.options.interpolation = "femm"
    m.options.mesher = "delaunay"
    for id, p in pairs(m.materials) do
        p.interpolation = "femm"
        f.material(m, id, p)
    end
    post.select(
        m,
        "armature",
        { material = "iron", contours = { { { -41, 32 }, { 41, 32 }, { 41, 49 }, { -41, 49 } } } }
    )
    local r = post.measure(m, { name = "reference", selection = "armature", quantity = "force" })
    local file = "build/interop/post-" .. h .. ".ans"
    f.export_ans(m, file)
    out:write(
        string.format(
            "%s\t0\t40\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
            file,
            r.energy,
            r.coenergy,
            r.Fx,
            r.Fy,
            r.torque
        )
    )
end
out:close()
