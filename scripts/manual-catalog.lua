-- SPDX-License-Identifier: LPPL-1.3c
-- Copyright (C) 2026 Christophe Jorssen
-- Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
-- LPPL maintenance status: maintained
-- This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-- Generate the typeset catalogue from the same records used at runtime.
dofile("tests/bootstrap.lua")
local out = assert(io.open("build/doc/materials-catalog.tex", "w"))
local function escape(text)
    return (text:gsub("[%%&#_$]", function(character)
        return "\\" .. character
    end))
end
out:write(
    "\\begingroup\\small\n\\begin{longtable}{p{.51\\linewidth}p{.39\\linewidth}}\n",
    "\\toprule Identifier & Original FEMM name \\\\ \\midrule\\endhead\n"
)
for _, record in ipairs(require("luafemm-materials").list()) do
    out:write("\\nolinkurl{", record.id, "} & ", escape(record.name), " \\\\\n")
end
out:write("\\bottomrule\n\\end{longtable}\n\\endgroup\n")
out:close()
