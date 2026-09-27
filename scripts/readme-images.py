#!/usr/bin/env python3
# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
"""Render README previews from the existing, executable TikZ examples.

Requires the usual documentation TeX tools and Poppler's pdftoppm. Temporary
standalone wrappers and PDFs stay in build/readme; tracked PNGs have a white
background so labels remain legible in either GitHub theme. No solver or
geometry is duplicated here.
"""
import shutil

import build


def main():
    if not shutil.which("pdftoppm"):
        raise SystemExit("Install Poppler's pdftoppm to render README previews")
    output = build.ROOT / "docs/images"
    output.mkdir(parents=True, exist_ok=True)
    work = build.BUILD / "readme"
    work.mkdir(parents=True, exist_ok=True)
    examples = (
        ("curved-core", "tilted-u", "\n\\begin{center}"),
        ("tutorial-machine-sine", "slotted-machine-field", "\\begin{document}"),
        ("machine-nodes", "slotted-rotor", None),
    )
    for example, name, settings_end in examples:
        source = (build.ROOT / "examples" / (example + ".tex")).read_text()
        start = source.index(r"\begin{tikzpicture}")
        stop = source.index(r"\end{tikzpicture}", start) + len(r"\end{tikzpicture}")
        settings = ""
        if settings_end:
            first = source.index(r"\tikzset{")
            settings = source[first:source.index(settings_end, first)]
        wrapper = work / (name + ".tex")
        wrapper.write_text(r"""\documentclass[tikz,border=3mm]{standalone}
\usepackage{fontspec}
\usepackage{tikz,pgfplots}
\usetikzlibrary{femm,calc}
\usepgfplotslibrary{femm}
\pgfplotsset{compat=1.18}
""" + settings + "\n\\begin{document}\n" + source[start:stop]
                           + "\n\\end{document}\n")
        build.compile_tex(wrapper, directory=work)
        build.run(["pdftoppm", "-singlefile", "-scale-to", "1280", "-png",
                   str(wrapper.with_suffix(".pdf")), str(output / name)],
                  log=work / (name + "-render.log"))
        print("Rendered", output / (name + ".png"), flush=True)


if __name__ == "__main__":
    main()
