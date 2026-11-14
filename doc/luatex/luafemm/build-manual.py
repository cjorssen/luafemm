#!/usr/bin/env python3
# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
"""Rebuild the manual and execute its displayed examples without shell escape.

Run from the browsable source archive root. This is a documentation builder,
not a runtime dependency or a general release/test tool. Generated fragments,
private caches, native-format checks and the PDF are written under build/.
"""
from pathlib import Path
import json
import os
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[3]
BUILD = ROOT / "build"
RUNTIME = ROOT / "tex/luatex/luafemm"

def environment(cwd=ROOT):
    env = os.environ.copy()
    env["TEXMFHOME"] = os.environ.get("LUAFEMM_TEXMF", str(ROOT))
    cache = BUILD / "font-cache" if cwd == ROOT else cwd / "font-cache"
    cache.mkdir(parents=True, exist_ok=True)
    # Relative paths also work under restrictive TeX file-opening rules.
    env["TEXMFCACHE"] = str(cache.relative_to(cwd))
    if cwd.name == "context":
        # The MkIV format builder changes directory; its cache must be absolute.
        cache = BUILD / "context-cache"
        cache.mkdir(parents=True, exist_ok=True)
        env["TEXMFCACHE"] = str(cache)
    env.pop("LUA_PATH", None)
    selected = Path(env["TEXMFHOME"]) / "tex"
    env["LUAINPUTS"] = str(selected) + "//:"
    env["TEXINPUTS"] = str(selected) + "//:"
    return env

def run(command, cwd=ROOT, log=None):
    print("+", " ".join(map(str, command)), flush=True)
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w") as stream:
            result = subprocess.run(command, cwd=cwd, env=environment(cwd), stdout=stream,
                                    stderr=subprocess.STDOUT)
        if result.returncode:
            print(log.read_text(errors="replace")[-7000:])
            raise SystemExit(f"Build failed; see {log}")
    else:
        subprocess.run(command, cwd=cwd, env=environment(cwd), check=True)

def compile_tex(source, engine="lualatex", directory=None):
    directory = directory or BUILD / "examples"
    directory.mkdir(parents=True, exist_ok=True)
    run([engine, "--no-shell-escape", "-interaction=nonstopmode", "-halt-on-error",
         f"-output-directory={directory}", str(source)], log=directory / (source.stem + "-console.log"))
    check_log(directory / (source.stem + ".log"))

def check_log(log):
    if not log.is_file():
        raise SystemExit(f"No TeX log produced; inspect {log.parent / 'console.log'}")
    text = log.read_text(errors="replace")
    if re.search(r"^!|Missing character:|Overfull \\[hv]box|tex error\s+>", text, re.MULTILINE):
        raise SystemExit(f"TeX diagnostics require attention: {log}")
    return text

def context_command():
    # TeX Live's context/mtxrun launcher may itself be LuaMetaTeX. Using the
    # distributed mtxrun.lua under texlua guarantees a MkIV invocation and
    # also avoids launcher-script lookup failures on some TeX Live installs.
    scripts = []
    for name in ("mtxrun.lua", "mtx-context.lua"):
        path = subprocess.check_output(["kpsewhich", "-format=texmfscripts", name], text=True).strip()
        if not path:
            raise SystemExit("ConTeXt MkIV scripts not found; install ConTeXt first")
        scripts.append(path)
    return ["texlua", scripts[0], "--script", scripts[1], "--luatex",
            "--batchmode", "--nonstopmode", "--once"]

def check_manual_examples():
    """Every distributed example must appear with its source and actual output.

    scene.tex is shared input, not a standalone document; its three wrappers
    each have their own source/output pair. Fail on omissions and stale names
    so adding an example cannot silently leave the manual incomplete.
    """
    text = "\n".join(p.read_text() for p in (ROOT / "doc").rglob("*.tex"))
    complete = set(re.findall(r"\\luafemm(?:format)?example\{([a-z0-9/-]+)\}", text))
    shared = set(re.findall(r"\\luafemmsource\{([a-z0-9/-]+)\}", text))
    actual = {p.relative_to(ROOT / "examples").with_suffix("").as_posix()
              for p in (ROOT / "examples").rglob("*.tex")}
    expected_shared = {"formats/scene"}
    if shared != expected_shared or complete != actual - expected_shared:
        raise SystemExit("Manual example coverage mismatch: "
                         f"missing {sorted(actual - complete - shared)}, "
                         f"stale {sorted((complete | shared) - actual)}, "
                         f"shared listings {sorted(shared)}")
    print(f"PASS: {len(complete)} complete source/output pairs and the shared scene are documented.")

def prepare_manual_examples():
    """Execute all LaTeX bodies inside the manual, from the displayed sources.

    Only three native-format wrappers are exceptions. Generated fragments omit
    document wrappers and preamble-only package loading, never physical code.
    Their exact source digests are retained for auditing the transformation.
    """
    import hashlib
    bodies = BUILD / "doc/bodies"
    bodies.mkdir(parents=True, exist_ok=True)
    exceptions = {"plain": "luatex", "formats/plain": "luatex", "formats/context": None}
    digests = {}
    for source in sorted((ROOT / "examples").rglob("*.tex")):
        name = source.relative_to(ROOT / "examples").with_suffix("").as_posix()
        if name == "formats/scene":
            continue
        text = source.read_text()
        digests[name] = hashlib.sha256(source.read_bytes()).hexdigest()
        if name in exceptions:
            folder = BUILD / "doc/format-checks" / name
            folder.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, folder / source.name)
            shutil.copy2(ROOT / "examples/formats/scene.tex", folder / "scene.tex")
            engine = exceptions[name]
            command = ([engine, "--no-shell-escape", "-interaction=nonstopmode", "-halt-on-error"]
                       if engine else context_command())
            pdf = folder / (source.stem + ".pdf")
            pdf.unlink(missing_ok=True)
            run(command + [source.name], cwd=folder, log=folder / "console.log")
            check_log(folder / (source.stem + ".log"))
            if not pdf.is_file():
                raise SystemExit(f"Missing native-format output: {name}")
            output = BUILD / "doc/examples" / (name + ".pdf")
            output.parent.mkdir(parents=True, exist_ok=True)
            run(["pdfcrop", "--luatex", "--margins", "3", str(pdf), str(output)],
                log=output.with_suffix(".log"))
            continue
        preamble, body = text.split(r"\begin{document}", 1)
        body, tail = body.split(r"\end{document}", 1)
        if tail.strip():
            raise SystemExit(f"Unexpected trailing example source: {source}")
        preamble = re.sub(r"^\\(?:documentclass|usepackage)(?:\[[^\n]*?\])?\{[^\n]*?\}[^\n]*$", "", preamble, flags=re.M)
        preamble = re.sub(r"^\\pagestyle\{[^}]+\}[^\n]*$", "", preamble, flags=re.M)
        fragment = preamble + "\n" + body
        fragment = fragment.replace(r"\input scene.tex", r"\input examples/formats/scene.tex")
        output = bodies / (name + ".tex")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(fragment)
    (BUILD / "doc/example-sources.json").write_text(json.dumps(digests, indent=2) + "\n")

def input_digests():
    """Digest all executable/documentation inputs used to build the manual."""
    import hashlib
    files = [ROOT / "README.md"]
    for folder in (ROOT / "tex", ROOT / "examples", ROOT / "doc/luatex/luafemm"):
        files.extend(p for p in folder.rglob("*") if p.is_file() and p.suffix in (".tex", ".lua", ".sty", ".mkiv", ".py"))
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(files)}


def record_build():
    """Record an auditable successful build for local archive preparation."""
    output = BUILD / "doc"
    log = check_log(output / "luafemm-manual.log")
    if "undefined references" in log or "Label(s) may have changed" in log:
        raise SystemExit("Cannot record an unresolved manual")
    passed = json.loads((output / "snippets/passed.json").read_text())
    if not passed['passed'] or passed['readme'] != input_digests()['README.md']:
        raise SystemExit("Fragment/README checks are stale")
    (output / "manual-inputs.json").write_text(json.dumps(input_digests(), indent=2) + "\n")


def manual():
    check_manual_examples()
    missing = [name for name in ("pdfcrop", "gs") if shutil.which(name) is None]
    if missing:
        raise SystemExit("Manual builds need pdfcrop and Ghostscript; missing: "
                         + ", ".join(missing))
    prepare_manual_examples()
    import importlib.util
    import sys
    spec = importlib.util.spec_from_file_location('luafemm_fragments', ROOT / 'doc/luatex/luafemm/manual-fragments.py')
    fragments = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fragments)
    fragments.check(sys.modules[__name__])
    out = BUILD / "doc"
    out.mkdir(parents=True, exist_ok=True)
    run(["texlua", "doc/luatex/luafemm/manual-catalog.lua"])
    source = ROOT / 'doc/luatex/luafemm/luafemm-manual.tex'
    compile_tex(source, directory=out)
    # New chapters can change both the table of contents and index pagination.
    # Require stable index input as well as stable cross-references; do not keep
    # the page numbers from the first pass after the document has grown.
    for _ in range(5):
        index = (out / "luafemm-manual.idx").read_bytes()
        run(["makeindex", "luafemm-manual.idx"], cwd=out, log=out / "makeindex-console.log")
        compile_tex(source, directory=out)
        log = (out / "luafemm-manual.log").read_text()
        if (index == (out / "luafemm-manual.idx").read_bytes()
                and "undefined references" not in log and "Label(s) may have changed" not in log):
            break
    else:
        raise SystemExit("Manual references/index did not stabilize after five passes")
    record_build()

if __name__ == "__main__":
    os.chdir(ROOT)
    manual()
