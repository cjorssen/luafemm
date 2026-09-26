#!/usr/bin/env python3
# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
"""Local development builds; runtime documents do not require Python.

All TeX outputs stay in build/. Runtime files are found through TEXMFHOME,
so tests exercise the installed layout rather than a flat working directory.
ConTeXt is explicitly run on LuaTeX (MkIV), never silently on LuaMetaTeX.
"""
from pathlib import Path
import argparse
import json
import os
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent.parent
BUILD = Path(os.environ.get("LUAFEMM_BUILD", ROOT / "build")).resolve()
RUNTIME = ROOT / "tex" / "generic" / "luafemm"


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
    env.pop("LUAINPUTS", None)
    env.pop("TEXINPUTS", None)
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


def formats():
    snapshots = {}
    for name, engine in (("plain", "luatex"), ("latex", "lualatex"), ("context", None)):
        folder = BUILD / "formats" / name
        folder.mkdir(parents=True, exist_ok=True)
        # Some ConTeXt launchers return zero after a format/TeX failure.
        # Require fresh outputs rather than accepting an older successful PDF.
        for extension in ("log", "pdf"):
            (folder / f"{name}.{extension}").unlink(missing_ok=True)
        for filename in (name + ".tex", "scene.tex"):
            shutil.copy2(ROOT / "tests" / "formats" / filename, folder / filename)
        if engine:
            cmd = [engine, "--no-shell-escape", "-interaction=nonstopmode", "-halt-on-error", name + ".tex"]
        else:
            cmd = context_command() + [name + ".tex"]
        console = folder / "console.log"
        run(cmd, cwd=folder, log=console)
        log = check_log(folder / (name + ".log"))
        if not (folder / (name + ".pdf")).is_file():
            raise SystemExit(f"No PDF produced by {name}")
        snapshots[name] = {k: float(v) for k, v in re.findall(r"LUAFEMM-SAMPLE-(\w+): ([+\d.eE-]+)", log)}
        if len(snapshots[name]) != 3:
            raise SystemExit(f"Missing numerical assertions in {name}")
        if name == "context":
            output = console.read_text(errors="replace")
            if "engine: lua 5.3" not in output or "used engine: luatex" not in output:
                raise SystemExit("The ConTeXt test did not run on LuaTeX/MkIV")
    reference = snapshots["plain"]
    for name, values in snapshots.items():
        for key, expected in reference.items():
            if abs(values[key] - expected) > 1e-10:
                raise SystemExit(f"Format mismatch: {name} {key}: {values[key]} != {expected}")
    (BUILD / "formats" / "results.json").write_text(json.dumps(snapshots, indent=2) + "\n")
    print("PASS: plain, LaTeX and ConTeXt MkIV agree to 1e-10.")


def tests():
    for path in sorted((ROOT / "tests").glob("test_*.lua")):
        run(["texlua", str(path)], log=BUILD / "tests" / (path.stem + ".log"))
    run(["python3", "tests/check_predicates.py"], log=BUILD / "tests" / "predicates.log")
    for path in (ROOT / "tests").glob("tikz-*.tex"):
        compile_tex(path, "luatex", BUILD / "tests")
    formats()


def examples():
    for path in sorted((ROOT / "examples").glob("*.tex")):
        compile_tex(path, "luatex" if path.stem == "plain" else "lualatex")
    for name, engine in (("plain", "luatex"), ("latex", "lualatex"), ("context", None)):
        folder = BUILD / "examples/formats" / name
        folder.mkdir(parents=True, exist_ok=True)
        for filename in (name + ".tex", "scene.tex"):
            shutil.copy2(ROOT / "examples/formats" / filename, folder / filename)
        for extension in ("log", "pdf"):
            (folder / f"{name}.{extension}").unlink(missing_ok=True)
        command = ([engine, "--no-shell-escape", "-interaction=nonstopmode", "-halt-on-error"]
                   if engine else context_command())
        run(command + [name + ".tex"], cwd=folder, log=folder / "console.log")
        check_log(folder / (name + ".log"))
        if not (folder / (name + ".pdf")).is_file():
            raise SystemExit(f"No PDF produced by the {name} example")


def manual():
    # The tutorial embeds the exact output of its independently usable sources.
    # Compile each model once here, not once per manual/index pass.
    for name in ("coil", "alnico"):
        compile_tex(ROOT / "examples" / f"tutorial-{name}.tex")
    out = BUILD / "doc"
    out.mkdir(parents=True, exist_ok=True)
    run(["texlua", "scripts/manual-catalog.lua"])
    source = ROOT / "doc" / "luafemm-manual.tex"
    compile_tex(source, directory=out)
    run(["makeindex", "luafemm-manual.idx"], cwd=out, log=out / "makeindex-console.log")
    compile_tex(source, directory=out)
    compile_tex(source, directory=out)
    log = (out / "luafemm-manual.log").read_text()
    if "undefined references" in log or "Label(s) may have changed" in log:
        raise SystemExit("Manual references need another pass")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("test", "formats", "examples", "manual", "all"))
    args = parser.parse_args()
    BUILD.mkdir(exist_ok=True)
    targets = {"test": tests, "formats": formats, "examples": examples, "manual": manual}
    if args.target == "all":
        for function in (tests, examples, manual):
            function()
    else:
        targets[args.target]()
