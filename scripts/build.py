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
        for filename in (name + ".tex", "scene.tex", "components.tex", "boundaries.tex"):
            shutil.copy2(ROOT / "tests" / "formats" / filename, folder / filename)
        if engine:
            cmd = [engine, "--no-shell-escape", "-interaction=nonstopmode", "-halt-on-error", name + ".tex"]
        else:
            cmd = context_command() + [name + ".tex"]
        cache_file = folder / "format-cache.lfc"
        cache_file.unlink(missing_ok=True)
        component_file = folder / "component-cache.lfc"
        component_file.unlink(missing_ok=True)
        boundary_file = folder / "boundary-cache.lfc"
        boundary_file.unlink(missing_ok=True)
        passes = {}
        for stage, mode, expected in (("cold", "auto", "miss"),
                                      ("warm", "auto", "solution"),
                                      ("frozen", "frozen", "solution")):
            # Warm runs change only presentation and profile sampling. Frozen
            # runs consume the cache originally written by plain LuaTeX.
            if stage == "frozen":
                cache_file.write_bytes(shared_cache)
                component_file.write_bytes(shared_component)
                boundary_file.write_bytes(shared_boundary)
            settings = {"CacheMode": mode, "ExpectedCache": expected,
                        "CacheScale": "1" if stage == "cold" else "2",
                        "CacheColor": "blue" if stage == "cold" else "red",
                        "CacheLines": "9" if stage == "cold" else "11",
                        "CacheSamples": "17" if stage == "cold" else "23"}
            (folder / "cache-settings.tex").write_text("".join(
                "\\def\\" + key + "{" + value + "}\n" for key, value in settings.items()))
            # Some ConTeXt launchers return zero after a format/TeX failure.
            # Every pass must produce fresh outputs, including frozen reads.
            for extension in ("log", "pdf"):
                (folder / f"{name}.{extension}").unlink(missing_ok=True)
            console = folder / (stage + "-console.log")
            run(cmd, cwd=folder, log=console)
            log = check_log(folder / (name + ".log"))
            shutil.copy2(folder / (name + ".log"), folder / (stage + ".log"))
            if not (folder / (name + ".pdf")).is_file():
                raise SystemExit(f"No PDF produced by {name}/{stage}")
            if f"LUAFEMM-CACHE: {expected}" not in log:
                raise SystemExit(f"Wrong cache outcome in {name}/{stage}")
            if stage != "cold" and "luafemm: Newton" in log:
                raise SystemExit(f"Cache hit reran the solver in {name}/{stage}")
            values = {k: float(v) for k, v in re.findall(
                r"LUAFEMM-SAMPLE-(\w+): ([+\d.eE-]+)", log)}
            if len(values) != 3:
                raise SystemExit(f"Missing numerical assertions in {name}/{stage}")
            component = re.search(r"LUAFEMM-COMPONENT: ([+\d.eE-]+)", log)
            if not component:
                raise SystemExit(f"Missing component result in {name}/{stage}")
            values["COMPONENT"] = float(component.group(1))
            boundary = re.search(r"LUAFEMM-BOUNDARY: ([+\d.eE-]+)", log)
            if not boundary:
                raise SystemExit(f"Missing boundary result in {name}/{stage}")
            values["BOUNDARY"] = float(boundary.group(1))
            passes[stage] = values
            if name == "plain" and stage == "cold":
                shared_cache = cache_file.read_bytes()
                shared_component = component_file.read_bytes()
                shared_boundary = boundary_file.read_bytes()
            if stage == "frozen" and (cache_file.read_bytes() != shared_cache
                                      or component_file.read_bytes() != shared_component
                                      or boundary_file.read_bytes() != shared_boundary):
                raise SystemExit("Frozen mode modified its input cache")
            if name == "context":
                output = console.read_text(errors="replace")
                if "engine: lua 5.3" not in output or "used engine: luatex" not in output:
                    raise SystemExit("The ConTeXt test did not run on LuaTeX/MkIV")
        snapshots[name] = passes["cold"]
        for stage in ("warm", "frozen"):
            for key, expected in passes["cold"].items():
                if abs(passes[stage][key] - expected) > 1e-10:
                    raise SystemExit(f"Cached result mismatch: {name}/{stage}/{key}")
    reference = snapshots["plain"]
    for name, values in snapshots.items():
        for key, expected in reference.items():
            if abs(values[key] - expected) > 1e-10:
                raise SystemExit(f"Format mismatch: {name} {key}: {values[key]} != {expected}")
    (BUILD / "formats" / "results.json").write_text(json.dumps(snapshots, indent=2) + "\n")
    print("PASS: cold, warm and shared frozen caches agree in all three formats to 1e-10.")


def invalid_components():
    """Reject unsupported placements and invalid declarations before solving."""
    folder = BUILD / "tests"
    folder.mkdir(parents=True, exist_ok=True)
    expected = ("require leg thickness", "is unavailable on u core",
                "component inside matrix", "declare components before meshing",
                "Unknown luafemm material")
    for case, message in enumerate(expected):
        source = folder / f"component-invalid-{case}.tex"
        source.write_text("\\def\\InvalidComponentCase{" + str(case) + "}\n"
                          + "\\input{" + str(ROOT / "tests/components-invalid.tex") + "}\n")
        command = ["luatex", "--no-shell-escape", "-interaction=nonstopmode",
                   "-halt-on-error", f"-output-directory={folder}", str(source)]
        result = subprocess.run(command, cwd=ROOT, env=environment(),
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (folder / (source.stem + "-console.log")).write_text(result.stdout)
        if result.returncode == 0 or message not in result.stdout:
            raise SystemExit(f"Expected component rejection missing: {message}")
    print("PASS: invalid component geometry, anchors, materials and placement are rejected.")


def tests():
    for path in sorted((ROOT / "tests").glob("test_*.lua")):
        run(["texlua", str(path)], log=BUILD / "tests" / (path.stem + ".log"))
    run(["python3", "tests/check_predicates.py"], log=BUILD / "tests" / "predicates.log")
    for path in (ROOT / "tests").glob("tikz-*.tex"):
        compile_tex(path, "luatex", BUILD / "tests")
    invalid_components()
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


def check_manual_examples():
    """Every distributed example must appear with its source and actual output.

    scene.tex is shared input, not a standalone document; its three wrappers
    each have their own source/output pair. Fail on omissions and stale names
    so adding an example cannot silently leave the manual incomplete.
    """
    text = "\n".join(p.read_text() for p in (ROOT / "doc").rglob("*.tex"))
    complete = set(re.findall(r"\\luafemmexample\{([a-z0-9/-]+)\}", text))
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


def crop_manual_examples():
    """Remove paper margins, preserving every page and its vector contents.

    The original example PDFs remain untouched. pdfcrop and Ghostscript are
    documentation-build tools only; runtime documents do not invoke them.
    """
    for source in sorted((ROOT / "examples").rglob("*.tex")):
        name = source.relative_to(ROOT / "examples").with_suffix("")
        if name.as_posix() == "formats/scene":
            continue
        original = BUILD / "examples" / name
        if name.parts[0] == "formats":
            original = original / name.name
        output = BUILD / "doc/examples" / name.with_suffix(".pdf")
        output.parent.mkdir(parents=True, exist_ok=True)
        run(["pdfcrop", "--luatex", "--margins", "3",
             str(original.with_suffix(".pdf")), str(output)],
            log=output.with_suffix(".log"))


def manual():
    check_manual_examples()
    missing = [name for name in ("pdfcrop", "gs") if shutil.which(name) is None]
    if missing:
        raise SystemExit("Manual builds need pdfcrop and Ghostscript; missing: "
                         + ", ".join(missing))
    # Compile every complete example once, then embed its actual output and
    # source with PGF's documentation machinery. Do not solve large models on
    # every manual/index pass or depend on a PDF left by an older build.
    examples()
    crop_manual_examples()
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
        for function in (tests, manual):
            function()
    else:
        targets[args.target]()
