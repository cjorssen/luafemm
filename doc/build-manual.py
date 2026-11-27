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

DOC = Path(__file__).resolve().parent
# Public sources are shallow; a development checkout may use TDS directories.
ROOT = DOC.parent if DOC.name == "doc" else DOC.parents[2]
BUILD = ROOT / "build"
RUNTIME = ROOT / "tex"
if (RUNTIME / "luatex/luafemm").is_dir():
    RUNTIME /= "luatex/luafemm"

def environment(cwd=ROOT, *, context=False):
    env = os.environ.copy()
    env["TEXMFHOME"] = os.environ.get("LUAFEMM_TEXMF", str(ROOT))
    # LuaLaTeX respects TeX's file-opening restrictions, so its cache stays
    # relative to the job directory. ConTeXt's format builder changes directory:
    # every MkIV job (including loading fragments) needs one absolute cache.
    cache = BUILD / "font-cache" if cwd == ROOT else cwd / "font-cache"
    cache.mkdir(parents=True, exist_ok=True)
    env["TEXMFCACHE"] = str(cache.relative_to(cwd))
    if context:
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
    env = environment(cwd, context=Path(command[0]).name == "texlua" and "--luatex" in command)
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w") as stream:
            result = subprocess.run(command, cwd=cwd, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT)
        if result.returncode:
            print(log.read_text(errors="replace")[-7000:])
            raise SystemExit(f"Build failed; see {log}")
    else:
        subprocess.run(command, cwd=cwd, env=env, check=True)

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
    """Inventory displayed fragments independently of downloadable documents.

    The only whole-document listings are the three format wrappers. Ordinary
    reference/tutorial examples must use actual PGF source/result blocks.
    Code-only fragments keep their separately checked, explicit contexts in
    manual-fragments.py; setup blocks execute visibly in the surrounding case.
    """
    import hashlib
    sections = DOC / "sections"
    entries = []
    cases = []
    for path in sorted(sections.glob("*.tex")):
        text = path.read_text()
        if not path.stem.startswith("tutorial") and re.search(r"\b(?:Kamil|Maya)\b", text):
            raise SystemExit(f"Narrative character outside a tutorial: {path}")
        if path.stem != "format-examples" and re.search(
                r"\\(?:luafemmsource|luafemmformatexample|luafemmlatexexample)\b", text):
            raise SystemExit(f"Unreviewed complete-source exception: {path}")
        if re.search(r"\\luafemmexample\{", text):
            raise SystemExit(f"Whole-document example outside the format appendix: {path}")
        begins = re.findall(r"\\luafemmexamplebegin\{([^}]+)\}", text)
        ends = re.findall(r"\\luafemmexampleend\{([^}]+)\}", text)
        if begins != ends:
            raise SystemExit(f"Unbalanced worked-example contexts: {path}")
        cases.extend(begins)
        contexts = list(re.finditer(r"\\luafemmexamplebegin\{([^}]+)\}(.*?)"
                                   r"\\luafemmexampleend\{\1\}", text, re.S))
        prerequisites = {}
        for i, match in enumerate(re.finditer(
                r"\\begin\{codeexample\}\[([^\]]*)\](.*?)\\end\{codeexample\}", text, re.S), 1):
            options, code = match.groups()
            if any(v in options for v in ("render instead", "hidden", "execute code=false", "from file")):
                raise SystemExit(f"Undisclosed example substitution: {path}:{i}")
            if path.stem not in ("getting-started", "tutorial") and re.search(
                    r"\\(?:documentclass|begin\{document\})", code):
                raise SystemExit(f"Document wrapper in ordinary reference: {path}:{i}")
            kind = "context fragment" if "code only" in options else "setup" if "setup code" in options else "rendered"
            context = next((c[1] for c in contexts if c.start() < match.start() < c.end()), None)
            identity = f"{path.stem}:{i}"
            entries.append({"id": identity, "kind": kind, "context": context,
                            "source": path.relative_to(ROOT).as_posix(),
                            "line": text.count("\n", 0, match.start()) + 1,
                            "requires": list(prerequisites.get(context, [])) if context else [],
                            "sha256": hashlib.sha256(code.encode()).hexdigest()})
            if context:
                prerequisites.setdefault(context, []).append(identity)
    if len(cases) != len(set(cases)):
        raise SystemExit("Duplicate worked-example context")
    for name in cases:
        if not (ROOT / "examples" / (name + ".tex")).is_file():
            raise SystemExit("Missing downloadable example: " + name)
    appendix = (sections / "format-examples.tex").read_text()
    native = re.findall(r"\\luafemmformatexample\{([^}]+)\}", appendix)
    if native != ["formats/plain", "formats/context"] or r"\luafemmlatexexample" not in appendix:
        raise SystemExit("Review the named native-format exceptions")
    inventory = {"blocks": entries, "worked_cases": cases,
        "exceptions": {"formats/plain": "plain wrapper requires its native engine",
                       "formats/context": "MkIV wrapper requires its native engine",
                       "formats/latex": "comparative loading wrapper; body executes in the manual"},
        "downloadable_documents": sorted(p.relative_to(ROOT / "examples").as_posix()
                                         for p in (ROOT / "examples").rglob("*.tex")),
        # Three wrappers plus the displayed shared scene use from-file listings.
        "expected_pgf_blocks": len(entries) + 4}
    out = BUILD / "doc"
    out.mkdir(parents=True, exist_ok=True)
    (out / "manual-examples.json").write_text(json.dumps(inventory, indent=2) + "\n")
    print(f"PASS: {len(entries)} displayed PGF blocks, {len(cases)} worked cases, three named format exceptions.")
    return inventory

def prepare_manual_examples():
    """Prepare the named format exceptions and explicit snippet contexts.

    Ordinary worked examples live directly in the section sources. Only the
    comparative LaTeX wrapper is unwrapped for execution in the format appendix.
    Two independent snippet contexts use downloadable sources; their output is
    discarded by the fragment checker, never inserted into the manual.
    """
    import hashlib
    bodies = BUILD / "doc/bodies"
    bodies.mkdir(parents=True, exist_ok=True)
    exceptions = {"formats/plain": "luatex", "formats/context": None}
    digests = {}
    for source in sorted((ROOT / "examples").rglob("*.tex")):
        name = source.relative_to(ROOT / "examples").with_suffix("").as_posix()
        # Two full bodies supply the explicit contexts of code-only snippets;
        # they are checked separately and are never used as manual figures.
        if name not in (*exceptions, "formats/latex", "sliding-band", "forces"):
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
    for folder in (ROOT / "tex", ROOT / "examples", DOC):
        files.extend(p for p in folder.rglob("*") if p.is_file() and p.suffix in (".tex", ".lua", ".sty", ".mkiv", ".py"))
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(files)}


def check_execution(log, inventory):
    """Reject omitted, repeated or out-of-order displayed PGF blocks."""
    markers = re.findall(r"LUAFEMM-MANUAL-EXAMPLE: (\d+)", log)
    if list(map(int, markers)) != list(range(1, inventory["expected_pgf_blocks"] + 1)):
        raise SystemExit("Displayed PGF example execution differs from the source inventory")


def record_build(expected_inputs):
    """Record an auditable successful build for local archive preparation."""
    current_inputs = input_digests()
    if current_inputs != expected_inputs:
        raise SystemExit("Manual inputs changed during compilation; rebuild before recording evidence")
    output = BUILD / "doc"
    prepared = json.loads((output / "example-sources.json").read_text())
    if any(current_inputs.get("examples/" + name + ".tex") != digest
           for name, digest in prepared.items()):
        raise SystemExit("Prepared examples differ from displayed sources; rebuild the manual")
    log = check_log(output / "luafemm-manual.log")
    inventory = json.loads((output / "manual-examples.json").read_text())
    check_execution(log, inventory)
    if "undefined references" in log or "Label(s) may have changed" in log:
        raise SystemExit("Cannot record an unresolved manual")
    passed = json.loads((output / "snippets/passed.json").read_text())
    if not passed['passed'] or passed['readme'] != input_digests()['README.md']:
        raise SystemExit("Fragment/README checks are stale")
    (output / "manual-inputs.json").write_text(json.dumps(input_digests(), indent=2) + "\n")


def stabilize_manual(source, out):
    """Require stable index, contents and bookmarks, including a fresh index.

    Adding the Index entry changes the .toc without necessarily producing a
    LaTeX rerun warning. Comparing only .idx or labels can accept a stale PDF.
    """
    def state():
        return {suffix: (out / ("luafemm-manual" + suffix)).read_bytes()
                if (out / ("luafemm-manual" + suffix)).exists() else b""
                for suffix in (".idx", ".toc", ".out")}
    for _ in range(5):
        before = state()
        run(["makeindex", "luafemm-manual.idx"], cwd=out, log=out / "makeindex-console.log")
        compile_tex(source, directory=out)
        log = (out / "luafemm-manual.log").read_text()
        if (before == state() and "undefined references" not in log
                and "Label(s) may have changed" not in log):
            return
    raise SystemExit("Manual references/index/contents/bookmarks did not stabilize after five passes")


def manual():
    expected_inputs = input_digests()
    check_manual_examples()
    missing = [name for name in ("pdfcrop", "gs") if shutil.which(name) is None]
    if missing:
        raise SystemExit("Manual builds need pdfcrop and Ghostscript; missing: "
                         + ", ".join(missing))
    prepare_manual_examples()
    import importlib.util
    import sys
    spec = importlib.util.spec_from_file_location('luafemm_fragments', DOC / 'manual-fragments.py')
    fragments = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fragments)
    fragments.check(sys.modules[__name__])
    out = BUILD / "doc"
    out.mkdir(parents=True, exist_ok=True)
    run(["texlua", str(DOC / "manual-catalog.lua"), str(RUNTIME)])
    source = DOC / 'luafemm-manual.tex'
    compile_tex(source, directory=out)
    stabilize_manual(source, out)
    record_build(expected_inputs)

if __name__ == "__main__":
    os.chdir(ROOT)
    manual()
