#!/usr/bin/env python3
# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
"""Prepare local source/TDS archives and test the extracted installation.

No VCS, network or publishing operation is performed. ZIP entry order, dates
and permissions are fixed; embedded PDFs are those produced by the local build.
"""
from pathlib import Path
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
import zipfile

import build

ROOT = build.ROOT
OUT = build.BUILD / "dist"


def version():
    loader = (build.RUNTIME / "luafemm.tex").read_text()
    value = re.search(r"\\def\\luafemmversion\{([^}]+)\}", loader).group(1)
    date = re.search(r"\\def\\luafemmversiondate\{([^}]+)\}", loader).group(1)
    lua = (build.RUNTIME / "luafemm.lua").read_text()
    wrapper = (ROOT / "tex/latex/luafemm/luafemm.sty").read_text()
    if f'version = "{value}"' not in lua or value not in wrapper:
        raise SystemExit("Inconsistent Lua/TeX release versions")
    for filename in ("README.md", "CHANGELOG.md"):
        if value not in (ROOT / filename).read_text():
            raise SystemExit(f"Missing release version in {filename}")
    return value, tuple(map(int, date.split("-"))) + (0, 0, 0)


def verify_catalogue():
    folder = OUT / "regenerated"
    folder.mkdir(parents=True, exist_ok=True)
    targets = (folder / "materials.lua", folder / "catalog.md")
    build.run(["texlua", "scripts/import-materials.lua", "data/femm/matlib.dat", *map(str, targets)])
    originals = (build.RUNTIME / "luafemm-materials-data.lua", ROOT / "docs/materials-catalog.md")
    for generated, original in zip(targets, originals):
        if generated.read_bytes() != original.read_bytes():
            raise SystemExit(f"Generated data are stale: {original}")


def source_files():
    top = ("README.md", "LICENSE", "LICENSES.md", "CHANGELOG.md", "CONTRIBUTING.md",
           "MANIFEST.md", "Makefile", ".editorconfig", ".gitattributes", ".gitignore", ".stylua.toml",
           ".styluaignore")
    files = [ROOT / name for name in top]
    for name in ("tex", "doc", "docs", "examples", "tests", "scripts", "data", ".github"):
        files.extend(p for p in (ROOT / name).rglob("*")
                     if p.is_file() and "__pycache__" not in p.parts and p.name != ".DS_Store")
    return sorted(files)


def archive(path, entries, timestamp):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as output:
        for name, source in sorted(entries):
            info = zipfile.ZipInfo(name, timestamp)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            output.writestr(info, source.read_bytes())
    print(f"Created {path}")


def check_installation(tds):
    # No current-directory or checkout runtime can participate in this lookup.
    with tempfile.TemporaryDirectory(prefix="luafemm-install-") as temporary:
        folder = Path(temporary).resolve()
        texmf = folder / "texmf"
        with zipfile.ZipFile(tds) as package:
            package.extractall(texmf)
        env = os.environ.copy()
        env["LUAFEMM_TEXMF"] = str(texmf)
        env["LUAFEMM_BUILD"] = str(folder / "build")
        env["TEXMFHOME"] = str(texmf)
        for name in ("TEXINPUTS", "LUAINPUTS", "LUA_PATH"):
            env.pop(name, None)
        resolved = subprocess.check_output(["kpsewhich", "luafemm.tex"], cwd=folder, env=env,
                                           text=True).strip()
        if Path(resolved).resolve() != texmf / "tex/generic/luafemm/luafemm.tex":
            raise SystemExit(f"Installed loader not selected: {resolved}")
        subprocess.run(["python3", str(ROOT / "scripts/build.py"), "formats"],
                       cwd=folder, env=env, check=True)
        shutil.copy2(folder / "build/formats/results.json", OUT / "installation.json")
        for name in ("plain", "latex", "context"):
            dest = OUT / "installation-logs" / name
            dest.mkdir(parents=True, exist_ok=True)
            for extension in ("log", "pdf"):
                shutil.copy2(folder / f"build/formats/{name}/{name}.{extension}", dest)
    print("PASS: extracted TDS installation works in all three formats.")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    release, timestamp = version()
    verify_catalogue()
    # Build the current source, not a possibly stale PDF left by an older run.
    build.manual()
    pdf = build.BUILD / "doc/luafemm-manual.pdf"
    files = source_files()
    inventory = set(re.findall(r"^- `([^`]+)`$", (ROOT / "MANIFEST.md").read_text(), re.MULTILINE))
    actual = {p.relative_to(ROOT).as_posix() for p in files}
    if inventory != actual:
        raise SystemExit(f"Update MANIFEST.md: missing {sorted(actual - inventory)}, "
                         f"stale {sorted(inventory - actual)}")
    prefix = f"luafemm-{release}"
    source = OUT / (prefix + ".zip")
    tds = OUT / (prefix + ".tds.zip")
    archive(source, [(f"{prefix}/{p.relative_to(ROOT).as_posix()}", p) for p in files]
            + [(f"{prefix}/doc/luafemm-manual.pdf", pdf)], timestamp)
    entries = []
    for path in files:
        relative = path.relative_to(ROOT)
        if relative.parts[0] == "tex":
            entries.append((relative.as_posix(), path))
        elif relative.parts[0] in ("docs", "examples", "data") or path.name in (
                "README.md", "LICENSE", "LICENSES.md", "CHANGELOG.md", "MANIFEST.md"):
            entries.append(("doc/generic/luafemm/" + relative.as_posix(), path))
    entries.append(("doc/generic/luafemm/luafemm-manual.pdf", pdf))
    archive(tds, entries, timestamp)
    check_installation(tds)
    sums = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}" for p in (source, tds)]
    (OUT / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    print("Local release files verified; no commit, tag or publication performed.")


if __name__ == "__main__":
    main()
