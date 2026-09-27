#!/usr/bin/env python3
# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
"""Optional independent xfemm validation; never used by installed TeX documents.

Use --build to copy a local xfemm checkout and build the native oracle in build/.
CMake may fetch xfemm's Tangle dependency. The upstream checkout is never edited.
Disk ANS loading is exercised through a small C++ FPProc client, since xfemm's
Lua open() command in the audited revision does not accept answer files.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'build/interop'
SOURCE = OUT / 'xfemm-source'
BUILD = OUT / 'xfemm-build'


def run(command, log=None):
    if log:
        with (OUT / log).open('w') as stream:
            subprocess.run(list(map(str, command)), cwd=ROOT, check=True,
                           stdout=stream, stderr=subprocess.STDOUT)
    else:
        return subprocess.check_output(list(map(str, command)), cwd=ROOT, text=True)


def build(checkout):
    revision = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    stamp = OUT / 'source-revision.txt'
    if stamp.exists() and stamp.read_text().strip() != revision:
        raise SystemExit('The copied xfemm revision differs; use a fresh build/interop directory.')
    if not SOURCE.exists():
        shutil.copytree(checkout, SOURCE, ignore=shutil.ignore_patterns('.git', 'build', '__pycache__'))
    stamp.write_text(revision + '\n')
    flags = []
    if platform.system() == 'Darwin':
        # Compatibility adaptations to the disposable upstream copy only.
        for p in (SOURCE / 'cfemm').rglob('*'):
            if p.suffix in ('.h', '.cpp', '.c'):
                text = p.read_text(errors='surrogateescape')
                fixed = text.replace('#include <malloc.h>', '#include <stdlib.h>')
                fixed = fixed.replace('#include "malloc.h"', '#include <stdlib.h>')
                fixed = fixed.replace('std::not1(std::ptr_fun<int, int>(std::isspace))',
                                      '[](unsigned char c) { return !std::isspace(c); }')
                if fixed != text:
                    p.write_text(fixed, errors='surrogateescape')
        flags = ['-DEXTRA_CMAKE_CXX_FLAGS=-Dsincos=__sincos']
    run(['cmake', '-S', SOURCE / 'cfemm', '-B', BUILD, '-DCMAKE_BUILD_TYPE=Release', *flags],
        'xfemm-configure.log')
    run(['cmake', '--build', BUILD, '-j', '4', '--target', 'fmesher', 'fsolver', 'fpproc'],
        'xfemm-build.log')
    include = [SOURCE / 'cfemm' / p for p in ('fpproc', 'libfemm', 'libfemm/liblua', 'fsolver', 'fmesher')]
    libs = ['fpproc/libfpproc.a', 'fsolver/libfsolver.a', 'fmesher/libfmesher.a',
            'libfemm/libfemm.a', 'libfemm/libluacomplex.a', 'fmesher/libtriangle.a', 'fmesher/libtangle.a']
    run(['c++', '-std=c++17', *['-I' + str(p) for p in include], ROOT / 'tests/xfemm_probe.cpp',
         *[BUILD / p for p in libs], '-o', OUT / 'xfemm-probe'], 'probe-build.log')
    revision = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
    (OUT / 'native-build.json').write_text(json.dumps(dict(xfemm_revision=revision,
        platform=platform.platform(), compiler=run(['c++', '--version']).splitlines()[0],
        macos_adaptations=bool(flags)), indent=2) + '\n')


def verify():
    probe = OUT / 'xfemm-probe'
    if not probe.is_file():
        raise SystemExit('Native oracle missing: rerun with --build --xfemm ../xfemm')
    run(['texlua', 'tests/test_interchange.lua'], 'file-tests.log')
    run(['texlua', 'tests/interop_cases.lua'], 'cases.log')
    bins = SOURCE / 'cfemm/bin'
    results = []
    for row in (OUT / 'samples.tsv').read_text().splitlines():
        parts = row.split('\t')
        path = ROOT / parts[0]
        reference = list(map(float, parts[3:]))
        snapshot = path.name.endswith('-lua.ans')
        if not snapshot:
            run([bins / 'fmesher', path.with_suffix('.fem')], path.stem + '-mesh.log')
            run([bins / 'fsolver', path.with_suffix('')], path.stem + '-solve.log')
        text = run([probe, path, parts[1], parts[2]])
        actual = list(map(float, text.split('VALUES ')[1].split()))[:len(reference)]
        error = max(abs(a-b) / max(1, abs(b)) for a, b in zip(actual, reference))
        # Same mesh/potential: serialization and constitutive agreement.
        # Independent solve: a 2% relative B_y gate on the finest draft U mesh.
        if snapshot and error > 1e-7:
            raise SystemExit(f'ANS mismatch: {path}: {error}')
        if not snapshot and abs(actual[2] / reference[2] - 1) > .02:
            raise SystemExit(f'Independent solve mismatch: {path}: {actual} != {reference}')
        results.append(dict(file=parts[0],x=float(parts[1]),y=float(parts[2]),
                            reference=reference,xfemm=actual,scaled_error=error,
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    for name,x,y in [('uniform-lua',8,9),('annulus-lua',5,0),('annulus-warm',5,0)]:
        actual=list(map(float,run([probe,OUT/(name+'.ans'),x,y]).split('VALUES ')[1].split()))
        if abs(actual[2]+.25)>1e-8 or abs(actual[1])>1e-8:
            raise SystemExit(f'Affine patch mismatch: {name}')
    # Independent solution of constant real type-2 mixed conditions.
    run([bins / 'fmesher', OUT / 'mixed.fem'], 'mixed-mesh.log')
    run([bins / 'fsolver', OUT / 'mixed'], 'mixed-solve.log')
    for filename in ('mixed.ans', 'mixed-lua.ans'):
        values=list(map(float,run([probe,OUT/filename,8,9]).split('VALUES ')[1].split()))
        if abs(values[0]-.002)>1e-8 or abs(values[1])+abs(values[2])>1e-7:
            raise SystemExit(f'Mixed-boundary mismatch: {filename}')
    # Slotted machine snapshots and independently meshed native solves.
    run(['texlua', 'tests/interop_machines.lua'], 'machine-cases.log')
    solved = set()
    machine_results = []
    machine_rows = (OUT / 'machine-samples.tsv').read_text().splitlines()
    field_scales = {}
    for row in machine_rows:
        name, x, y, az, bx, by, *rest = row.split('\t')
        field_scales[name] = max(field_scales.get(name, 0), (float(bx)**2 + float(by)**2)**.5)
    for row in machine_rows:
        name, x, y, *numbers = row.split('\t')
        reference = list(map(float, numbers))
        if name not in solved:
            run([bins / 'fmesher', OUT / (name + '.fem')], name + '-mesh.log')
            run([bins / 'fsolver', OUT / name], name + '-solve.log')
            solved.add(name)
        for suffix in ('-lua.ans', '.ans'):
            actual = list(map(float, run([probe, OUT / (name + suffix), x, y]).split('VALUES ')[1].split()))
            if suffix == '-lua.ans':
                error = max(abs(a-b)/max(1,abs(b)) for a,b in zip(actual,reference))
                if error > 1e-7:
                    raise SystemExit(f'Machine snapshot mismatch: {name}: {error}')
            else:
                # Normalize by the case's peak sampled field: a pointwise
                # relative error is undefined at angular field zeros.
                scale = field_scales[name]
                error = ((actual[1]-reference[1])**2 + (actual[2]-reference[2])**2)**.5 / scale
                if error > .03:
                    raise SystemExit(f'Independent machine B mismatch: {name}: {error}')
            machine_results.append(dict(file=name+suffix,x=float(x),y=float(y),
                                        reference=reference,xfemm=actual,error=error))
    (OUT / 'machine-results.json').write_text(json.dumps(machine_results,indent=2)+'\n')
    print(f'PASS: {len(machine_results)} machine snapshot/native comparisons.')
    (OUT / 'results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(f'PASS: {len(results)} independent A/B/H samples plus three affine patches; see build/interop/results.json')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--xfemm',type=Path,default=ROOT.parent/'xfemm')
    parser.add_argument('--build',action='store_true')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if args.build:
        build(args.xfemm.resolve())
    verify()
