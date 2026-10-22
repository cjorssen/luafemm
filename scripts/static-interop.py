#!/usr/bin/env python3
# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
"""Check static extensions against the optional, previously built xfemm oracle.

ANS checks compare fields on the identical mesh. FEM checks independently mesh
and solve the exported problem. Native tools are developer-only dependencies.
"""
from pathlib import Path
import json
import re
import subprocess

import interop


def main():
    root = interop.ROOT
    out = root / 'build/static'
    probe = interop.OUT / 'xfemm-probe'
    bins = interop.SOURCE / 'cfemm/bin'
    if not probe.is_file():
        raise SystemExit('Build the optional oracle with scripts/interop.py --build first')
    include = [interop.SOURCE / 'cfemm' / p for p in
               ('fpproc', 'libfemm', 'libfemm/liblua', 'fsolver', 'fmesher')]
    libs = ['fpproc/libfpproc.a', 'fsolver/libfsolver.a', 'fmesher/libfmesher.a',
            'libfemm/libfemm.a', 'libfemm/libluacomplex.a',
            'fmesher/libtriangle.a', 'fmesher/libtangle.a']
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'probe-build.log').open('w') as stream:
        subprocess.run(['c++', '-std=c++17', *['-I' + str(p) for p in include],
                        str(root / 'tests/xfemm_probe.cpp'),
                        *[str(interop.BUILD / p) for p in libs], '-o', str(probe)],
                       cwd=root, check=True, stdout=stream, stderr=subprocess.STDOUT)
    subprocess.run(['texlua', 'tests/static_reference.lua'], cwd=root, check=True)
    solved, report = set(), []
    for row in (out / 'reference.tsv').read_text().splitlines():
        name, x, y, *numbers = row.split('\t')
        reference = list(map(float, numbers))
        if name not in solved:
            for program, path in [('fmesher', out / (name + '.fem')), ('fsolver', out / name)]:
                with (out / (name + '-' + program + '.log')).open('w') as stream:
                    subprocess.run([str(bins / program), str(path)], cwd=root,
                                   check=True, stdout=stream, stderr=subprocess.STDOUT)
            solved.add(name)
        for suffix in ['-lua.ans', '.ans']:
            filename = out / (name + suffix)
            text = subprocess.check_output([str(probe), "--current", str(filename), x, y], cwd=root, text=True, stderr=subprocess.STDOUT)
            match = re.search(r'^VALUES (.+)$', text, re.MULTILINE)
            if not match:
                raise SystemExit('No native field sample: ' + text)
            values = list(map(float, match.group(1).split()))
            bscale = max(1e-8, (reference[1] ** 2 + reference[2] ** 2) ** .5)
            hscale = max(1, (reference[3] ** 2 + reference[4] ** 2) ** .5)
            if suffix == '-lua.ans':
                limits = [1e-11, bscale * 1e-7, bscale * 1e-7, hscale * 1e-7, hscale * 1e-7]
            else:
                # Distinct meshes/source quadrature: compare B vector only.
                limits = [float('inf'), bscale * .06, bscale * .06, float('inf'), float('inf')]
            limits.append(max(1e-7, abs(reference[5]) * 1e-8))
            passed = all(abs(a-b) <= tol for a, b, tol in zip(values, reference, limits))
            report.append(dict(file=str(filename.relative_to(root)), x=float(x), y=float(y),
                               luafemm=reference, xfemm=values, passed=passed))
    (out / 'xfemm-results.json').write_text(json.dumps(report, indent=2) + '\n')
    if not all(r['passed'] for r in report):
        raise SystemExit('Static interoperability mismatch: see build/static/xfemm-results.json')
    print('PASS: four independent static FEM solves and identical-mesh ANS A/B/H/J comparisons')


if __name__ == '__main__':
    main()
