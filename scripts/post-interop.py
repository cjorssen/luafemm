#!/usr/bin/env python3
# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
"""Compare regional energy/coenergy and forces with the optional xfemm oracle.

First run scripts/interop.py --build --xfemm ../xfemm. No native dependency
is introduced into the installed package. This runner rebuilds only the probe.
"""
import json
import re
import subprocess

import interop


def main():
    source, build, root, out = interop.SOURCE, interop.BUILD, interop.ROOT, interop.OUT
    include = [source / 'cfemm' / p for p in
               ('fpproc', 'libfemm', 'libfemm/liblua', 'fsolver', 'fmesher')]
    libs = ['fpproc/libfpproc.a', 'fsolver/libfsolver.a', 'fmesher/libfmesher.a',
            'libfemm/libfemm.a', 'libfemm/libluacomplex.a',
            'fmesher/libtriangle.a', 'fmesher/libtangle.a']
    if not all((build / p).is_file() for p in libs):
        raise SystemExit('Build the optional xfemm oracle with scripts/interop.py --build first')
    probe = out / 'xfemm-probe'
    interop.run(['c++', '-std=c++17', *['-I' + str(p) for p in include],
                 root / 'tests/xfemm_probe.cpp', *[build / p for p in libs],
                 '-o', probe], 'post-probe-build.log')
    interop.run(['texlua', 'tests/post_reference.lua'], 'post-fixtures.log')
    report = []
    for row in (out / 'post.tsv').read_text().splitlines():
        file, x, y, *values = row.split('\t')
        text = subprocess.check_output([str(probe), '--post', file, x, y], text=True, cwd=root)
        match = re.search(r'^POST (.+)$', text, re.MULTILINE)
        if not match:
            raise SystemExit('No native postprocessing results: ' + text)
        lua, native = list(map(float, values)), list(map(float, match.group(1).split()))
        # Both readers use the same mesh and nonlinear table. Compare vector
        # errors against the force scale, including nearly zero components.
        limits = [max(1e-10, abs(native[i]) * 1e-6) for i in range(2)]
        force_scale = max(abs(native[2]), abs(native[3]), 1e-8)
        limits += [force_scale * 1e-7, force_scale * 1e-7, force_scale * .08 * 1e-7]
        passed = all(abs(a-b) <= tol for a, b, tol in zip(lua, native, limits))
        report.append(dict(file=file, lua=lua, xfemm=native, tolerances=limits, passed=passed))
    (out / 'post-results.json').write_text(json.dumps(report, indent=2) + '\n')
    if not all(r['passed'] for r in report):
        raise SystemExit('Postprocessing comparison failed; see build/interop/post-results.json')
    print('PASS: xfemm regional energy, coenergy, force and torque on two nonlinear U meshes')


if __name__ == '__main__':
    main()
