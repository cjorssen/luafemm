# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
"""Independent rational oracle for the exact sign of binary64 determinants."""
from fractions import Fraction as F
from pathlib import Path
import math
import random
import subprocess
import tempfile
import os

rng = random.Random(20260925)
cases = []
def orient(p):
    a, b, c = [[F(x) for x in v] for v in p]
    return (a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0])
def incircle(p):
    a, b, c, d = [[F(x) for x in v] for v in p]
    x, y, z = [[v[i]-d[i] for i in range(2)] for v in (a,b,c)]
    cross = lambda a,b: a[0]*b[1]-a[1]*b[0]
    lift = lambda a: a[0]*a[0]+a[1]*a[1]
    return lift(x)*cross(y,z)+lift(y)*cross(z,x)+lift(z)*cross(x,y)
def case(kind, points):
    val = (orient if kind == 'orient' else incircle)(points)
    cases.append((kind, points, (val>0)-(val<0)))
for _ in range(1500):
    a, b = [[rng.uniform(-1,1) for _ in range(2)] for _ in range(2)]
    t = rng.random()
    c = [a[i]+t*(b[i]-a[i]) for i in range(2)]
    c[0] = math.nextafter(c[0], rng.choice([-math.inf, math.inf]))
    case('orient', [a,b,c])
    angle = rng.random()*2*math.pi
    points = [[math.cos(angle+j*math.pi/2), math.sin(angle+j*math.pi/2)] for j in range(4)]
    points[3][0] = math.nextafter(points[3][0],rng.choice([-math.inf,math.inf]))
    case('incircle', points)
for exponent in (-30,-10,0,10,30):
    s = 2.0**exponent
    case('orient', [[0.,0.],[s,s],[2*s,2*s]])
    case('incircle', [[s,0.],[0.,s],[-s,0.],[0.,-s]])
    case('orient', [[0.,0.],[s,s],[2*s,math.nextafter(2*s,math.inf)]])
with tempfile.TemporaryDirectory(prefix='luafemm-predicates-') as temp:
    script = Path(temp)/'check.lua'
    lines = ["local p=require('luafemm-predicates')"]
    for kind, points, sign in cases:
        coords = ','.join('{'+','.join(float(x).hex() for x in v)+'}' for v in points)
        lines.append(f'do local d=p.{kind}({coords}); print(d>0 and 1 or (d<0 and -1 or 0)) end')
    lines.append("io.stderr:write('exact fallbacks: '..p.stats.orientation..' / '..p.stats.incircle..'\\n')")
    script.write_text('\n'.join(lines))
    run = subprocess.run(['texlua',str(script)],cwd=Path(__file__).resolve().parent.parent,
                         env={**os.environ, "LUA_PATH": "./tex/generic/luafemm/?.lua;;"},
                         capture_output=True,text=True,check=True)
    actual = [int(v) for v in run.stdout.splitlines()]
    assert len(actual)==len(cases)
    for i,(value,(_,points,expected)) in enumerate(zip(actual,cases)):
        assert value==expected, (i,points,value,expected)
    print(f'{len(cases)} exact predicate signs agree with Python Fraction; '+run.stderr.strip())
