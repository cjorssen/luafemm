<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Numerical formulation

The complete reference is in `doc/luafemm-manual.tex` (build with `make manual`).

The planar potential is **A = Az ez**, with **B = (∂y Az, −∂x Az)**.
The outer rectangular boundary carries Az = 0 in the TikZ interface.
The weak equation is discretized with continuous linear triangular elements;
B is constant in each triangle. Coordinates are converted from model units
to SI before assembly. The default model unit is 0.001 m (one millimetre).

For isotropic media, **H = ν(|B|) B − Hc m**, where m is the imposed
magnetization direction. The B–H table is interpolated linearly and must be
monotone. Beyond its last point, the slope becomes the vacuum slope 1/μ0, keeping the magnetization saturated. The tangent includes the radial derivative
of reluctivity. Newton steps use line search; each symmetric positive system
uses conjugate gradients with an incomplete-Cholesky preconditioner and
controlled diagonal-shift retries when needed.

The positive coercive source contribution on a triangle is proportional to
`Hc * (mx * ∂y Ni − my * ∂x Ni)`. This sign is tested with a piecewise-affine
magnet/air solution and magnetization reversal. Remanence and coercivity are
alternative declarations; a nonlinear permanent magnet uses the shifted FEMM
table with coercivity. See [permanent magnets](permanent-magnets.md).

Field lines are level sets of Az. Triangle intersections are connected through
mesh edges and oriented so their tangents follow B. The rendered paths are
piecewise linear without smoothing, preserving the computed potential levels.
Arrow decorations use PGF's floating `veclen` locally to avoid fixed-point
overflow on very short segments.

Analytic affine and manufactured Poisson solutions test signs and convergence.
A nonlinear residual measures algebraic convergence only: domain truncation,
mesh spacing, curve approximation, material data and interpolation all affect
the physical answer. This release does not claim numerical equivalence with
FEMM/xfemm; a direct cross-solver comparison remains future work.
