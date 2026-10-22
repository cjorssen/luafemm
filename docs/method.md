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
The outer rectangular boundary carries Az = 0 by default. `\femmboundary`
sets Dirichlet, Neumann or nonnegative Robin data on each complete side.
With the counterclockwise tangent, Neumann specifies Ht = g and Robin
specifies Ht = c Az + g. The weak residual adds the integral of `(c Az + g) Ni`
on these edges, including the coercive contribution in the definition of H.
Affine boundary data use SI coordinates and consistent P1 edge integration.
An all-Neumann problem requires net current equal to the boundary integral
of Ht; only after checking this does the solver fix one gauge node.
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
the physical answer. Selected planar DC cases have been compared with independent
xfemm solves and postprocessing; see [the validation record](interoperability.md).
Those checks do not establish universal numerical equivalence with FEMM/xfemm.


## Static tensor, boundary and circuit extensions

The linear constitutive law uses a rotated positive symmetric reluctivity
tensor. Both residual and tangent use the same tensor; energy is
`B^T nu B/2 - Hc m.B`. Isotropic nonlinear laws retain their radial tangent.
Spatial coercive directions are integrated using the symmetric three-point
triangle rule, consistently with reported H and magnet energy.

Native domains and tagged paths compile through the same planar arrangement
as FEM import/export. Paired exterior traces receive the union of their
transformed boundary vertices. A signed union-find supplies the exact linear
constraint `A_slave = sign*A_master`; assembly is reduced as `T^T K T` and
`T^T R`. The nullspace check distinguishes a genuine additive constant from
an antiperiodic constraint that already removes it. No penalty stiffness is
introduced. Drawing repetitions use the same rigid transform and reverse
antiperiodic arrows.

DC circuit currents are resolved from active leg areas and explicit conductor
resistances before magnetic assembly. Parallel branches split the current by
conductance. Closed-winding linkage integrates the P1 potential exactly and
uses signed turns; the tests check its independent relation to magnetic energy.
These real static circuits do not model induction or dynamic circuit unknowns.
