<!--
SPDX-License-Identifier: LPPL-1.3c
Copyright (C) 2026 Christophe Jorssen
Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
LPPL maintenance status: maintained
This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
-->

# Meshing and robustness

The native-path interface defaults to a constrained Delaunay mesher written
entirely in Lua. Triangle is an algorithm reference, not a linked dependency.
The complete key reference and limitations are in the manual.

## Construction

Evaluated TikZ soft paths are flattened into closed polygons. The rectangular
outer boundary and region edges form a planar straight-line graph. The mesher
splits intersections, shared sides and T junctions, merges vertices under the
configured tolerance, inserts interior points at global/local spacing,
constructs a triangulation, recovers constraints and performs legal edge flips.
The last declared region containing a triangle's interior determines its medium.
An explicit air region can therefore cut an opening in a preceding material.

Orientation and incircle predicates use floating filters, then exact expansion
arithmetic when the sign is uncertain. Exactness is with respect to the supplied
binary floating-point coordinates. Intersection construction, tolerance-based
merging and path flattening still use finite precision; predicates do not make
these constructions exact.

Every result is audited for positive orientation, domain area/coverage, recovered
constraints, manifold edge incidence and the local constrained-Delaunay property.
Errors abort rather than returning an incomplete-looking plausible mesh.

## Controls and limits

`femm/problem={mesh size=2}` sets the nominal spacing in millimetres.
`femm/region={material=air,mesh size=.25}` requests denser local points, for
example in a narrow air gap. `mesh tolerance=0` selects an automatic geometric
tolerance; a positive value overrides it. `curve tolerance` controls path
flattening separately. The manual documents their units and inheritance.

This is not Triangle's full quality-refinement algorithm. No minimum angle or
maximum aspect ratio is guaranteed. Narrow gaps, acute corners and local density
transitions may produce thin triangles, and refinement can lower the minimum
angle. Path topology must survive flattening and vertex merging; features below
these tolerances cannot be recovered later. Extremely large or ill-scaled models
can exceed the practical range of the construction or its finite safety limits.

The legacy `grid` mesher accepts axis-aligned region edges only. Use the native
Delaunay path interface for rotated or curved media.

## Evidence

`tests/test_mesh.lua` covers constraint recovery, areas, overlaps, holes, concavity,
oblique geometry and narrow gaps. `tests/test_mesh_stress.lua` adds 49 deterministic
cases including scale/translation changes, near degeneracies and intersections.
`tests/check_predicates.py` compares 3,015 signs against exact rational arithmetic
independent of the Lua expansion implementation. Passing these tests is evidence
for the tested cases, not a proof for every planar input.
