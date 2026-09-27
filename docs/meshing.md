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

Evaluated TikZ soft paths are flattened into closed polygonal contours.
Compound regions use native PGF winding or parity fill semantics. Every contour
is constrained; material cavities retain earlier media or air and remain meshed.
Refinement is requested inside the filled area and along all contour edges,
not throughout the interior of a hole. The rectangular
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

## Optional quality refinement

`femm/problem={minimum angle=20,maximum area=0,max Steiner points=5000}`
requests an angle target in degrees, with no area limit in this example.
Both targets default to zero (disabled); a positive maximum area uses model
units squared and applies to every triangle, including background air.
The maximum accepted angle target is 30 degrees. These keys are ordinary
problem options and compose in PGF styles before the model is meshed.

The Lua refinement module bisects encroached constrained segments first, then
processes deficient triangles in shortest-edge order. A filtered, exact-fallback
dot product classifies points against a segment's open diametral disk. Proposed
circumcentres that encroach a constraint or cross one during point location
request a segment split instead. Accepted vertices trigger local Lawson flips;
new triangles and adjacent segments return to their queues. Exterior tags are
inherited by both children of a split segment, preserving boundary assembly.
The final audit independently checks the quality targets as well as topology.

`steiner_points`, `segment_splits`, `circumcenters` and `exact_diametral` expose
the work performed. The insertion budget counts both splits and circumcentres,
but excludes the initial segment subdivision and lattice points. Budget or
precision exhaustion fails explicitly; it never relaxes the request or saves
a partial mesh as a success. The total 40,000-vertex cap still applies.
Quality options and statistics are included in cache revision 4.

## Controls and limits

`femm/problem={mesh size=2}` sets the nominal spacing in millimetres.
`femm/region={material=air,mesh size=.25}` requests denser local points, for
example in a narrow air gap. `mesh tolerance=0` selects an automatic geometric
tolerance; a positive value overrides it. `curve tolerance` controls path
flattening separately. The manual documents their units and inheritance.

This is not Triangle's full quality-refinement algorithm. There is no specialised
protection for acute input corners, coarsening or smooth-grading guarantee.
A wedge below the requested angle cannot satisfy that angle everywhere; very
narrow gaps can exceed the budget. Candidate disk checks are conservative across
all constraints and may split more than a visibility-aware algorithm. Successful
quality refinement meets the requested angle/area criteria within audit margins
(1e-8 degrees and 1e-10 relative area); successful refinement for arbitrary inputs
is not guaranteed. With the targets disabled, spacing refinement alone may lower the
minimum angle. Path topology must survive flattening and vertex merging; features below
these tolerances cannot be recovered later. Extremely large or ill-scaled models
can exceed the practical range of the construction or its finite safety limits.

The legacy `grid` mesher accepts axis-aligned region edges only. Use the native
Delaunay path interface for rotated or curved media.

## Evidence

`tests/test_mesh.lua` covers constraint recovery, areas, overlaps, holes, concavity,
oblique geometry and narrow gaps. `tests/test_mesh_stress.lua` adds 49 deterministic
cases including scale/translation changes, near degeneracies and intersections.
`tests/test_quality.lua` independently recomputes angles and areas after vertex
renumbering, integrates materials, exercises compound cavities and a thin gap,
checks scale/translation changes and an affine Robin patch, and verifies cache
invalidation and explicit budget failure. `examples/mesh-quality.tex` compares
the same U device with its quality target disabled/enabled: approximately
7.29 to 20.00 degrees with 148 additional points in this revision.
`tests/check_predicates.py` compares 4,530 signs against exact rational arithmetic
independent of the Lua expansion implementation. Passing these tests is evidence
for the tested cases, not a proof for every planar input.
