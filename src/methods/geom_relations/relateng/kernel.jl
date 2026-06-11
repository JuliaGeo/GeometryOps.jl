# # RelateKernel API contract

#=
The geometry layer of RelateNG (design doc D1/D2): every coordinate-level
question the topology layer may ask, answered with exact predicates and no
constructed coordinates. Each function takes the manifold as its first
argument; the planar implementations live in `kernel_planar.jl`. A future
`Spherical` kernel implements the same functions and must pass the same
conformance testset (Task 9).

All kernel functions are prefixed `rk_` (RelateKernel) and are internal —
nothing here is exported. Points are coordinate tuples (typically
`Tuple{Float64, Float64}`) obtained via `_tuple_point`; the `exact` flag is a
keyword taking `True()`/`False()` (GeometryOpsCore BoolsAsTypes), threaded
exactly like `Predicates.orient`.

The contract — what every manifold implementation must provide:

    rk_orient(m, a, b, c; exact)::Integer

Orientation of point `c` relative to the oriented segment `(a, b)`:
`> 0` if `c` is to the left, `< 0` if to the right, `0` if collinear
(or `a == b`). With `exact = True()` the sign must be correct even for
adversarial near-collinear inputs.

    rk_point_on_segment(m, p, q0, q1; exact)::Bool

Whether point `p` lies on the closed segment `[q0, q1]`, endpoints included.

    rk_point_in_ring(m, p, ring; exact)::Int8

Location of point `p` relative to the area enclosed by the closed `ring`
(a GeoInterface linestring/linearring, assumed closed regardless of a
repeated last point): one of `LOC_INTERIOR`, `LOC_BOUNDARY`, `LOC_EXTERIOR`.

    rk_interaction_bounds(m, geom)::Extents.Extent

The bounding region within which `geom` can interact with another geometry.
On the plane this is the ordinary extent; other manifolds may need to widen
it (e.g. great-circle edges bulge outside the coordinate box of their
endpoints).

    rk_bounds_disjoint(extA, extB)::Bool
    rk_bounds_covers(extA, extB)::Bool

Conservative interaction-bounds tests used for short-circuiting:
`rk_bounds_disjoint` must only return `true` when no interaction is possible;
`rk_bounds_covers` must only return `true` when `extA` covers `extB` in the
X/Y dimensions. These operate on the extents produced by
`rk_interaction_bounds` and are manifold-generic.
=#

# Manifold-generic helpers

# Whether `p` lies within the coordinate bounding box of segment `(q0, q1)`.
# Valid as an on-segment test only when `p` is already known collinear with
# `(q0, q1)`; shared by manifolds whose segments are coordinate-monotone.
@inline function _collinear_between(p, q0, q1)
    (min(GI.x(q0), GI.x(q1)) <= GI.x(p) <= max(GI.x(q0), GI.x(q1))) &&
    (min(GI.y(q0), GI.y(q1)) <= GI.y(p) <= max(GI.y(q0), GI.y(q1)))
end
