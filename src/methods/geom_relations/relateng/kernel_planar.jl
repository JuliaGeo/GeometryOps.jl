# # Planar RelateKernel

#=
Planar implementation of the RelateKernel contract declared in `kernel.jl`.
Orientation goes through `Predicates.orient` (AdaptivePredicates when
`exact = True()`), point-in-ring reuses the existing Hao–Sun ray-crossing
machinery (`_point_filled_curve_orientation`), and bounds are plain extents.
No coordinates are ever constructed.
=#

rk_orient(::Planar, a, b, c; exact) = Predicates.orient(a, b, c; exact)

function rk_point_on_segment(m::Planar, p, q0, q1; exact)
    rk_orient(m, q0, q1, p; exact) == 0 || return false
    return _collinear_between(p, q0, q1)
end

function rk_point_in_ring(m::Planar, p, ring; exact)
    o = _point_filled_curve_orientation(m, p, ring; in = point_in, on = point_on, out = point_out, exact)
    o == point_in && return LOC_INTERIOR
    o == point_on && return LOC_BOUNDARY
    return LOC_EXTERIOR
end

rk_interaction_bounds(::Planar, geom) = GI.extent(geom, fallback = true)

rk_bounds_disjoint(extA, extB) = !Extents.intersects(extA, extB)

function rk_bounds_covers(extA, extB)
    (extA.X[1] <= extB.X[1] && extB.X[2] <= extA.X[2]) &&
    (extA.Y[1] <= extB.Y[1] && extB.Y[2] <= extA.Y[2])
end
