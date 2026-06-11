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

# Exact coordinate equality of two points.
_equals2(p, q) = GI.x(p) == GI.x(q) && GI.y(p) == GI.y(q)

function rk_classify_intersection(m::Planar, a0, a1, b0, b1; exact)
    oa0 = rk_orient(m, b0, b1, a0; exact)
    oa1 = rk_orient(m, b0, b1, a1; exact)
    ob0 = rk_orient(m, a0, a1, b0; exact)
    ob1 = rk_orient(m, a0, a1, b1; exact)
    # fully collinear configuration (handles zero-length segments too)
    if oa0 == 0 && oa1 == 0 && ob0 == 0 && ob1 == 0
        a0_on_b = _collinear_between(a0, b0, b1)
        a1_on_b = _collinear_between(a1, b0, b1)
        b0_on_a = _collinear_between(b0, a0, a1)
        b1_on_a = _collinear_between(b1, a0, a1)
        n_inc = a0_on_b + a1_on_b + b0_on_a + b1_on_a
        n_inc == 0 && return SegSegClass(SS_DISJOINT, false, false, false, false)
        # single shared endpoint counts twice (one endpoint of each on the other)
        shared_endpoint_only = n_inc == 2 &&
            ((a0_on_b || a1_on_b) && (b0_on_a || b1_on_a)) &&
            (_equals2(a0, b0) || _equals2(a0, b1) || _equals2(a1, b0) || _equals2(a1, b1))
        kind = shared_endpoint_only ? SS_TOUCH : SS_COLLINEAR
        # zero-length degenerate: a point on a segment is a touch, not an overlap
        if _equals2(a0, a1) || _equals2(b0, b1)
            kind = SS_TOUCH
        end
        return SegSegClass(kind, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    a0_on_b = oa0 == 0 && _collinear_between(a0, b0, b1)
    a1_on_b = oa1 == 0 && _collinear_between(a1, b0, b1)
    b0_on_a = ob0 == 0 && _collinear_between(b0, a0, a1)
    b1_on_a = ob1 == 0 && _collinear_between(b1, a0, a1)
    if a0_on_b || a1_on_b || b0_on_a || b1_on_a
        return SegSegClass(SS_TOUCH, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    if (oa0 > 0) != (oa1 > 0) && oa0 != 0 && oa1 != 0 &&
       (ob0 > 0) != (ob1 > 0) && ob0 != 0 && ob1 != 0
        return SegSegClass(SS_PROPER, false, false, false, false)
    end
    return SegSegClass(SS_DISJOINT, false, false, false, false)
end
