# # Shared NG segment primitives

"""
    NGPrecisionModel

Dispatch hook for coordinate precision policies used by NG primitives.
"""
abstract type NGPrecisionModel end

"""
    FixedPrecisionModel(scale; offset = (0, 0))

Round coordinates to a JTS-style fixed precision grid before NG primitive decisions.
Negative scales are treated as exact grid sizes, matching JTS `PrecisionModel`.
"""
struct FixedPrecisionModel{T} <: NGPrecisionModel
    scale::T
    offset::Tuple{T,T}
end

function FixedPrecisionModel(scale::Real; offset = (0, 0))
    iszero(scale) && throw(ArgumentError("FixedPrecisionModel scale must be non-zero."))
    T = promote_type(typeof(scale), Float64)
    jts_scale = scale < 0 ? inv(abs(T(scale))) : abs(T(scale))
    return FixedPrecisionModel(jts_scale, (T(offset[1]), T(offset[2])))
end

"""
    NGSegmentIntersection

Stable wrapper for segment-intersection kind, points, and segment fractions.
"""
struct NGSegmentIntersection{T}
    orientation::LineOrientation
    point1::Tuple{T,T}
    fraction1::Tuple{T,T}
    point2::Tuple{T,T}
    fraction2::Tuple{T,T}
    is_degenerate_a::Bool
    is_degenerate_b::Bool
end

function NGSegmentIntersection(
    orientation::LineOrientation,
    intr1::Tuple,
    intr2::Tuple,
    is_degenerate_a::Bool,
    is_degenerate_b::Bool,
)
    point1, fraction1 = intr1
    point2, fraction2 = intr2
    T = promote_type(eltype(point1), eltype(point2), eltype(fraction1), eltype(fraction2))
    return NGSegmentIntersection{T}(
        orientation,
        point1,
        fraction1,
        point2,
        fraction2,
        is_degenerate_a,
        is_degenerate_b,
    )
end

function _ng_no_segment_intersection(::Type{T}, is_degenerate_a::Bool, is_degenerate_b::Bool) where {T}
    intr = ((zero(T), zero(T)), (zero(T), zero(T)))
    return NGSegmentIntersection(line_out, intr, intr, is_degenerate_a, is_degenerate_b)
end

"""
    ng_orientation(a, b, c; exact = True())

Return the robust orientation sign of point `c` relative to directed segment `a => b`.
"""
ng_orientation(a, b, c; exact = True()) = Predicates.orient(a, b, c; exact)

"""
    ng_jts_orientation_index(p1, p2, q)

Return JTS `Orientation.index` semantics for point `q` relative to `p1 => p2`.
"""
function ng_jts_orientation_index(p1, p2, q)
    p1x, p1y = _tuple_point(p1, Float64)
    p2x, p2y = _tuple_point(p2, Float64)
    qx, qy = _tuple_point(q, Float64)

    index = _ng_jts_orientation_index_filter(p1x, p1y, p2x, p2y, qx, qy)
    index <= 1 && return index

    return _ng_jts_orientation_index_dd(p1x, p1y, p2x, p2y, qx, qy)
end

function _ng_jts_orientation_index_filter(p1x, p1y, p2x, p2y, qx, qy)
    detleft = (p1x - qx) * (p2y - qy)
    detright = (p1y - qy) * (p2x - qx)
    det = detleft - detright

    if detleft > 0.0
        detright <= 0.0 && return _ng_jts_signum(det)
        detsum = detleft + detright
    elseif detleft < 0.0
        detright >= 0.0 && return _ng_jts_signum(det)
        detsum = -detleft - detright
    else
        return _ng_jts_signum(det)
    end

    errbound = 1e-15 * detsum
    (det >= errbound || -det >= errbound) && return _ng_jts_signum(det)
    return 2
end

function _ng_jts_orientation_index_dd(p1x, p1y, p2x, p2y, qx, qy)
    dx1 = _ng_jts_dd_add(_ng_jts_dd(p2x), -p1x)
    dy1 = _ng_jts_dd_add(_ng_jts_dd(p2y), -p1y)
    dx2 = _ng_jts_dd_add(_ng_jts_dd(qx), -p2x)
    dy2 = _ng_jts_dd_add(_ng_jts_dd(qy), -p2y)
    det = _ng_jts_dd_subtract(
        _ng_jts_dd_multiply(dx1, dy2),
        _ng_jts_dd_multiply(dy1, dx2),
    )
    return _ng_jts_dd_signum(det)
end

_ng_jts_signum(value) = value > zero(value) ? 1 : value < zero(value) ? -1 : 0

const _NG_JTS_DD_SPLIT = 134217729.0

_ng_jts_dd(value) = (Float64(value), 0.0)

function _ng_jts_dd_add((hi, lo)::Tuple{Float64,Float64}, y::Float64)
    S = hi + y
    e = S - hi
    s = S - e
    s = (y - e) + (hi - s)
    f = s + lo
    H = S + f
    h = f + (S - H)
    zhi = H + h
    zlo = h + (H - zhi)
    return zhi, zlo
end

function _ng_jts_dd_subtract(
    (hi, lo)::Tuple{Float64,Float64},
    (yhi0, ylo0)::Tuple{Float64,Float64},
)
    yhi, ylo = -yhi0, -ylo0
    S = hi + yhi
    T = lo + ylo
    e = S - hi
    f = T - lo
    s = S - e
    t = T - f
    s = (yhi - e) + (hi - s)
    t = (ylo - f) + (lo - t)
    e = s + T
    H = S + e
    h = e + (S - H)
    e = t + h
    zhi = H + e
    zlo = e + (H - zhi)
    return zhi, zlo
end

function _ng_jts_dd_multiply(
    (hi, lo)::Tuple{Float64,Float64},
    (yhi, ylo)::Tuple{Float64,Float64},
)
    C = _NG_JTS_DD_SPLIT * hi
    hx = C - hi
    c = _NG_JTS_DD_SPLIT * yhi
    hx = C - hx
    tx = hi - hx
    hy = c - yhi
    C = hi * yhi
    hy = c - hy
    ty = yhi - hy
    c = ((((hx * hy - C) + hx * ty) + tx * hy) + tx * ty) + (hi * ylo + lo * yhi)
    zhi = C + c
    hx = C - zhi
    zlo = c + hx
    return zhi, zlo
end

function _ng_jts_dd_signum((hi, lo)::Tuple{Float64,Float64})
    hi > 0 && return 1
    hi < 0 && return -1
    lo > 0 && return 1
    lo < 0 && return -1
    return 0
end

"""
    ng_jts_point_on_segment(point, p0, p1)

Return whether `point` lies on a segment using JTS `PointLocation.isOnSegment`.
"""
function ng_jts_point_on_segment(point, p0, p1)
    point = _tuple_point(point, Float64)
    p0 = _tuple_point(p0, Float64)
    p1 = _tuple_point(p1, Float64)

    _point_in_extent(point, ng_segment_extent((p0, p1), Float64)) || return false
    point == p0 && return true
    return ng_jts_orientation_index(p0, p1, point) == 0
end

"""
    ng_jts_locate_point_in_ring(point, ring)

Locate `point` against a ring using JTS `RayCrossingCounter` semantics.
"""
function ng_jts_locate_point_in_ring(point, ring)
    npoints = GI.npoint(ring)
    npoints < 2 && return loc_exterior

    point = _tuple_point(point, Float64)
    crossing_count = 0
    for i in 2:npoints
        p1 = _tuple_point(GI.getpoint(ring, i), Float64)
        p2 = _tuple_point(GI.getpoint(ring, i - 1), Float64)
        count, is_on_segment = _ng_jts_count_ray_segment(point, p1, p2)
        is_on_segment && return loc_boundary
        crossing_count += count
    end
    return isodd(crossing_count) ? loc_interior : loc_exterior
end

function _ng_jts_count_ray_segment(point, p1, p2)
    px, py = point
    p1x, p1y = p1
    p2x, p2y = p2

    p1x < px && p2x < px && return 0, false
    (px == p2x && py == p2y) && return 0, true

    if p1y == py && p2y == py
        minx, maxx = minmax(p1x, p2x)
        return 0, px >= minx && px <= maxx
    end

    if ((p1y > py && p2y <= py) || (p2y > py && p1y <= py))
        orient = ng_jts_orientation_index(p1, p2, point)
        orient == 0 && return 0, true
        p2y < p1y && (orient = -orient)
        return orient == 1 ? 1 : 0, false
    end
    return 0, false
end

"""
    ng_cross(a, b; exact = True())

Return the robust sign of the 2D cross product `a × b`.
"""
ng_cross(a, b; exact = True()) = Predicates.cross(a, b; exact)

"""
    apply_ng_precision(model, point, [T])

Apply an NG precision policy to a point.  The default policy only converts type.
"""
apply_ng_precision(model, point, ::Type{T} = Float64) where {T} = _tuple_point(point, T)

function apply_ng_precision(model::FixedPrecisionModel, point, ::Type{T} = Float64) where {T}
    x = _jts_precision_round((T(GI.x(point)) - T(model.offset[1])) * T(model.scale)) / T(model.scale) +
        T(model.offset[1])
    y = _jts_precision_round((T(GI.y(point)) - T(model.offset[2])) * T(model.scale)) / T(model.scale) +
        T(model.offset[2])
    return (x, y)
end

_jts_precision_round(value) = floor(value + one(value) / 2)

function _apply_ng_precision_segment(model, segment::Tuple{<:Any,<:Any}, ::Type{T}) where {T}
    return (
        apply_ng_precision(model, segment[1], T),
        apply_ng_precision(model, segment[2], T),
    )
end

function _apply_ng_precision_intr(model, intr::Tuple, ::Type{T}) where {T}
    point, fraction = intr
    return apply_ng_precision(model, point, T), fraction
end

"""
    ng_is_degenerate_segment(segment)

Return `true` when a segment has identical endpoints.
"""
ng_is_degenerate_segment((p1, p2)::Tuple{<:Any,<:Any}) = _tuple_point(p1) == _tuple_point(p2)

"""
    ng_segment_extent(segment)

Return the 2D extent of a segment.
"""
function ng_segment_extent((p1, p2)::Tuple{<:Any,<:Any}, ::Type{T} = Float64) where {T}
    p1t, p2t = _tuple_point(p1, T), _tuple_point(p2, T)
    return Extents.Extent(X = minmax(p1t[1], p2t[1]), Y = minmax(p1t[2], p2t[2]))
end

"""
    ng_segments_maybe_intersect(segment_a, segment_b)

Cheap extent test for whether two segments may intersect.
"""
function ng_segments_maybe_intersect(segment_a, segment_b, ::Type{T} = Float64) where {T}
    return Extents.intersects(ng_segment_extent(segment_a, T), ng_segment_extent(segment_b, T))
end

function _ng_order_overlap_intrs(intr1::Tuple, intr2::Tuple)
    _, (α1, β1) = intr1
    _, (α2, β2) = intr2
    (α1, β1) <= (α2, β2) && return intr1, intr2
    return intr2, intr1
end

"""
    ng_segment_intersection(segment_a, segment_b, [T]; kwargs...)

Intersect two non-degenerate segments using GeometryOps' robust segment primitive.
"""
function ng_segment_intersection(
    segment_a::Tuple{<:Any,<:Any},
    segment_b::Tuple{<:Any,<:Any},
    ::Type{T} = Float64;
    manifold::Manifold = Planar(),
    exact = True(),
    precision_model = nothing,
) where {T}
    segment_a = _apply_ng_precision_segment(precision_model, segment_a, T)
    segment_b = _apply_ng_precision_segment(precision_model, segment_b, T)
    degenerate_a = ng_is_degenerate_segment(segment_a)
    degenerate_b = ng_is_degenerate_segment(segment_b)
    (degenerate_a || degenerate_b) &&
        return _ng_no_segment_intersection(T, degenerate_a, degenerate_b)

    orientation, intr1, intr2 = _intersection_point(manifold, T, segment_a, segment_b; exact)
    if orientation == line_cross
        intr1 = _ng_stable_cross_intr(segment_a, segment_b, T)
    end
    if orientation == line_over
        intr1, intr2 = _ng_order_overlap_intrs(intr1, intr2)
    end
    if orientation != line_out
        intr1 = _apply_ng_precision_intr(precision_model, intr1, T)
        intr2 = _apply_ng_precision_intr(precision_model, intr2, T)
    end
    return NGSegmentIntersection(orientation, intr1, intr2, false, false)
end

"""
    _ng_stable_cross_intr(segment_a, segment_b, T)

Compute a crossing point in high precision so equivalent NG nodes share a key.
"""
function _ng_stable_cross_intr(segment_a, segment_b, ::Type{T}) where {T}
    return setprecision(BigFloat, 256) do
        (a1x, a1y), (a2x, a2y) = _ng_big_segment(segment_a)
        (b1x, b1y), (b2x, b2y) = _ng_big_segment(segment_b)
        adx, ady = a2x - a1x, a2y - a1y
        bdx, bdy = b2x - b1x, b2y - b1y
        bax, bay = b1x - a1x, b1y - a1y
        denom = adx * bdy - ady * bdx
        alpha = (bax * bdy - bay * bdx) / denom
        beta = (bax * ady - bay * adx) / denom
        x = abs(adx) >= abs(bdx) ? a1x + alpha * adx : b1x + beta * bdx
        y = abs(ady) >= abs(bdy) ? a1y + alpha * ady : b1y + beta * bdy
        return (T(x), T(y)), (
            clamp(T(alpha), eps(T), one(T) - eps(T)),
            clamp(T(beta), eps(T), one(T) - eps(T)),
        )
    end
end

function _ng_big_segment((p1, p2))
    p1 = _tuple_point(p1)
    p2 = _tuple_point(p2)
    return (BigFloat(p1[1]), BigFloat(p1[2])), (BigFloat(p2[1]), BigFloat(p2[2]))
end

"""
    ng_has_intersection(intersection)

Return `true` when a segment-intersection result has at least one point.
"""
ng_has_intersection(intersection::NGSegmentIntersection) = intersection.orientation != line_out

"""
    ng_intersection_points(intersection)

Return the concrete point or points represented by an NG segment intersection.
"""
function ng_intersection_points(intersection::NGSegmentIntersection)
    intersection.orientation == line_out && return typeof(intersection.point1)[]
    intersection.orientation == line_over && return [intersection.point1, intersection.point2]
    return [intersection.point1]
end

"""
    ng_is_closed_segment_string(segment_string)

Return `true` when a segment string's first and last coordinates are equal.
"""
function ng_is_closed_segment_string(segment_string::NGSegmentString)
    length(segment_string.points) < 2 && return false
    return first(segment_string.points) == last(segment_string.points)
end

ng_is_closed_segment_string(points::AbstractVector{<:Tuple}) =
    length(points) >= 2 && first(points) == last(points)

"""
    ng_is_containing_segment(segment_string, segment_index, point)

RelateNG half-closed ownership test for assigning vertex hits to one segment.
"""
function ng_is_containing_segment(segment_string::NGSegmentString, segment_index::Integer, point)
    return ng_is_containing_segment(segment_string.points, segment_index, point)
end

function ng_is_containing_segment(points::AbstractVector{<:Tuple}, segment_index::Integer, point)
    1 <= segment_index < length(points) ||
        throw(BoundsError(points, segment_index))

    point = _tuple_point(point, eltype(first(points)))
    segment_start = points[segment_index]
    segment_end = points[segment_index + 1]

    point == segment_start && return true
    if point == segment_end
        is_final_segment = segment_index == length(points) - 1
        return !ng_is_closed_segment_string(points) && is_final_segment
    end
    return true
end
