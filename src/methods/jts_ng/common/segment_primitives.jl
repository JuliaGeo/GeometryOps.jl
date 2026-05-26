# # Shared NG segment primitives

"""
    NGPrecisionModel

Dispatch hook for coordinate precision policies used by NG primitives.
"""
abstract type NGPrecisionModel end

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
    ng_cross(a, b; exact = True())

Return the robust sign of the 2D cross product `a × b`.
"""
ng_cross(a, b; exact = True()) = Predicates.cross(a, b; exact)

"""
    apply_ng_precision(model, point, [T])

Apply an NG precision policy to a point.  The default policy only converts type.
"""
apply_ng_precision(model, point, ::Type{T} = Float64) where {T} = _tuple_point(point, T)

function _apply_ng_precision_segment(model, segment::Tuple{<:Any,<:Any}, ::Type{T}) where {T}
    return (
        apply_ng_precision(model, segment[1], T),
        apply_ng_precision(model, segment[2], T),
    )
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
    if orientation == line_over
        intr1, intr2 = _ng_order_overlap_intrs(intr1, intr2)
    end
    return NGSegmentIntersection(orientation, intr1, intr2, false, false)
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
