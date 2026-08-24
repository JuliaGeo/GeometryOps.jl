# # Geometry Intersection
export intersection, intersection_points

"""
    Enum LineOrientation
Enum for the orientation of a line with respect to a curve. A line can be
`line_cross` (crossing over the curve), `line_hinge` (crossing the endpoint of the curve),
`line_over` (collinear with the curve), or `line_out` (not interacting with the curve).
"""
@enum LineOrientation line_cross=1 line_hinge=2 line_over=3 line_out=4

"""
    intersection(geom_a, geom_b, [T::Type]; target::Type, fix_multipoly = UnionIntersectingPolygons())

Return the intersection between two geometries as a list of geometries. Return an empty list
if none are found. The type of the list will be constrained as much as possible given the
input geometries. Furthermore, the user can provide a `target` type as a keyword argument and
a list of target geometries found in the intersection will be returned. The user can also
provide a float type that they would like the points of returned geometries to be. If the
user is taking a intersection involving one or more multipolygons, and the multipolygon
might be comprised of polygons that intersect, if `fix_multipoly` is set to an
`IntersectingPolygons` correction (the default is `UnionIntersectingPolygons()`), then the
needed multipolygons will be fixed to be valid before performing the intersection to ensure
a correct answer. Only set `fix_multipoly` to nothing if you know that the multipolygons are
valid, as it will avoid unneeded computation. 

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
inter_points = GO.intersection(line1, line2; target = GI.PointTrait())
GI.coordinates.(inter_points)

# output
1-element Vector{Vector{Float64}}:
 [125.58375366067548, -14.83572303404496]
```
"""
function intersection(
    alg::FosterHormannClipping, geom_a, geom_b, ::Type{T}=Float64; target=nothing, kwargs...
) where {T<:AbstractFloat}
    return _intersection(
        alg, TraitTarget(target), T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b;
        exact = True(), kwargs...,
    )
end
# fallback definitions
# if no manifold - assume planar (until we have best_manifold)
function intersection(
    geom_a, geom_b, ::Type{T}=Float64; target=nothing, kwargs...
) where {T<:AbstractFloat}
    return intersection(FosterHormannClipping(Planar()), geom_a, geom_b; target, kwargs...)
end
# if manifold but no algorithm - assume FosterHormannClipping with provided manifold.
function intersection(m::Manifold, geom_a, geom_b, ::Type{T}=Float64; target=nothing, kwargs...) where {T<:AbstractFloat}
    return intersection(FosterHormannClipping(m), geom_a, geom_b; target, kwargs...)
end

# Curve-Curve Intersections with target Point
_intersection(
    alg::FosterHormannClipping, ::TraitTarget{GI.PointTrait}, ::Type{T},
    trait_a::Union{GI.LineTrait, GI.LineStringTrait, GI.LinearRingTrait}, geom_a,
    trait_b::Union{GI.LineTrait, GI.LineStringTrait, GI.LinearRingTrait}, geom_b;
    kwargs...,
) where T = _intersection_points(alg.manifold, alg.accelerator, T, trait_a, geom_a, trait_b, geom_b)

#= Polygon-Polygon Intersections with target Polygon
The algorithm to determine the intersection was adapted from "Efficient clipping
of efficient polygons," by Greiner and Hormann (1998).
DOI: https://doi.org/10.1145/274363.274364 =#
function _intersection(
    alg::FosterHormannClipping, ::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b;
    exact, kwargs...,
) where {T}
    # First we get the exteriors of 'poly_a' and 'poly_b'
    ext_a = GI.getexterior(poly_a)
    ext_b = GI.getexterior(poly_b)
    # Then we find the intersection of the exteriors
    a_list, b_list, a_idx_list = _build_ab_list(alg, T, ext_a, ext_b, _inter_delay_cross_f, _inter_delay_bounce_f; exact)
    polys = _trace_polynodes(alg, T, a_list, b_list, a_idx_list, _inter_step, poly_a, poly_b)
    if isempty(polys) # no crossing points, determine if either poly is inside the other
        a_in_b, b_in_a = _find_non_cross_orientation(alg, a_list, b_list, ext_a, ext_b; exact)
        if a_in_b
            push!(polys, GI.Polygon([tuples(ext_a)]))
        elseif b_in_a
            push!(polys, GI.Polygon([tuples(ext_b)]))
        end
    end
    remove_idx = falses(length(polys))
    # If the original polygons had holes, take that into account.
    if GI.nhole(poly_a) != 0 || GI.nhole(poly_b) != 0
        hole_iterator = Iterators.flatten((GI.gethole(poly_a), GI.gethole(poly_b)))
        _add_holes_to_polys!(alg, T, polys, hole_iterator, remove_idx; exact)
    end
    # Remove unneeded collinear points on same edge
    _remove_collinear_points!(alg, polys, remove_idx, poly_a, poly_b)
    return polys
end
# # Helper functions for Intersections with Greiner and Hormann Polygon Clipping

#= When marking the crossing status of a delayed crossing, the chain start point is bouncing
when the start point is a entry point and is a crossing point when the start point is an
exit point. The end of the chain has the opposite crossing / bouncing status. x is the 
entry/exit status. =#
_inter_delay_cross_f(x) = (!x, x)
#= When marking the crossing status of a delayed bouncing, the chain start and end points
are crossing if the current polygon's adjacent edges are within the non-tracing polygon. If
the edges are outside then the chain endpoints are marked as bouncing. x is a boolean
representing if the edges are inside or outside of the polygon. =#
_inter_delay_bounce_f(x, _) = x
#= When tracing polygons, step forward if the most recent intersection point was an entry
point, else step backwards where x is the entry/exit status. =#
_inter_step(x, _) =  x ? 1 : (-1)

#= Polygon with multipolygon intersection - note that all intersection regions between
`poly_a` and any of the sub-polygons of `multipoly_b` are counted as intersection polygons.
Unless specified with `fix_multipoly = nothing`, `multipolygon_b` will be validated using
the given (default is `UnionIntersectingPolygons()`) correction. =#
function _intersection(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_b to prevent duplicated intersection regions
        multipoly_b = fix_multipoly(multipoly_b)
    end
    polys = Vector{_get_poly_type(T)}()
    for poly_b in GI.getpolygon(multipoly_b)
        append!(polys, intersection(alg, poly_a, poly_b; target))
    end
    return polys
end

#= Multipolygon with polygon intersection is equivalent to taking the intersection of the
polygon with the multipolygon and thus simply switches the order of operations and calls the
above method. =#
_intersection(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.PolygonTrait, poly_b;
    kwargs...,
) where T = intersection(alg, poly_b, multipoly_a; target , kwargs...)

#= Multipolygon with multipolygon intersection - note that all intersection regions between
any sub-polygons of `multipoly_a` and any of the sub-polygons of `multipoly_b` are counted
as intersection polygons. Unless specified with `fix_multipoly = nothing`, both 
`multipolygon_a` and `multipolygon_b` will be validated using the given (default is
`UnionIntersectingPolygons()`) correction. =#
function _intersection(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix both multipolygons to prevent duplicated regions
        multipoly_a = fix_multipoly(multipoly_a)
        multipoly_b = fix_multipoly(multipoly_b)
        fix_multipoly = nothing
    end
    polys = Vector{_get_poly_type(T)}()
    for poly_a in GI.getpolygon(multipoly_a)
        append!(polys, intersection(alg, poly_a, multipoly_b; target, fix_multipoly))
    end
    return polys
end
# catch-all method for multipolygontraits
function _intersection(
    alg::FosterHormannClipping, ::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    trait_a::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, polylike_a,
    trait_b::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, polylike_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...
) where T
    polys = _intersection(alg, TraitTarget(GI.PolygonTrait()), T, trait_a, polylike_a, trait_b, polylike_b; kwargs...)
    if isnothing(fix_multipoly)
        return GI.MultiPolygon(polys)
    else
        return fix_multipoly(GI.MultiPolygon(polys))
    end
end


# Many type and target combos aren't implemented
function _intersection(
    alg::GeometryOpsCore.Algorithm, target::TraitTarget{Target}, ::Type{T},
    trait_a::GI.AbstractTrait, geom_a,
    trait_b::GI.AbstractTrait, geom_b;
    kwargs...,
) where {Target, T}
    @assert(
        false,
        "Intersection between $trait_a and $trait_b with target $Target and algorithm $alg isn't implemented yet.",
    )
    return nothing
end

"""
    intersection_points(geom_a, geom_b, [T::Type])

Return a list of intersection tuple points between two geometries. If no intersection points
exist, returns an empty list.

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
inter_points = GO.intersection_points(line1, line2)

# output
1-element Vector{Tuple{Float64, Float64}}:
 (125.58375366067548, -14.83572303404496)
"""
intersection_points(geom_a, geom_b, ::Type{T} = Float64) where T <: AbstractFloat = intersection_points(FosterHormannClipping(Planar()), geom_a, geom_b, T)
function intersection_points(alg::FosterHormannClipping{M, A}, geom_a, geom_b, ::Type{T} = Float64) where {M, A, T <: AbstractFloat}
    return _intersection_points(alg.manifold, alg.accelerator, T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)
end

function intersection_points(m::Manifold, a::IntersectionAccelerator, geom_a, geom_b, ::Type{T} = Float64) where T <: AbstractFloat
    return _intersection_points(m, a, T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)
end


#= Calculates the list of intersection points between two geometries, including line
segments, line strings, linear rings, polygons, and multipolygons. =#
function _intersection_points(manifold::M, accelerator::A, ::Type{T}, ::GI.AbstractTrait, a, ::GI.AbstractTrait, b; exact = True()) where {M <: Manifold, A <: IntersectionAccelerator, T}
    # Initialize an empty list of points
    result = Tuple{T, T}[]
    # Check if the geometries extents even overlap
    Extents.intersects(GI.extent(a), GI.extent(b)) || return result
    # Create a list of edges from the two input geometries
    # edges_a, edges_b = map(sort! ∘ to_edges, (a, b))
    # Loop over pairs of edges and add any unique intersection points to results
    # TODO: add intersection acceleration here.

    function f_on_each_maybe_intersect((a_edge, a_idx), (b_edge, b_idx))
        line_orient, intr1, intr2 = _intersection_point(manifold, T, a_edge, b_edge; exact)
        line_orient == line_out && return LoopStateMachine.Action(:continue) # use LoopStateMachine.Continue() to skip this edge - in this case it doesn't matter but you could use it to e.g. break once you found the first intersecting point.
        pt1, _ = intr1
        push!(result, pt1)  # if not line_out, there is at least one intersection point
        if line_orient == line_over # if line_over, there are two intersection points
            pt2, _ = intr2
            push!(result, pt2)
        end
    end

    # iterate over each pair of intersecting edges only,
    # calling `f_on_each_maybe_intersect` for each pair 
    # that may intersect.
    foreach_pair_of_maybe_intersecting_edges_in_order(
        manifold, accelerator, 
        nothing, # f_on_each_a
        nothing, # f_after_each_a
        f_on_each_maybe_intersect, # f_on_each_maybe_intersect
        a,
        b,
        T
    )
    
    #= TODO: We might be able to just add unique points with checks on the α and β values
    returned from `_intersection_point`, but this would be different for curves vs polygons
    vs multipolygons depending on if the shape is closed. This then wouldn't allow using the
    `to_edges` functionality.  =# 
    unique!(sort!(result))
    return result
end

#= Calculates the intersection points between two lines if they exists and the fractional
component of each line from the initial end point to the intersection point where α is the
fraction along (a1, a2) and β is the fraction along (b1, b2).

Note that the first return is the type of intersection (line_cross, line_hinge, line_over,
or line_out). The type of intersection determines how many intersection points there are.
If the intersection is line_out, then there are no intersection points and the two
intersections aren't valid and shouldn't be used. If the intersection is line_cross or
line_hinge then the lines meet at one point and the first intersection is valid, while the
second isn't. Finally, if the intersection is line_over, then both points are valid and they
are the two points that define the endpoints of the overlapping region between the two
lines.

Also note again that each intersection is a tuple of two tuples. The first is the
intersection point (x,y) while the second is the ratio along the initial lines (α, β) for
that point. 

Calculation derivation can be found here: https://stackoverflow.com/questions/563198/ =#
function _intersection_point(manifold::M, ::Type{T}, (a1, a2)::Edge, (b1, b2)::Edge; exact) where {M <: Manifold, T}
    # Default answer for no intersection
    line_orient = line_out
    intr1 = ((zero(T), zero(T)), (zero(T), zero(T)))
    intr2 = intr1
    no_intr_result = (line_orient, intr1, intr2)
    # Seperate out line segment points
    (a1x, a1y), (a2x, a2y) = _tuple_point(a1, T), _tuple_point(a2, T)
    (b1x, b1y), (b2x, b2y) = _tuple_point(b1, T), _tuple_point(b2, T)
    # Check if envelopes of lines intersect
    a_ext = Extent(X = minmax(a1x, a2x), Y = minmax(a1y, a2y))
    b_ext = Extent(X = minmax(b1x, b2x), Y = minmax(b1y, b2y))
    !Extents.intersects(a_ext, b_ext) && return no_intr_result
    # Check orientation of two line segments with respect to one another
    a1_orient = Predicates.orient(b1, b2, a1; exact)
    a2_orient = Predicates.orient(b1, b2, a2; exact)
    a1_orient != 0 && a1_orient == a2_orient && return no_intr_result  # α < 0 or α > 1
    b1_orient = Predicates.orient(a1, a2, b1; exact)
    b2_orient = Predicates.orient(a1, a2, b2; exact)
    b1_orient != 0 && b1_orient == b2_orient && return no_intr_result  # β < 0 or β > 1
    # Determine intersection type and intersection point(s)
    if a1_orient == a2_orient == b1_orient == b2_orient == 0
        # Intersection is collinear if all endpoints lie on the same line
        line_orient, intr1, intr2 = _find_collinear_intersection(manifold, T, a1, a2, b1, b2, a_ext, b_ext, no_intr_result)
    elseif a1_orient == 0 || a2_orient == 0 || b1_orient == 0 || b2_orient == 0
        # Intersection is a hinge if the intersection point is an endpoint
        line_orient = line_hinge
        intr1 = _find_hinge_intersection(T, a1, a2, b1, b2, a1_orient, a2_orient, b1_orient)
    else
        # Intersection is a cross if there is only one non-endpoint intersection point
        line_orient = line_cross
        intr1 = _find_cross_intersection(T, a1, a2, b1, b2, a_ext, b_ext)
    end
    return line_orient, intr1, intr2
end

# TODO: deprecate this
_intersection_point(::Type{T}, (a1, a2)::Edge, (b1, b2)::Edge; exact) where T = _intersection_point(Planar(), T, (a1, a2), (b1, b2); exact)

#=
The same question between two great-circle arcs.

The method above takes a manifold but decides everything planar: its envelope test is a
lon/lat box, its four orientations are `Predicates.orient`, and the points it builds are
segment crossings in the chart. On the sphere all three are wrong for the same reason — an
edge of a DGG cell is a great-circle arc that bulges off the chart line joining its
endpoints — so this method replaces the whole decision rather than patching parts of it.

The classification is `_sph_arc_arc_class`, which is the identical `LineOrientation` split
the planar method produces (`line_cross` / `line_hinge` / `line_over` / `line_out`) and
additionally names which of the four endpoints lie on the other arc. No lon/lat envelope
pre-filter is applied: a great-circle arc leaves the box that bounds its endpoints, so such
a filter rejects real intersections near the poles, and it is unsound across the
antimeridian in either direction.

## Which coordinate is returned

Only `line_cross` builds a new point. Every other case meets at an endpoint that is already
a vertex of one of the two rings, and there the *input* tuple is returned unchanged rather
than a `lon/lat → xyz → lon/lat` round trip of it. That is load-bearing, not tidiness:
`_build_b_list` matches an intersection to a `b` vertex by `fracs[2] == 0`, and `equals`
compares `PolyNode` points for bit equality. A round trip moves a vertex by an ulp or two,
which would leave the match to fail and open a zero-length edge at exactly the shared
vertices a tiling is made of.

`exact` selects the orientation predicate the arc classification runs on.

Points come back in the representation they arrived in. On `Spherical` the edge iterator
hands FH `UnitSphericalPoint`s, so the whole path is xyz end to end and an input vertex is
returned as the very object that came in; the lon/lat method is kept for callers outside
the clipper and converts the one genuinely computed point on the way out.
=#
const _USPEdge{T} = Tuple{UnitSpherical.UnitSphericalPoint{T}, UnitSpherical.UnitSphericalPoint{T}}

_intersection_point(m::Spherical, ::Type{T}, a::_USPEdge, b::_USPEdge; exact) where {T} =
    _sph_intersection_point(m, T, a, b; exact)
_intersection_point(m::Spherical, ::Type{T}, a::Edge, b::Edge; exact) where {T} =
    _sph_intersection_point(m, T, a, b; exact)

#-- Give a computed crossing back in the representation the edge arrived in.
_as_ingested(x, ::Tuple, ::Type{T}) where {T} = _sph_lonlat(T, x)
_as_ingested(x, ::UnitSpherical.UnitSphericalPoint, ::Type{T}) where {T} = x

function _sph_intersection_point(m::Spherical, ::Type{T}, (a1, a2), (b1, b2); exact) where {T}
    #= The sentinel point is never read: `line_out` carries no point, and `intr2` is only
    destructured on the `line_over` branch, which fills it. Reusing `a1` keeps the returned
    tuple concretely typed in whichever representation came in. =#
    zero_intr = (a1, (zero(T), zero(T)))
    no_intr_result = (line_out, zero_intr, zero_intr)

    A0 = _spherical_kernel_point(a1); A1 = _spherical_kernel_point(a2)
    B0 = _spherical_kernel_point(b1); B1 = _spherical_kernel_point(b2)

    #= Clipping asks for the exact predicate. The banded `spherical_orient` reports 0
    inside an eps*16 window, which at DGG cell scale is wider than the determinant it is
    judging: two arcs that genuinely cross are then classified as a hinge, the entry/exit
    alternation collapses, and the tracer emits the whole subject ring. =#
    orient, a0_on_b, a1_on_b, b0_on_a, b1_on_a = _sph_arc_arc_class(A0, A1, B0, B1, exact)
    orient === line_out && return no_intr_result

    if orient === line_cross
        x = _arc_crossing_point(A0, A1, B0, B1)
        #= A proper crossing cannot produce parallel normals, but the constructor guards
        the normalization anyway and so does this. =#
        x === nothing && return no_intr_result
        α = _sph_arc_frac(T, A0, A1, x)
        β = _sph_arc_frac(T, B0, B1, x)
        return line_cross, (_as_ingested(x, a1, T), (α, β)), zero_intr
    end

    #= Both remaining cases meet at named endpoints. Building the four candidates up front
    keeps one spelling of "which point, and where does it sit on each arc"; the tuple is
    homogeneous and stack-allocated, so this costs nothing. =#
    cands = (
        (a0_on_b, _tuple_point(a1, T), zero(T), _sph_arc_frac(T, B0, B1, A0)),
        (a1_on_b, _tuple_point(a2, T), one(T),  _sph_arc_frac(T, B0, B1, A1)),
        (b0_on_a, _tuple_point(b1, T), _sph_arc_frac(T, A0, A1, B0), zero(T)),
        (b1_on_a, _tuple_point(b2, T), _sph_arc_frac(T, A0, A1, B1), one(T)),
    )

    if orient === line_hinge
        #= One distinct meeting point, possibly named by several incidences (a shared
        vertex is on both arcs and is reported four times). Any of them is that point. =#
        @inbounds for k in 1:4
            on, p, α, β = cands[k]
            on && return line_hinge, (p, (α, β)), zero_intr
        end
        return no_intr_result
    end

    #= `line_over`: the arcs share a great circle and overlap with positive length. The
    overlap runs between the two extreme incident endpoints, ordered along `a` — the same
    ordering convention the planar method's `intr1`/`intr2` follow. =#
    lo_k, hi_k = 0, 0
    lo_α, hi_α = T(Inf), T(-Inf)
    @inbounds for k in 1:4
        on, _, α, _ = cands[k]
        on || continue
        if α < lo_α; lo_α = α; lo_k = k end
        if α > hi_α; hi_α = α; hi_k = k end
    end
    (lo_k == 0 || lo_k == hi_k) && return no_intr_result
    @inbounds begin
        _, p_lo, α_lo, β_lo = cands[lo_k]
        _, p_hi, α_hi, β_hi = cands[hi_k]
    end
    return line_over, (p_lo, (α_lo, β_lo)), (p_hi, (α_hi, β_hi))
end

#= Arc-length fraction of `x` along the arc `p0 → p1`, clamped to `[0, 1]`.

## The endpoints are settled by identity, never by arithmetic

`_build_a_list` claims a shared vertex with exact comparisons — `α == 0`, `β == 0`,
`0 ≤ β < 1` — so a fraction that is *meant* to be an endpoint has to be exactly `0` or
exactly `1`, not within a tolerance of it. Recovering that by dividing two separately
computed angles does not deliver it: for a vertex shared between two rings the two `atan`
calls take bit-identical inputs, yet inlined into the four-candidate construction they do
not contract identically, and the quotient lands an ulp below one. `0.9999999999999998`
passes `β < 1`, so the vertex is claimed a second time and the traversal is handed a
duplicate intersection — which is what `_trace_polynodes!` reports as a `TracingError`.

Comparing against the endpoints first removes the arithmetic from the case that has to be
exact. `_spherical_kernel_point` is deterministic and signed-zero normalized, so the same
lon/lat vertex always converts to the same bits, and shared vertices — the common case in a
tiling — take this path rather than the quotient.

A zero-length arc is caught by the same test: `x == p0` fires, and `0` is the only answer
that keeps the caller's endpoint tests meaningful.

## Why `atan`, not `acos`

`atan(‖a × b‖, a ⋅ b)` stays accurate across the whole range, where `acos(a ⋅ b)` loses
about half its bits for the small angles a cell-scale edge subtends — a HEALPix level-18
edge is ~4e-6 rad. =#
@inline function _sph_arc_frac(::Type{T}, p0, p1, x) where {T}
    x == p0 && return zero(T)
    x == p1 && return one(T)
    total = atan(norm(cross(p0, p1)), dot(p0, p1))
    total == 0 && return zero(T)
    part = atan(norm(cross(p0, x)), dot(p0, x))
    return T(clamp(part / total, 0, 1))
end

@inline _sph_lonlat(::Type{T}, u) where {T} = ((ll = _usp_to_lonlat(u)); (T(ll[1]), T(ll[2])))

#= If lines defined by (a1, a2) and (b1, b2) are collinear, find endpoints of overlapping
region if they exist. This could result in three possibilities. First, there could be no
overlapping region, in which case, the default 'no_intr_result' intersection information is
returned. Second, the two regions could just meet at one shared endpoint, in which case it
is a hinge intersection with one intersection point. Otherwise, it is a overlapping
intersection defined by two of the endpoints of the line segments. =#
function _find_collinear_intersection(manifold::M, ::Type{T}, a1, a2, b1, b2, a_ext, b_ext, no_intr_result) where {M <: Manifold, T}
    # Define default return for no intersection points
    line_orient, intr1, intr2 = no_intr_result
    # Determine collinear line overlaps
    a1_in_b = _point_in_extent(a1, b_ext)
    a2_in_b = _point_in_extent(a2, b_ext)
    b1_in_a = _point_in_extent(b1, a_ext)
    b2_in_a = _point_in_extent(b2, a_ext)
    # Determine line distances
    a_dist, b_dist = distance(a1, a2, T), distance(b1, b2, T)
    # Set collinear intersection points if they exist
    if a1_in_b && a2_in_b      # 1st vertex of a and 2nd vertex of a form overlap
        line_orient = line_over
        β1 = _clamped_frac(distance(a1, b1, T), b_dist)
        β2 = _clamped_frac(distance(a2, b1, T), b_dist)
        intr1 = (_tuple_point(a1, T), (zero(T), β1))
        intr2 = (_tuple_point(a2, T), (one(T), β2))
    elseif b1_in_a && b2_in_a  # 1st vertex of b and 2nd vertex of b form overlap
        line_orient = line_over
        α1 = _clamped_frac(distance(b1, a1, T), a_dist)
        α2 = _clamped_frac(distance(b2, a1, T), a_dist)
        intr1 = (_tuple_point(b1, T), (α1, zero(T)))
        intr2 = (_tuple_point(b2, T), (α2, one(T)))
    elseif a1_in_b && b1_in_a  # 1st vertex of a and 1st vertex of b form overlap
        if equals(a1, b1)
            line_orient = line_hinge
            intr1 = (_tuple_point(a1, T), (zero(T), zero(T)))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(T, a1, b1, zero(T), zero(T), a1, b1, a_dist, b_dist)
        end
    elseif a1_in_b && b2_in_a  # 1st vertex of a and 2nd vertex of b form overlap
        if equals(a1, b2)
            line_orient = line_hinge
            intr1 = (_tuple_point(a1, T), (zero(T), one(T)))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(T, a1, b2, zero(T), one(T), a1, b1, a_dist, b_dist) 
        end
    elseif a2_in_b && b1_in_a  # 2nd vertex of a and 1st vertex of b form overlap
        if equals(a2, b1)
            line_orient = line_hinge
            intr1 = (_tuple_point(a2, T), (one(T), zero(T)))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(T, a2, b1, one(T), zero(T), a1, b1, a_dist, b_dist)
        end
    elseif a2_in_b && b2_in_a  # 2nd vertex of a and 2nd vertex of b form overlap
        if equals(a2, b2)
            line_orient = line_hinge
            intr1 = (_tuple_point(a2, T), (one(T), one(T)))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(T, a2, b2, one(T), one(T), a1, b1, a_dist, b_dist)
        end
    end
    return line_orient, intr1, intr2
end

#= Determine intersection points and segment fractions when overlap is made up one one 
endpoint of segment (a1, a2) and one endpoint of segment (b1, b2). =#
_set_ab_collinear_intrs(::Type{T}, a_pt, b_pt, a_pt_α, b_pt_β, a1, b1, a_dist, b_dist) where T =
    (
        (_tuple_point(a_pt, T), (a_pt_α, _clamped_frac(distance(a_pt, b1, T), b_dist))),
        (_tuple_point(b_pt, T), (_clamped_frac(distance(b_pt, a1, T), a_dist), b_pt_β))
    )

#= If lines defined by (a1, a2) and (b1, b2) are just touching at one of those endpoints and
are not collinear, then they form a hinge, with just that one shared intersection point.
Point equality is checked before segment orientation to have maximal accurary on fractions
to avoid floating point errors. If the points are not equal, we know that the hinge does not
take place at an endpoint and the fractions must be between 0 or 1 (exclusive). =#
function _find_hinge_intersection(::Type{T}, a1, a2, b1, b2, a1_orient, a2_orient, b1_orient) where T
    pt, α, β = if equals(a1, b1)
        _tuple_point(a1, T), zero(T), zero(T)
    elseif equals(a1, b2)
        _tuple_point(a1, T), zero(T), one(T)
    elseif equals(a2, b1)
        _tuple_point(a2, T), one(T), zero(T)
    elseif equals(a2, b2)
        _tuple_point(a2, T), one(T), one(T)
    elseif a1_orient == 0
        β_val = _clamped_frac(distance(b1, a1, T), distance(b1, b2, T), eps(T))
        _tuple_point(a1, T), zero(T), β_val
    elseif a2_orient == 0
        β_val = _clamped_frac(distance(b1, a2, T), distance(b1, b2, T), eps(T))
        _tuple_point(a2, T), one(T), β_val
    elseif b1_orient == 0
        α_val = _clamped_frac(distance(a1, b1, T), distance(a1, a2, T), eps(T))
        _tuple_point(b1, T), α_val, zero(T)
    else  # b2_orient == 0
        α_val = _clamped_frac(distance(a1, b2, T), distance(a1, a2, T), eps(T))
        _tuple_point(b2, T), α_val, one(T)
    end
    return pt, (α, β)
end

#= If lines defined by (a1, a2) and (b1, b2) meet at one point that is not an endpoint of
either segment, they form a crossing intersection with a singular intersection point. That 
point is calculated by finding the fractional distance along each segment the point occurs
at (α, β). If the point is too close to an endpoint to be distinct, the point shares a value
with the endpoint, but with a non-zero and non-one fractional value. If the intersection
point calculated is outside of the envelope of the two segments due to floating point error,
it is set to the endpoint of the two segments that is closest to the other segment.
Regardless of point value, we know that it does not actually occur at an endpoint so the
fractions must be between 0 or 1 (exclusive). =#
function _find_cross_intersection(::Type{T}, a1, a2, b1, b2, a_ext, b_ext) where T
    # First line runs from a to a + Δa
    (a1x, a1y), (a2x, a2y) = _tuple_point(a1, T), _tuple_point(a2, T)
    Δax, Δay = a2x - a1x, a2y - a1y
    # Second line runs from b to b + Δb 
    (b1x, b1y), (b2x, b2y) = _tuple_point(b1, T), _tuple_point(b2, T)
    Δbx, Δby = b2x - b1x, b2y - b1y
    # Differences between starting points
    Δbax = b1x - a1x
    Δbay = b1y - a1y
    a_cross_b = Δax * Δby - Δay * Δbx
    # Determine α value where 0 < α < 1 and β value where 0 < β < 1
    α = _clamped_frac(Δbax * Δby - Δbay * Δbx, a_cross_b, eps(T))
    β = _clamped_frac(Δbax * Δay - Δbay * Δax, a_cross_b, eps(T))

    #= Intersection will be where a1 + α * Δa = b1 + β * Δb. However, due to floating point
    inaccuracies, α and β calculations may yield different intersection points. Average
    both points together to minimize difference from real value, as long as segment isn't 
    vertical or horizontal as this will almost certainly lead to the point being outside the
    envelope due to floating point error. Also note that floating point limitations could
    make intersection be endpoint if α≈0 or α≈1.=#
    x = if Δax == 0
        a1x
    elseif Δbx == 0
        b1x
    else
        (a1x + α * Δax + b1x + β * Δbx) / 2
    end
    y = if Δay == 0
        a1y
    elseif Δby == 0
        b1y
    else
        (a1y + α * Δay + b1y + β * Δby) / 2
    end
    pt = (x, y)
    # Check if point is within segment envelopes and adjust to endpoint if not
    if !_point_in_extent(pt, a_ext) || !_point_in_extent(pt, b_ext)
        pt, α, β = _nearest_endpoint(T, a1, a2, b1, b2)
    end
    return (pt, (α, β))
end

# Find endpoint of either segment that is closest to the opposite segment
function _nearest_endpoint(::Type{T}, a1, a2, b1, b2) where T
    # Create lines from segments and calculate segment length
    a_line, a_dist = GI.Line(StaticArrays.SVector(a1, a2)), distance(a1, a2, T)
    b_line, b_dist = GI.Line(StaticArrays.SVector(b1, b2)), distance(b1, b2, T)
    # Determine distance from a1 to segment b
    min_pt, min_dist = a1, distance(a1, b_line, T)
    α, β = eps(T), _clamped_frac(distance(min_pt, b1, T), b_dist, eps(T))
    # Determine distance from a2 to segment b
    dist = distance(a2, b_line, T)
    if dist < min_dist
        min_pt, min_dist = a2, dist
        α, β = one(T) - eps(T), _clamped_frac(distance(min_pt, b1, T), b_dist, eps(T))
    end
    # Determine distance from b1 to segment a
    dist = distance(b1, a_line, T)
    if dist < min_dist
        min_pt, min_dist = b1, dist
        α, β = _clamped_frac(distance(min_pt, a1, T), a_dist, eps(T)), eps(T)
    end
    # Determine distance from b2 to segment a
    dist = distance(b2, a_line, T)
    if dist < min_dist
        min_pt, min_dist = b2, dist
        α, β = _clamped_frac(distance(min_pt, a2, T), a_dist, eps(T)), one(T) - eps(T)
    end
    # Return point with smallest distance
    return _tuple_point(min_pt, T), α, β
end

# Return value of x/y clamped between ϵ and 1 - ϵ
_clamped_frac(x::T, y::T, ϵ = zero(T)) where T = clamp(x / y, ϵ, one(T) - ϵ)
