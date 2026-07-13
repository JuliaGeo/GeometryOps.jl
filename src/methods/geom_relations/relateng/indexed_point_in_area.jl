# # RelateNG indexed point-in-area location
#
# Prepared-mode point-in-area locator for RelateNG (Task 22). This file holds
# the ports of the two JTS classes behind prepared-mode point location
# (JTS file boundaries preserved as clearly marked sections):
#
# 1. `RayCrossingCounter`         (JTS algorithm/RayCrossingCounter.java)
# 2. `IndexedPointInAreaLocator`  (JTS algorithm/locate/IndexedPointInAreaLocator.java)
#
# `RelatePointLocator` (point_locator.jl) swaps this locator in for the
# SimplePointInAreaLocator ring loop when `is_prepared` is set, mirroring
# JTS `RelatePointLocator.getLocator`.
#
# JTS backs the locator with its 1D `SortedPackedIntervalRTree` over
# segment y-intervals; here that role is played by the shared
# `RTree(STR(), ...)` over 1-D `(Y,)` extents — sort-tile-recursive in one
# dimension IS the midpoint sort of JTS's `NodeComparator`, so the packed
# layout is the same idea with a wider fanout. The query is inherently
# 1-dimensional: the horizontal ray from the test point must visit *every*
# segment whose y-interval contains `p.y`, regardless of x (segments wholly
# left of the point are rejected inside `count_segment!`, exactly as in
# JTS), so a 2D index could not prune more candidates without changing the
# ray-crossing counting contract. The midpoint sort earns its cost in this
# index's only (prepared, build-once-query-forever) use: it groups same-`y`
# segments so a point query descends few subtrees, where insertion (ring)
# order recrosses the query `y` in many separated runs — measured on
# Natural Earth 10m Canada, ring-order queries are ~3× slower.

#==========================================================================
## RayCrossingCounter (port of JTS RayCrossingCounter.java)
==========================================================================#

"""
    RayCrossingCounter(m::Manifold, p; exact)

Counts the number of segments crossed by a horizontal ray extending to the
right from a given point, in an incremental fashion. This can be used to
determine whether a point lies in a polygonal geometry. The class determines
the situation where the point lies exactly on a segment. This handles
polygonal geometries with any number of shells and holes; ring orientation
is unimportant. In order to compute a correct location for a given polygonal
geometry, it is essential that **all** segments are counted which touch the
ray or lie in any ring which may contain the point — which is what allows
optimization by y-interval indexing, since segments whose y-extent misses
the ray a priori cannot touch it.

The manifold `m` and the `exact` flag are stored in the struct (consistent
with [`AdjacentEdgeLocator`](@ref)); the orientation test goes through
`rk_orient` (JTS uses the extended-precision `Orientation.index`, matching
`exact = True()`). The horizontal-ray sweep itself is coordinate-plane
logic, exactly as in JTS.
"""
mutable struct RayCrossingCounter{M <: Manifold, E}
    const m::M
    const exact::E
    const p::Tuple{Float64, Float64}
    crossing_count::Int
    # true if the test point lies on an input segment
    is_point_on_segment::Bool
end

RayCrossingCounter(m::Manifold, p; exact) =
    RayCrossingCounter(m, exact, _node_point(p), 0, false)

"""
    count_segment!(rcc::RayCrossingCounter, p1, p2)

Counts a segment with endpoints `p1`, `p2`. Port of
`RayCrossingCounter.countSegment`.
"""
function count_segment!(rcc::RayCrossingCounter, p1, p2)
    px, py = rcc.p
    p1x, p1y = GI.x(p1), GI.y(p1)
    p2x, p2y = GI.x(p2), GI.y(p2)
    #=
    For each segment, check if it crosses a horizontal ray running from the
    test point in the positive x direction.
    =#
    #-- check if the segment is strictly to the left of the test point
    (p1x < px && p2x < px) && return nothing
    #-- check if the point is equal to the current ring vertex
    if px == p2x && py == p2y
        rcc.is_point_on_segment = true
        return nothing
    end
    #=
    For horizontal segments, check if the point is on the segment.
    Otherwise, horizontal segments are not counted.
    =#
    if p1y == py && p2y == py
        minx, maxx = minmax(p1x, p2x)
        if minx <= px <= maxx
            rcc.is_point_on_segment = true
        end
        return nothing
    end
    #=
    Evaluate all non-horizontal segments which cross a horizontal ray to the
    right of the test pt. To avoid double-counting shared vertices, we use
    the convention that
    - an upward edge includes its starting endpoint, and excludes its final
      endpoint
    - a downward edge excludes its starting endpoint, and includes its final
      endpoint
    =#
    if (p1y > py && p2y <= py) || (p2y > py && p1y <= py)
        orient = rk_orient(rcc.m, (p1x, p1y), (p2x, p2y), rcc.p; exact = rcc.exact)
        if orient == 0
            rcc.is_point_on_segment = true
            return nothing
        end
        #-- re-orient the result if needed to ensure effective segment direction is upwards
        if p2y < p1y
            orient = -orient
        end
        #-- the upward segment crosses the ray if the test point lies to the left (CCW) of the segment
        if orient > 0
            rcc.crossing_count += 1
        end
    end
    return nothing
end

# Port of RayCrossingCounter.isOnSegment: whether the point lies exactly on
# one of the supplied segments. May be checked at any time as segments are
# processed; once true, the result never changes again.
is_on_segment(rcc::RayCrossingCounter) = rcc.is_point_on_segment

# Port of RayCrossingCounter.getLocation: the `LOC_*` location of the point
# relative to the ring, polygon or multipolygon from which the processed
# segments were provided. Only correct once all relevant segments have been
# processed.
function rcc_location(rcc::RayCrossingCounter)
    rcc.is_point_on_segment && return LOC_BOUNDARY
    #-- the point is in the interior of the ring if the number of X-crossings is odd
    return isodd(rcc.crossing_count) ? LOC_INTERIOR : LOC_EXTERIOR
end

#==========================================================================
## IndexedPointInAreaLocator (port of JTS IndexedPointInAreaLocator.java)
==========================================================================#

# Leaf item of the segment index: a ring segment as a pair of node points.
const _PIASegment = Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}
# The segment index: a midpoint-sorted packed tree over the segments'
# y-intervals (see the header note on how this maps to JTS's
# SortedPackedIntervalRTree).
const _PIAExtent = Extents.Extent{(:Y,), Tuple{NTuple{2, Float64}}}
const _PIAIndex = RTree{STR, _PIAExtent, Vector{_PIASegment}, Vector{Int}}

"""
    IndexedPointInAreaLocator(m::Manifold, geom; exact)

Determines the location (`LOC_*` code) of points relative to an areal
geometry, using indexing for efficiency. This algorithm is suitable for use
in cases where many points will be tested against a given area. The location
is computed precisely: points located on the geometry boundary or segments
return `LOC_BOUNDARY`.

Port of JTS `IndexedPointInAreaLocator` together with its private
`IntervalIndexedGeometry` (the y-interval segment index; its `isEmpty`
flag is `index === nothing` here, since a recursively empty polygonal
geometry contributes no rings, hence no segments). JTS lazy-loads the
index on the first `locate`; here the index is built in the constructor,
since `RelatePointLocator` already creates the locator itself lazily on
the first use per polygonal element (`_get_poly_locator`, the port of
`RelatePointLocator.getLocator`).
"""
struct IndexedPointInAreaLocator{M <: Manifold, E}
    m::M
    exact::E
    index::Union{Nothing, _PIAIndex}
end

function IndexedPointInAreaLocator(m::Manifold, geom; exact)
    exts = _PIAExtent[]
    segs = _PIASegment[]
    n = GI.npoint(geom)
    sizehint!(exts, n); sizehint!(segs, n)
    _interval_index_add_geom!(exts, segs, GI.trait(geom), geom)
    index = isempty(segs) ? nothing : RTree(STR(), segs; extents = exts)
    return IndexedPointInAreaLocator(m, exact, index)
end

"""
    locate(loc::IndexedPointInAreaLocator, p)

The location (`LOC_*` code) of point `p` in the locator's areal geometry.
Port of `IndexedPointInAreaLocator.locate`.
"""
function locate(loc::IndexedPointInAreaLocator, p)
    index = loc.index
    index === nothing && return LOC_EXTERIOR   #-- IntervalIndexedGeometry.isEmpty
    rcc = RayCrossingCounter(loc.m, p; exact = loc.exact)
    y = rcc.p[2]
    ray = Extents.Extent(Y = (y, y))
    #-- SegmentVisitor: count every segment whose y-interval touches the ray
    SpatialTreeInterface.depth_first_search(Base.Fix1(Extents.intersects, ray), index) do i
        seg = index.data[i]
        count_segment!(rcc, seg[1], seg[2])
    end
    return rcc_location(rcc)
end

# Port of IntervalIndexedGeometry.init. JTS extracts all linear components
# (LinearComponentExtracter) and keeps the closed ones; here only polygonal
# elements ever reach this locator (RelatePointLocator extracts Polygon /
# MultiPolygon elements), so the rings are iterated directly.
function _interval_index_add_geom!(exts, segs, ::GI.PolygonTrait, poly)
    GI.isempty(poly) && return nothing
    _interval_index_add_line!(exts, segs, GI.getexterior(poly))
    for hole in GI.gethole(poly)
        _interval_index_add_line!(exts, segs, hole)
    end
    return nothing
end

function _interval_index_add_geom!(exts, segs, ::GI.MultiPolygonTrait, mp)
    for poly in GI.getgeom(mp)
        _interval_index_add_geom!(exts, segs, GI.trait(poly), poly)
    end
    return nothing
end

# Port of IntervalIndexedGeometry.addLine: index each ring segment on its
# y-interval, streaming the points directly (no `_node_points` copy).
# GI rings may be implicitly closed (no repeated end point); the
# SimplePointInAreaLocator ring loop (`rk_point_in_ring`) treats rings as
# closed regardless, so the closing segment is added here too.
function _interval_index_add_line!(exts, segs, ring)
    n = GI.npoint(ring)
    n < 2 && return nothing
    first_pt = _node_point(GI.getpoint(ring, 1))
    prev = first_pt
    for i in 2:n
        pt = _node_point(GI.getpoint(ring, i))
        _interval_index_add_segment!(exts, segs, prev, pt)
        prev = pt
    end
    if prev != first_pt
        _interval_index_add_segment!(exts, segs, prev, first_pt)
    end
    return nothing
end

function _interval_index_add_segment!(exts, segs, p0, p1)
    push!(exts, Extents.Extent(Y = minmax(p0[2], p1[2])))
    push!(segs, (p0, p1))
    return nothing
end
