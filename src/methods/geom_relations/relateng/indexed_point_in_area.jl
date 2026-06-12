# # RelateNG indexed point-in-area location
#
# Prepared-mode point-in-area locator for RelateNG (Task 22). This file holds
# the ports of the three JTS classes behind prepared-mode point location, in
# this order (JTS file boundaries preserved as clearly marked sections):
#
# 1. `SortedPackedIntervalRTree`  (JTS index/intervalrtree/SortedPackedIntervalRTree.java)
# 2. `RayCrossingCounter`         (JTS algorithm/RayCrossingCounter.java)
# 3. `IndexedPointInAreaLocator`  (JTS algorithm/locate/IndexedPointInAreaLocator.java)
#
# `RelatePointLocator` (point_locator.jl) swaps this locator in for the
# SimplePointInAreaLocator ring loop when `is_prepared` is set, mirroring
# JTS `RelatePointLocator.getLocator`.
#
# Indexing choice: this ports the JTS 1D `SortedPackedIntervalRTree` over
# segment y-intervals rather than reusing the existing 2D `STRtree`
# machinery (`_make_prepared_edge_index` in relate_ng.jl). The query here is
# inherently 1-dimensional: the horizontal ray from the test point must
# visit *every* segment whose y-interval contains `p.y`, regardless of x
# (segments wholly left of the point are rejected inside `count_segment!`,
# exactly as in JTS), so a 2D index could not prune more candidates without
# changing the ray-crossing counting contract — and the packed 1D tree is
# smaller, cheaper to build, and queries without allocating.

#==========================================================================
## SortedPackedIntervalRTree (port of JTS SortedPackedIntervalRTree.java)
==========================================================================#

"""
    SortedPackedIntervalRTree(mins, maxs, items)

A static index on a set of 1-dimensional intervals, using an R-Tree packed
based on the order of the interval midpoints. It supports range searching,
where the range is an interval of the real line (which may be a single
point). A common use is to index 1-dimensional intervals which are the
projection of 2-D objects onto an axis of the coordinate system.

Port of JTS `SortedPackedIntervalRTree`, with two representation changes
(behavior, tree shape and query order are identical):

- JTS builds the tree lazily from incremental `insert` calls on the first
  query; the index is static once queried, so here the constructor takes all
  the intervals at once and packs eagerly.
- JTS builds an object tree of branch/leaf nodes (`IntervalRTreeNode` and
  subclasses); an abstractly-typed node field would box in Julia, so the
  packed tree is stored as flat per-level extent arrays instead: level 1 is
  the leaves, and node `j` of level `k + 1` covers nodes `2j - 1` and `2j`
  of level `k` (an unpaired trailing node is carried up unchanged, as in
  `buildLevel`). The last level is the root.
- JTS always sorts the leaves by interval midpoint (`NodeComparator`)
  before packing. Here `sort_leaves = false` skips that and keeps insertion
  (ring) order — the `NaturalIndexing` observation. Query results are
  order-independent (every visit is extent-checked), but the layouts trade
  off differently: midpoint order groups same-`y` segments so a point query
  descends few subtrees, while a long coastline ring in natural order
  recrosses the query `y` in many separated runs. Measured on Natural Earth
  10m Canada: the sort is ~4× the rest of the build, and natural-order
  queries are ~3× slower. So prepared mode sorts (build once, query
  forever) and the lazily indexed unprepared path doesn't (its query count
  is at most a few hundred, far below the ~1000-query crossover; see
  `locate_on_polygonal`).
"""
struct SortedPackedIntervalRTree{I}
    # leaf items: midpoint-sorted (`sort_leaves = true`) or insertion order
    items::Vector{I}
    # level_min[1][i] / level_max[1][i] is the interval of leaf item i;
    # level k > 1 holds the pairwise-combined extents of level k - 1
    level_min::Vector{Vector{Float64}}
    level_max::Vector{Vector{Float64}}
end

# Port of insert + init/buildRoot/buildTree/buildLevel, packed eagerly.
# Without `sort_leaves` the leaf arrays are taken over by the tree, not
# copied.
function SortedPackedIntervalRTree(mins::Vector{Float64}, maxs::Vector{Float64},
        items::Vector{I}; sort_leaves::Bool = true) where {I}
    if sort_leaves
        #-- sort the leaf nodes (IntervalRTreeNode.NodeComparator: by
        #-- midpoint; sortperm is stable, matching Collections.sort)
        n = length(items)
        perm = sortperm(Float64[(mins[i] + maxs[i]) / 2 for i in 1:n])
        mins = mins[perm]
        maxs = maxs[perm]
        items = items[perm]
    end
    level_min = [mins]
    level_max = [maxs]
    #-- now group nodes into blocks of two and build tree up recursively
    while length(level_min[end]) > 1
        src_min = level_min[end]
        src_max = level_max[end]
        nsrc = length(src_min)
        ndest = cld(nsrc, 2)
        dest_min = Vector{Float64}(undef, ndest)
        dest_max = Vector{Float64}(undef, ndest)
        for j in 1:ndest
            i = 2j - 1
            if i + 1 <= nsrc
                #-- IntervalRTreeBranchNode.buildExtent
                dest_min[j] = min(src_min[i], src_min[i + 1])
                dest_max[j] = max(src_max[i], src_max[i + 1])
            else
                #-- unpaired trailing node is carried up unchanged
                dest_min[j] = src_min[i]
                dest_max[j] = src_max[i]
            end
        end
        push!(level_min, dest_min)
        push!(level_max, dest_max)
    end
    return SortedPackedIntervalRTree{I}(items, level_min, level_max)
end

"""
    query_interval(f, tree::SortedPackedIntervalRTree, qmin, qmax)

Search for intervals in the index which intersect the given closed interval
`[qmin, qmax]` and apply the function `f` to each matched item. Port of
`SortedPackedIntervalRTree.query` with the `ItemVisitor` replaced by a
function (typically a `do`-block closure).
"""
function query_interval(f::F, tree::SortedPackedIntervalRTree, qmin::Float64, qmax::Float64) where {F}
    #-- if there are no leaves the tree is empty (Java: root == null)
    isempty(tree.items) && return nothing
    _irt_query(f, tree, length(tree.level_min), 1, qmin, qmax)
    return nothing
end

# Port of IntervalRTreeBranchNode.query / IntervalRTreeLeafNode.query over
# the packed levels: node `i` of `level`, recursing down to the leaves.
function _irt_query(f::F, tree::SortedPackedIntervalRTree, level::Int, i::Int, qmin::Float64, qmax::Float64) where {F}
    #-- IntervalRTreeNode.intersects
    (tree.level_min[level][i] > qmax || tree.level_max[level][i] < qmin) && return nothing
    if level == 1
        #-- leaf node: visit the item
        f(tree.items[i])
    else
        #-- branch node: query both children
        child = 2i - 1
        _irt_query(f, tree, level - 1, child, qmin, qmax)
        if child + 1 <= length(tree.level_min[level - 1])
            _irt_query(f, tree, level - 1, child + 1, qmin, qmax)
        end
    end
    return nothing
end

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

"""
    IndexedPointInAreaLocator(m::Manifold, geom; exact)

Determines the location (`LOC_*` code) of points relative to an areal
geometry, using indexing for efficiency. This algorithm is suitable for use
in cases where many points will be tested against a given area. The location
is computed precisely: points located on the geometry boundary or segments
return `LOC_BOUNDARY`.

Port of JTS `IndexedPointInAreaLocator` together with its private
`IntervalIndexedGeometry` (the `is_empty` flag and the y-interval segment
index). JTS lazy-loads the index on the first `locate`; here the index is
built in the constructor, since `RelatePointLocator` already creates the
locator itself lazily on the first use per polygonal element
(`_get_poly_locator`, the port of `RelatePointLocator.getLocator`).
"""
struct IndexedPointInAreaLocator{M <: Manifold, E}
    m::M
    exact::E
    index::SortedPackedIntervalRTree{_PIASegment}
    is_empty::Bool
end

function IndexedPointInAreaLocator(m::Manifold, geom; exact, sort_leaves::Bool = true)
    mins = Float64[]
    maxs = Float64[]
    segs = _PIASegment[]
    n = GI.npoint(geom)
    sizehint!(mins, n); sizehint!(maxs, n); sizehint!(segs, n)
    _iig_add_geom!(mins, maxs, segs, GI.trait(geom), geom)
    index = SortedPackedIntervalRTree(mins, maxs, segs; sort_leaves)
    #-- IntervalIndexedGeometry.isEmpty: a (recursively) empty polygonal
    #-- geometry contributes no rings, hence no segments
    return IndexedPointInAreaLocator(m, exact, index, isempty(segs))
end

"""
    locate(loc::IndexedPointInAreaLocator, p)

The location (`LOC_*` code) of point `p` in the locator's areal geometry.
Port of `IndexedPointInAreaLocator.locate`.
"""
function locate(loc::IndexedPointInAreaLocator, p)
    loc.is_empty && return LOC_EXTERIOR
    rcc = RayCrossingCounter(loc.m, p; exact = loc.exact)
    y = rcc.p[2]
    #-- SegmentVisitor: count every segment whose y-interval touches the ray
    query_interval(loc.index, y, y) do seg
        count_segment!(rcc, seg[1], seg[2])
    end
    return rcc_location(rcc)
end

# Port of IntervalIndexedGeometry.init. JTS extracts all linear components
# (LinearComponentExtracter) and keeps the closed ones; here only polygonal
# elements ever reach this locator (RelatePointLocator extracts Polygon /
# MultiPolygon elements), so the rings are iterated directly.
function _iig_add_geom!(mins, maxs, segs, ::GI.PolygonTrait, poly)
    GI.isempty(poly) && return nothing
    _iig_add_line!(mins, maxs, segs, GI.getexterior(poly))
    for hole in GI.gethole(poly)
        _iig_add_line!(mins, maxs, segs, hole)
    end
    return nothing
end

function _iig_add_geom!(mins, maxs, segs, ::GI.MultiPolygonTrait, mp)
    for poly in GI.getgeom(mp)
        _iig_add_geom!(mins, maxs, segs, GI.trait(poly), poly)
    end
    return nothing
end

# Port of IntervalIndexedGeometry.addLine: index each ring segment on its
# y-interval, streaming the points directly (no `_node_points` copy — this
# runs on the unprepared hot path). GI rings may be implicitly closed (no
# repeated end point); the SimplePointInAreaLocator ring loop
# (`rk_point_in_ring`) treats rings as closed regardless, so the closing
# segment is added here too.
function _iig_add_line!(mins, maxs, segs, ring)
    n = GI.npoint(ring)
    n < 2 && return nothing
    first_pt = _node_point(GI.getpoint(ring, 1))
    prev = first_pt
    for i in 2:n
        pt = _node_point(GI.getpoint(ring, i))
        _iig_add_seg!(mins, maxs, segs, prev, pt)
        prev = pt
    end
    if prev != first_pt
        _iig_add_seg!(mins, maxs, segs, prev, first_pt)
    end
    return nothing
end

function _iig_add_seg!(mins, maxs, segs, p0, p1)
    push!(mins, min(p0[2], p1[2]))
    push!(maxs, max(p0[2], p1[2]))
    push!(segs, (p0, p1))
    return nothing
end
