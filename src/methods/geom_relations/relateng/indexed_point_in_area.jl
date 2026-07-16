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

# One polygon's cached kernel rings (spherical): shell plus holes, in the
# even-odd composition order of `_locate_point_in_polygonal`.
struct _SphPolyRings
    shell::SphericalKernelRing
    holes::Vector{SphericalKernelRing}
end

# The spherical segment index (Layer 2 of the 2026-07-14
# spherical-indexed-locator design): ring edges as kernel-point pairs, on
# 1-D *longitude* intervals in radians — the spherical isomorph of the
# planar y-interval trick with `Y → X (lon)` and `ray to -∞ → meridian arc
# to a pole`. An edge can only meet the reference meridian arc if its
# longitude span contains the query longitude.
const _SphPIASegment = Tuple{UnitSphericalPoint{Float64}, UnitSphericalPoint{Float64}}
const _SphPIAExtent = Extents.Extent{(:X,), Tuple{NTuple{2, Float64}}}
const _SphPIAIndex = RTree{STR, _SphPIAExtent, Vector{_SphPIASegment}, Vector{Int}}

const _SPH_SOUTH_POLE = UnitSphericalPoint(0.0, 0.0, -1.0)
const _SPH_NORTH_POLE = UnitSphericalPoint(0.0, 0.0, 1.0)

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

On `Spherical`, the locator caches the element's rings in kernel space
(`polys` — [`SphericalKernelRing`](@ref)s; Layer 1 of the 2026-07-14
spherical-indexed-locator design), so queries never reconvert vertices, and
— with `indexed = true`, the default — builds the longitude-interval edge
index plus the parity anchor (Layer 2): the location of a reference pole,
computed once by the exact scan, from which each query is a 1-D stab at its
longitude and a crossing-parity count along its meridian arc to the pole.
With `indexed = false` (the unprepared arm, where a one-shot query cannot
amortize the build) it locates by the exact ring scan over the cached
rings. The `polys` vector is empty on `Planar`, where the corresponding
roles are played by `index` and the implicit EXTERIOR at the ray's far end.
"""
struct IndexedPointInAreaLocator{M <: Manifold, E}
    m::M
    exact::E
    index::Union{Nothing, _PIAIndex, _SphPIAIndex}
    polys::Vector{_SphPolyRings}
    #-- spherical parity anchor: the reference-arc far end (a pole whose
    #-- location the exact scan computed at build) and its location.
    #-- `anchor_loc == LOC_BOUNDARY` means both poles lie on the element
    #-- boundary: no index is built and every query takes the exact scan.
    anchor::UnitSphericalPoint{Float64}
    anchor_loc::Int8
end

function IndexedPointInAreaLocator(m::Planar, geom; exact)
    exts = _PIAExtent[]
    segs = _PIASegment[]
    n = GI.npoint(geom)
    sizehint!(exts, n); sizehint!(segs, n)
    _interval_index_add_geom!(exts, segs, GI.trait(geom), geom)
    index = isempty(segs) ? nothing : RTree(STR(), segs; extents = exts)
    return IndexedPointInAreaLocator(m, exact, index, _SphPolyRings[], _SPH_SOUTH_POLE, LOC_EXTERIOR)
end

function IndexedPointInAreaLocator(m::Spherical, geom; exact, indexed::Bool = true)
    polys = _SphPolyRings[]
    _sph_rings_add_geom!(polys, m, GI.trait(geom), geom; exact)
    anchor, anchor_loc = _SPH_SOUTH_POLE, LOC_EXTERIOR
    index = nothing
    if indexed && !isempty(polys)
        anchor, anchor_loc = _sph_pia_anchor(m, polys; exact)
        if anchor_loc != LOC_BOUNDARY
            index = _sph_lon_index(polys)
        end
    end
    return IndexedPointInAreaLocator(m, exact, index, polys, anchor, anchor_loc)
end

"""
    locate(loc::IndexedPointInAreaLocator, p)

The location (`LOC_*` code) of point `p` in the locator's areal geometry.
Port of `IndexedPointInAreaLocator.locate`.
"""
function locate(loc::IndexedPointInAreaLocator{<:Planar}, p)
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

function locate(loc::IndexedPointInAreaLocator{<:Spherical}, p)
    index = loc.index
    #-- signed-zero-normalize an already-kernel point: the index longitudes
    #-- come from normalized vertices, and atan(-0.0, x<0) is -π, not π
    q = p isa UnitSphericalPoint{Float64} ? _node_point(p) : _to_kernel_point(loc.m, p)
    #-- scan mode: unindexed (unprepared), empty, or both poles on boundary
    index isa _SphPIAIndex ||
        return _sph_scan_locate(loc.m, loc.polys, q; exact = loc.exact)
    if q[1] == 0.0 && q[2] == 0.0
        #-- on the polar axis: the anchor itself answers from the stored
        #-- bit; its antipode has no reference arc (antipodal pair), so it
        #-- takes the exact scan
        q == loc.anchor && return loc.anchor_loc
        return _sph_scan_locate(loc.m, loc.polys, q; exact = loc.exact)
    end
    acc = ArcCrossingCounter(loc.m, loc.exact, q, loc.anchor, 0, false)
    λ = atan(q[2], q[1])
    stab = Extents.Extent(X = (λ, λ))
    #-- count every edge whose longitude interval contains the query's
    SpatialTreeInterface.depth_first_search(Base.Fix1(Extents.intersects, stab), index) do i
        seg = index.data[i]
        count_arc_segment!(acc, seg[1], seg[2])
    end
    acc.is_point_on_segment && return LOC_BOUNDARY
    if isodd(acc.crossing_count)
        return loc.anchor_loc == LOC_INTERIOR ? LOC_EXTERIOR : LOC_INTERIOR
    end
    return loc.anchor_loc
end

#==========================================================================
## Spherical kernel-ring cache (Layer 1 of the spherical indexed locator)
==========================================================================#

# The exact ring scan of `_locate_point_in_polygonal` (point_locator.jl)
# over cached kernel rings: same shell-then-holes even-odd composition,
# with the conversion, dedup, and orientation bit precomputed.
function _sph_scan_locate(m::Spherical, polys::Vector{_SphPolyRings}, p; exact)
    for pr in polys
        l = _sph_locate_in_poly(m, pr, p; exact)
        l != LOC_EXTERIOR && return l
    end
    return LOC_EXTERIOR
end

function _sph_locate_in_poly(m::Spherical, pr::_SphPolyRings, p; exact)
    shell_loc = rk_point_in_ring(m, p, pr.shell; exact)
    shell_loc != LOC_INTERIOR && return shell_loc
    for hole in pr.holes
        hole_loc = rk_point_in_ring(m, p, hole; exact)
        hole_loc == LOC_BOUNDARY && return LOC_BOUNDARY
        hole_loc == LOC_INTERIOR && return LOC_EXTERIOR
        #-- if in EXTERIOR of this hole keep checking the other ones
    end
    return LOC_INTERIOR
end

function _sph_rings_add_geom!(polys, m::Spherical, ::GI.PolygonTrait, poly; exact)
    GI.isempty(poly) && return nothing
    shell = SphericalKernelRing(m, GI.getexterior(poly); exact)
    holes = [SphericalKernelRing(m, h; exact, is_hole = true)
             for h in GI.gethole(poly) if !GI.isempty(h)]
    push!(polys, _SphPolyRings(shell, holes))
    return nothing
end

function _sph_rings_add_geom!(polys, m::Spherical, ::GI.MultiPolygonTrait, mp; exact)
    for poly in GI.getgeom(mp)
        _sph_rings_add_geom!(polys, m, GI.trait(poly), poly; exact)
    end
    return nothing
end

#==========================================================================
## Spherical longitude-interval index (Layer 2)
==========================================================================#

#=
Parity anchor: in the plane the far end of the crossing ray is at infinity
and EXTERIOR by construction; on the sphere the far end is a pole, whose
location relative to the element is computed once by the exact scan. The
south pole is preferred; if it lies ON the boundary the reference arc
could graze the linework everywhere along it, so the north pole takes
over. Both poles on the boundary leaves `LOC_BOUNDARY`, which the
constructor maps to unindexed scan mode — pathological, and correctness
is cheap.
=#
function _sph_pia_anchor(m::Spherical, polys; exact)
    loc_s = _sph_scan_locate(m, polys, _SPH_SOUTH_POLE; exact)
    loc_s != LOC_BOUNDARY && return (_SPH_SOUTH_POLE, loc_s)
    loc_n = _sph_scan_locate(m, polys, _SPH_NORTH_POLE; exact)
    return (_SPH_NORTH_POLE, loc_n)
end

"""
    ArcCrossingCounter(m::Spherical, exact, p, anchor, 0, false)

Counts ring edges crossing the reference meridian arc from the query point
`p` to the `anchor` pole, in an incremental fashion — the spherical
counterpart of [`RayCrossingCounter`](@ref). As there, the location is only
correct once **all** edges whose longitude interval contains `p`'s
longitude have been counted, and a query point found to lie on an edge is
recorded in `is_point_on_segment` (final location `LOC_BOUNDARY`).

The final location is `anchor`'s location, flipped once per crossing.
Vertex grazing — an edge endpoint exactly on the reference arc — is
resolved symbolically, S2-`VertexCrossing` style: the edge counts iff its
off-arc endpoint lies strictly on the positive side of the meridian great
circle, so a crossing pair of incident edges counts once, a same-side pair
counts zero or twice (parity-equal), and edges collinear with the meridian
count never (their chain's terminal edges decide). With `exact = True()`
every branch below is decided by exact predicates.
"""
mutable struct ArcCrossingCounter{M <: Spherical, E}
    const m::M
    const exact::E
    const p::UnitSphericalPoint{Float64}
    const anchor::UnitSphericalPoint{Float64}
    crossing_count::Int
    is_point_on_segment::Bool
end

function count_arc_segment!(acc::ArcCrossingCounter, a, b)
    m, exact = acc.m, acc.exact
    q, s = acc.p, acc.anchor
    #-- the query point on the closed edge (vertex or interior, including
    #-- collinear overlap): boundary
    if rk_point_on_segment(m, q, a, b; exact)
        acc.is_point_on_segment = true
        return nothing
    end
    a == b && return nothing   # zero-length edge: only boundary-relevant
    sa = rk_orient(m, q, s, a; exact)
    sb = rk_orient(m, q, s, b; exact)
    if sa == 0 || sb == 0
        #-- edge collinear with the meridian circle: were q interior to it
        #-- the boundary test above caught it; the anchor is never on an
        #-- edge (anchor selection); otherwise its neighbors decide
        (sa == 0 && sb == 0) && return nothing
        #-- vertex grazing: an endpoint on the meridian great circle. Two
        #-- distinct great circles meet only at one antipodal point pair, so
        #-- the edge can touch the reference arc only at that endpoint —
        #-- count iff the endpoint is ON the arc and the other endpoint is
        #-- strictly on the positive side (see the docstring)
        von, s_off = sa == 0 ? (a, sb) : (b, sa)
        if s_off > 0 && rk_point_on_segment(m, von, q, s; exact)
            acc.crossing_count += 1
        end
        return nothing
    end
    #-- endpoints strictly on the same side: no crossing
    (sa > 0) == (sb > 0) && return nothing
    #-- q (or the anchor) on the edge's great circle but not on the edge:
    #-- the circles meet only at ±q (resp. ±anchor), out of reach of both
    #-- arcs — cf. `_arc_crossing_parity`
    sq = rk_orient(m, a, b, q; exact)
    sq == 0 && return nothing
    sm = rk_orient(m, a, b, s; exact)
    sm == 0 && return nothing
    (sq > 0) == (sm > 0) && return nothing
    #=
    Mutual strict straddle: the two great circles meet at ±x, x = n₁×n₂,
    and each arc contains exactly one of the pair — they cross iff it is
    the same one. Which one each arc contains is already encoded in the
    computed signs: traveling the reference arc q → s (tangent n₁×p), the
    plane-2 crossing is +x iff sq > 0; traveling the edge a → b (tangent
    n₂×p), the plane-1 crossing is −x iff sa > 0. So the arcs cross iff
    sa and sq differ — no constructed point, and exactly as exact as
    `rk_orient` (this replaces a `Rational{BigInt}` in-arc confirmation
    that dominated indexed query time).
    =#
    if (sa > 0) != (sq > 0)
        acc.crossing_count += 1
    end
    return nothing
end

# The longitude-interval edge index over every ring of the element, built
# from the same cached kernel points the exact scan walks.
function _sph_lon_index(polys)
    exts = _SphPIAExtent[]
    segs = _SphPIASegment[]
    for pr in polys
        _sph_lon_index_add_ring!(exts, segs, pr.shell)
        for hole in pr.holes
            _sph_lon_index_add_ring!(exts, segs, hole)
        end
    end
    return isempty(segs) ? nothing : RTree(STR(), segs; extents = exts)
end

# Index the boundary edge walk (consecutive vertex pairs plus the implicit
# closing edge) — the same edge set `rk_point_in_ring` tests for boundary
# and parity.
function _sph_lon_index_add_ring!(exts, segs, kr::SphericalKernelRing)
    pts = kr.pts
    n = length(pts)
    n < 2 && return nothing
    for i in 1:(n - 1)
        _sph_lon_entries!(exts, segs, pts[i], pts[i + 1])
    end
    if pts[n] != pts[1]
        _sph_lon_entries!(exts, segs, pts[n], pts[1])
    end
    return nothing
end

#=
The longitude interval(s) of one great-circle edge. Longitude is strictly
monotonic along a great-circle arc — the east component of the tangent at
`p` is `(n × p) ⋅ (ẑ × p) = n_z`, constant along the whole circle (`n` the
circle normal) — so an edge spans exactly the wrapped interval between its
endpoint longitudes; and since `n_z = cos(lat_a) cos(lat_b) sin(λ_b − λ_a)`
the sweep is always the SHORTER of the two candidate intervals. (The design
note's worry that longitude "bulges" like latitude does not arise.)
Intervals are padded a few ulps against `atan` rounding, as the arc extents
pad. Conservative full-interval fallbacks, where the float longitudes do
not determine the sweep: an endpoint exactly on the polar axis (undefined
longitude), or endpoint longitudes within ~1e-9 of half a turn apart
(pole-hugging edges — the short/long choice would hang on the sign of a
vanishing cross product). An interval crossing the antimeridian contributes
two entries; the two never overlap, so no edge is double-counted.
=#
function _sph_lon_entries!(exts, segs, a, b)
    seg = (a, b)
    halfturn = Float64(π)
    if (a[1] == 0.0 && a[2] == 0.0) || (b[1] == 0.0 && b[2] == 0.0)
        return _push_lon_entry!(exts, segs, seg, -halfturn, halfturn)
    end
    λa = atan(a[2], a[1])
    λb = atan(b[2], b[1])
    Δ = λb - λa
    Δ > halfturn && (Δ -= 2halfturn)
    Δ < -halfturn && (Δ += 2halfturn)
    if abs(Δ) > halfturn - 1e-9
        return _push_lon_entry!(exts, segs, seg, -halfturn, halfturn)
    end
    lo, hi = Δ >= 0 ? (λa, λa + Δ) : (λa + Δ, λa)
    lo = prevfloat(lo, 32)
    hi = nextfloat(hi, 32)
    #-- antimeridian wraparound: split the overflow back into [-π, π]
    if lo < -halfturn
        _push_lon_entry!(exts, segs, seg, lo + 2halfturn, halfturn)
        lo = -halfturn
    elseif hi > halfturn
        _push_lon_entry!(exts, segs, seg, -halfturn, hi - 2halfturn)
        hi = halfturn
    end
    return _push_lon_entry!(exts, segs, seg, lo, hi)
end

function _push_lon_entry!(exts, segs, seg, lo, hi)
    push!(exts, Extents.Extent(X = (lo, hi)))
    push!(segs, seg)
    return nothing
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
