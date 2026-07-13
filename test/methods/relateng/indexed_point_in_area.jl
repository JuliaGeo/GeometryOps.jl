# Tests for the prepared-mode indexed point-in-area locator
# (indexed_point_in_area.jl): the 1-D y-interval segment index, the
# RayCrossingCounter / IndexedPointInAreaLocator ports, and prepared- vs
# unprepared-mode agreement of RelatePointLocator point location. The
# unprepared SimplePointInAreaLocator ring loop is the oracle: prepared mode
# must locate every point identically. On-edge points are constructed on
# horizontal/vertical edges (exactly representable) and query points
# deliberately share y-coordinates with ring vertices — the classic
# RayCrossingCounter edge cases (vertex on ray, horizontal edge on ray).

using Test
import GeometryOps as GO
import GeometryOps: Planar, True
import GeoInterface as GI
import Extents

@testset "1-D y-interval stabbing" begin
    # the interval-index shape the locator builds: RTree(STR(), items;
    # extents = y-intervals), queried with a closed [qmin, qmax] extent
    interval_tree(mins, maxs, items) =
        GO.FlexibleRTrees.RTree(GO.FlexibleRTrees.STR(), items;
            extents = [Extents.Extent(Y = (mins[i], maxs[i])) for i in eachindex(mins)])
    collect_query(tree, qmin, qmax) = begin
        out = Int[]
        q = Extents.Extent(Y = (qmin, qmax))
        GO.SpatialTreeInterface.depth_first_search(Base.Fix1(Extents.intersects, q), tree) do i
            push!(out, tree.data[i])
        end
        sort!(out)
    end

    # single item
    one = interval_tree([1.0], [2.0], [1])
    @test collect_query(one, 1.5, 1.5) == [1]
    @test collect_query(one, 2.5, 3.0) == Int[]
    @test collect_query(one, 2.0, 3.0) == [1]   # closed-interval touch

    # several overlapping intervals, including duplicates and a point interval
    mins = [0.0, 1.0, 2.0, 2.0, 5.0, 5.0, -3.0]
    maxs = [1.0, 3.0, 4.0, 2.0, 9.0, 6.0, -1.0]
    items = collect(1:7)
    tree = interval_tree(mins, maxs, items)
    brute(qmin, qmax) = sort!([i for i in 1:7 if !(mins[i] > qmax || maxs[i] < qmin)])
    for (qmin, qmax) in [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (2.5, 2.5),
                         (-4.0, -3.5), (-2.0, 0.5), (3.5, 5.0), (10.0, 11.0),
                         (-10.0, 10.0), (4.5, 4.9)]
        @test collect_query(tree, qmin, qmax) == brute(qmin, qmax)
    end
end

@testset "RayCrossingCounter" begin
    m = Planar()
    square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]
    locate_in_ring(p, ring) = begin
        # port of RayCrossingCounter.locatePointInRing as the usage exemplar
        rcc = GO.RayCrossingCounter(m, p; exact = True())
        for i in 2:length(ring)
            GO.count_segment!(rcc, ring[i], ring[i - 1])
            GO.is_on_segment(rcc) && return GO.rcc_location(rcc)
        end
        return GO.rcc_location(rcc)
    end
    @test locate_in_ring((5.0, 5.0), square) == GO.LOC_INTERIOR
    @test locate_in_ring((15.0, 5.0), square) == GO.LOC_EXTERIOR
    @test locate_in_ring((-5.0, 5.0), square) == GO.LOC_EXTERIOR
    @test locate_in_ring((0.0, 0.0), square) == GO.LOC_BOUNDARY    # vertex
    @test locate_in_ring((5.0, 0.0), square) == GO.LOC_BOUNDARY    # horizontal edge
    @test locate_in_ring((10.0, 5.0), square) == GO.LOC_BOUNDARY   # vertical edge
    # ray passes exactly through vertices: a diamond, query at vertex height
    diamond = [(0.0, 0.0), (5.0, -5.0), (10.0, 0.0), (5.0, 5.0), (0.0, 0.0)]
    @test locate_in_ring((5.0, 0.0), diamond) == GO.LOC_INTERIOR   # ray exits through vertex (10,0)
    @test locate_in_ring((-1.0, 0.0), diamond) == GO.LOC_EXTERIOR  # ray enters AND exits through vertices
    @test locate_in_ring((11.0, 0.0), diamond) == GO.LOC_EXTERIOR
    @test locate_in_ring((0.0, 0.0), diamond) == GO.LOC_BOUNDARY
end

# -- prepared vs unprepared location agreement --------------------------------

# 10k-vertex circle
const N_CIRC = 10_000
circ = [(cos(t), sin(t)) for t in range(0.0, 2pi; length = N_CIRC)]
circ[end] = circ[1]
poly_circle = GI.Polygon([circ])

# polygon with two holes; shell/hole1 axis-aligned so on-edge points are exact
shell = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]
hole1 = [(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0), (2.0, 2.0)]
hole2 = [(6.0, 6.0), (8.0, 6.0), (7.0, 8.0), (6.0, 6.0)]
poly_holes = GI.Polygon([shell, hole1, hole2])

# multipolygon: two squares, second with a hole
mp = GI.MultiPolygon([
    [[(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (0.0, 5.0), (0.0, 0.0)]],
    [[(10.0, 0.0), (15.0, 0.0), (15.0, 5.0), (10.0, 5.0), (10.0, 0.0)],
     [(11.0, 1.0), (14.0, 1.0), (14.0, 4.0), (11.0, 4.0), (11.0, 1.0)]],
])

function check_prepared_agreement(geom, pts)
    m = Planar()
    loc_prep = GO.RelatePointLocator(m, geom; exact = True(), is_prepared = true)
    loc_unprep = GO.RelatePointLocator(m, geom; exact = True(), is_prepared = false)
    n_mismatch = 0
    for pt in pts
        # unprepared = direct ring loop, prepared = indexed locator — a true
        # indexed-vs-simple differential on every query
        GO.locate(loc_prep, pt) == GO.locate(loc_unprep, pt) || (n_mismatch += 1)
        GO.locate_with_dim(loc_prep, pt) == GO.locate_with_dim(loc_unprep, pt) || (n_mismatch += 1)
    end
    @test n_mismatch == 0
    # cache sanity: prepared mode built (at most) one locator per polygonal
    # element; unprepared mode never builds one
    @test length(loc_prep.poly_locator) == length(loc_prep.polygons)
    @test all(isnothing, loc_unprep.poly_locator)
end

function check_prepared_relate(geom, pts)
    alg = GO.RelateNG()
    prep = GO.prepare(alg, geom)
    n_mismatch = 0
    for pt in pts
        p = GI.Point(pt)
        GO.relate(prep, p) == GO.relate(alg, geom, p) || (n_mismatch += 1)
    end
    @test n_mismatch == 0
end

@testset "prepared vs unprepared: 10k circle" begin
    pts = Vector{Tuple{Float64, Float64}}()
    append!(pts, circ[1:97:end])                                       # exact vertices
    append!(pts, [(0.999c[1], 0.999c[2]) for c in circ[1:301:end]])    # just inside
    append!(pts, [(1.001c[1], 1.001c[2]) for c in circ[1:301:end]])    # just outside
    append!(pts, [(x, c[2]) for c in circ[1:211:end] for x in (-2.0, 0.0, 0.5, 2.0)])  # share y with vertices
    push!(pts, (0.0, 0.0))
    push!(pts, (1.0, 0.0))                                             # the t = 0 vertex
    push!(pts, (0.0, -1.5))
    check_prepared_agreement(poly_circle, pts)
    check_prepared_relate(poly_circle, pts)
end

@testset "prepared vs unprepared: polygon with holes" begin
    pts = Vector{Tuple{Float64, Float64}}()
    append!(pts, shell); append!(pts, hole1); append!(pts, hole2)      # exact vertices
    append!(pts, [(5.0, 0.0), (10.0, 5.0), (5.0, 10.0), (0.0, 5.0)])   # on shell edges (horiz + vert)
    append!(pts, [(3.0, 2.0), (4.0, 3.0), (3.0, 4.0), (2.0, 3.0)])     # on hole1 edges
    append!(pts, [(7.0, 6.0), (6.5, 7.0), (7.5, 7.0)])                 # on hole2 edges
    # rays through vertices and along horizontal edges, from inside,
    # in-hole, and outside positions
    append!(pts, [(1.0, 2.0), (3.0, 2.0), (5.0, 2.0), (1.0, 4.0), (5.0, 4.0),
                  (0.5, 0.0), (-1.0, 0.0), (-1.0, 2.0), (-1.0, 10.0), (5.0, 6.0),
                  (-1.0, 6.0), (9.0, 8.0)])
    append!(pts, [(3.0, 3.0), (7.0, 6.5), (1.0, 1.0), (5.0, 5.0)])     # hole interiors + interior
    append!(pts, [(-1.0, 5.0), (11.0, 5.0), (5.0, -1.0), (5.0, 11.0)]) # exterior
    # dense grid: hits vertices, edges, and every shared-y configuration
    append!(pts, vec([(x, y) for x in -1.0:0.5:11.0, y in -1.0:0.5:11.0]))
    check_prepared_agreement(poly_holes, pts)
    check_prepared_relate(poly_holes, pts)
end

@testset "prepared vs unprepared: multipolygon" begin
    pts = Vector{Tuple{Float64, Float64}}()
    append!(pts, vec([(x, y) for x in -1.0:0.5:16.0, y in -1.0:0.5:6.0]))
    append!(pts, [(2.5, 0.0), (12.5, 0.0), (12.5, 1.0), (12.5, 4.0),   # on edges
                  (11.0, 2.0), (14.0, 2.0), (7.5, 2.5), (12.5, 2.5)])  # in gap / in hole
    check_prepared_agreement(mp, pts)
    check_prepared_relate(mp, pts)
end

@testset "empty polygonal element" begin
    # the GI.Polygon wrapper cannot represent POLYGON EMPTY (zero rings), so
    # exercise the no-segments short-circuit on a directly constructed locator
    loc = GO.IndexedPointInAreaLocator(Planar(), True(), nothing)
    @test loc.index === nothing
    @test GO.locate(loc, (0.0, 0.0)) == GO.LOC_EXTERIOR
end

@testset "implicitly closed ring" begin
    # no repeated closing point: the indexed locator must close the ring,
    # matching rk_point_in_ring's assumed-closed semantics
    open_tri = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (0.0, 10.0)]])
    loc = GO.IndexedPointInAreaLocator(Planar(), open_tri; exact = True())
    @test GO.locate(loc, (1.0, 1.0)) == GO.LOC_INTERIOR
    @test GO.locate(loc, (5.0, 5.0)) == GO.LOC_BOUNDARY   # on the implicit closing edge
    @test GO.locate(loc, (6.0, 6.0)) == GO.LOC_EXTERIOR
    @test GO.locate(loc, (5.0, 0.0)) == GO.LOC_BOUNDARY
end
