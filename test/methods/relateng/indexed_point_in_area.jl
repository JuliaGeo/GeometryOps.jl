# Tests for the prepared-mode indexed point-in-area locator
# (indexed_point_in_area.jl): the 1-D y-interval segment index, the
# RayCrossingCounter / IndexedPointInAreaLocator ports, and prepared- vs
# unprepared-mode agreement of RelatePointLocator point location. The
# unprepared SimplePointInAreaLocator ring loop is the oracle: prepared mode
# must locate every point identically. On-edge points are constructed on
# horizontal/vertical edges (exactly representable) and query points
# deliberately share y-coordinates with ring vertices — the classic
# RayCrossingCounter edge cases (vertex on ray, horizontal edge on ray).
# The spherical analogue (longitude-interval index + meridian-arc crossing
# parity against a pole anchor) is tested the same way at the bottom, with
# the unindexed exact ring scan as the oracle.

using Test
using Random
using LinearAlgebra: normalize
import GeometryOps as GO
import GeometryOps: Planar, Spherical, True
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
    loc = GO.IndexedPointInAreaLocator(Planar(), True(), nothing, GO._SphPolyRings[],
        GO._SPH_SOUTH_POLE, GO.LOC_EXTERIOR)
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

# -- spherical indexed locator -------------------------------------------------

# The unindexed exact ring scan over the cached kernel rings (`indexed =
# false`, the unprepared arm) is the oracle: the longitude-interval stab plus
# meridian-arc crossing parity must locate every point identically. Clouds
# are seeded, and every ring runs in both windings — the enclosed region is
# winding-independent and the orientation bit lives in the ring cache.
@testset "spherical indexed locator" begin
    m = Spherical()
    kp(ll) = GO._to_kernel_point(m, GI.Point(ll))
    function check_sph_agreement(geom, lls; indexed_built = true)
        scan = GO.IndexedPointInAreaLocator(m, geom; exact = True(), indexed = false)
        idx = GO.IndexedPointInAreaLocator(m, geom; exact = True())
        @test scan.index === nothing
        @test (idx.index isa GO._SphPIAIndex) == indexed_built
        n_mismatch = count(ll -> GO.locate(idx, kp(ll)) != GO.locate(scan, kp(ll)), lls)
        @test n_mismatch == 0
    end
    rng = Xoshiro(11)
    cloud(n, lonr, latr) = [(lonr[1] + rand(rng) * (lonr[2] - lonr[1]),
                             latr[1] + rand(rng) * (latr[2] - latr[1])) for _ in 1:n]

    @testset "mid-latitude star, both windings" begin
        star = [(20.0 + (5 + 2rand(rng)) * cosd(t), 40.0 + (5 + 2rand(rng)) * sind(t))
                for t in 0:15:345]
        push!(star, star[1])
        for pts in (star, reverse(star))
            check_sph_agreement(GI.Polygon([GI.LinearRing(pts)]),
                vcat(cloud(400, (10, 30), (30, 50)), pts))
        end
    end

    @testset "antimeridian-crossing box, both windings" begin
        # edges straddle lon ±180°, so their intervals split into two index
        # entries; queries sit on both sides plus the seam itself
        am = [(170.0, -10.0), (-170.0, -10.0), (-170.0, 10.0), (170.0, 10.0), (170.0, -10.0)]
        for pts in (am, reverse(am))
            qs = vcat(cloud(200, (160, 180), (-20, 20)), cloud(200, (-180, -160), (-20, 20)), pts,
                      [(180.0, 0.0), (-180.0, 5.0), (175.0, 0.0), (-175.0, 0.0), (165.0, 0.0), (-165.0, 0.0)])
            check_sph_agreement(GI.Polygon([GI.LinearRing(pts)]), qs)
        end
    end

    @testset "polar caps, both poles, both windings" begin
        # a cap encloses its pole; the pole query exercises the polar-axis
        # fallback (no meridian reference arc from a pole), and for the south
        # cap the anchor itself is INTERIOR
        ncap = [(t, 80.0) for t in 0.0:30.0:330.0]; push!(ncap, ncap[1])
        scap = [(t, -80.0) for t in 0.0:30.0:330.0]; push!(scap, scap[1])
        for pts in (ncap, scap), w in (pts, reverse(pts))
            qs = vcat(cloud(200, (-180, 180), (60, 90)), cloud(100, (-180, 180), (-90, -60)), w,
                      [(0.0, 90.0), (0.0, -90.0), (45.0, 85.0), (45.0, -85.0)])
            check_sph_agreement(GI.Polygon([GI.LinearRing(w)]), qs)
        end
    end

    @testset "south pole on the boundary: anchor falls back to the north pole" begin
        spb = [(0.0, -90.0), (20.0, -60.0), (-20.0, -60.0), (0.0, -90.0)]
        p = GI.Polygon([GI.LinearRing(spb)])
        idx = GO.IndexedPointInAreaLocator(m, p; exact = True())
        @test idx.anchor == GO._SPH_NORTH_POLE
        @test idx.anchor_loc == GO.LOC_EXTERIOR
        @test GO.locate(idx, kp((0.0, -90.0))) == GO.LOC_BOUNDARY   # the pole itself
        check_sph_agreement(p, vcat(cloud(300, (-30, 30), (-90, -50)), spb, [(0.0, -90.0)]))
    end

    @testset "both poles on the boundary: unindexed scan mode" begin
        bpb = [(0.0, -90.0), (20.0, 0.0), (0.0, 90.0), (40.0, 0.0), (0.0, -90.0)]
        p = GI.Polygon([GI.LinearRing(bpb)])
        idx = GO.IndexedPointInAreaLocator(m, p; exact = True())
        @test idx.index === nothing            # no anchor: every query scans
        @test idx.anchor_loc == GO.LOC_BOUNDARY
        check_sph_agreement(p, vcat(cloud(300, (-10, 50), (-80, 80)), bpb);
            indexed_built = false)
    end

    @testset "boundary points locate as LOC_BOUNDARY" begin
        eqp = [(0.0, 0.0), (30.0, 0.0), (15.0, 20.0), (0.0, 0.0)]
        idx = GO.IndexedPointInAreaLocator(m, GI.Polygon([GI.LinearRing(eqp)]); exact = True())
        for lon in (5.0, 10.0, 15.0, 29.0)      # on the equator edge
            @test GO.locate(idx, kp((lon, 0.0))) == GO.LOC_BOUNDARY
        end
        @test GO.locate(idx, kp((0.0, 0.0))) == GO.LOC_BOUNDARY     # vertex
        @test GO.locate(idx, kp((15.0, 5.0))) == GO.LOC_INTERIOR
        @test GO.locate(idx, kp((15.0, -5.0))) == GO.LOC_EXTERIOR
    end

    @testset "polygon with a hole" begin
        hp = GI.Polygon([GI.LinearRing([(0., 0.), (20., 0.), (20., 20.), (0., 20.), (0., 0.)]),
                         GI.LinearRing([(5., 5.), (15., 5.), (15., 15.), (5., 15.), (5., 5.)])])
        check_sph_agreement(hp, vcat(cloud(400, (-2, 22), (-2, 22)),
            [(10.0, 10.0), (10.0, 5.0), (2.0, 2.0)]))
    end

    @testset "edge index: endpoint longitude rounding" begin
        # a vertex at a tiny non-zero longitude next to a long edge: an
        # interval end that is not the endpoint's own `atan` value can miss
        # it by ~ulp(edge span), thousands of ulps of the tiny end, so every
        # query on that side of the vertex's meridian lost the edge and
        # inverted its parity — an ordinary mid-latitude ring, no pole
        tri = [(0.0017, 30.0), (40.0, 30.0), (20.0, 60.0), (0.0017, 30.0)]
        v = kp(tri[1]); λ = atan(v[2], v[1])
        qs = [(rad2deg(k >= 0 ? nextfloat(λ, k) : prevfloat(λ, -k)), lat)
              for lat in 31.0:2.0:59.0 for k in -40:40]
        append!(qs, [(0.0017, lat) for lat in 31.0:1.0:59.0])
        check_sph_agreement(GI.Polygon([GI.LinearRing(tri)]), qs)
        # the same shape at a pole-hugging vertex, where the longitude
        # itself is ill-conditioned and the edge takes the full range
        ph = [(0.0033563130853053735, -89.99999999813261), (46.891183665388354, -89.9999999962315),
              (93.63815566658627, -89.99999999842129), (142.0182497462225, -89.99999993297193),
              (180.64483654387297, -72.73095921005729), (225.34348517929183, -89.9999999999999),
              (277.31417119248607, -65.32296008285503), (319.1523906912383, -54.67823057378641),
              (0.0033563130853053735, -89.99999999813261)]
        v = kp(ph[1]); λ = atan(v[2], v[1])
        qs = [(rad2deg(k >= 0 ? nextfloat(λ, k) : prevfloat(λ, -k)), -lat)
              for lat in (30.0, 60.0, 85.0, 89.9, 89.999999) for k in (-200, -64, -33, -1, 0, 1, 33, 64, 200)]
        check_sph_agreement(GI.Polygon([GI.LinearRing(ph)]), qs)
    end

    @testset "pole-hugging rings: seeded fuzz against the exact scan" begin
        # vertices within 1e-6..1e-12 degrees of a pole mixed with ordinary
        # ones, queries random plus stabs a few-to-many ulps around each
        # vertex longitude — the index must never drop or double-count an edge
        frng = Xoshiro(3)
        for trial in 1:60
            pole = rand(frng, (-1, 1))
            pts = Tuple{Float64, Float64}[]
            for t in 0:45:315
                lat = rand(frng, Bool) ? pole * (90 - 10.0^(-rand(frng, 6:12)) * rand(frng)) :
                                         pole * (50 + 30rand(frng))
                push!(pts, (t + 10rand(frng), lat))
            end
            push!(pts, pts[1])
            geom = GI.Polygon([GI.LinearRing(pts)])
            qs = [(360rand(frng) - 180, pole * (40 + 50rand(frng))) for _ in 1:100]
            for ll in pts
                v = kp(ll); λ = atan(v[2], v[1])
                for k in (-200, -33, 0, 33, 200), lat in (30.0, 85.0, 89.999999)
                    push!(qs, (rad2deg(k >= 0 ? nextfloat(λ, k) : prevfloat(λ, -k)), pole * lat))
                end
            end
            check_sph_agreement(geom, qs)
        end
    end

    @testset "ordinary mid-latitude rings: seeded fuzz against the exact scan" begin
        frng = Xoshiro(7)
        for trial in 1:20
            c = (360rand(frng) - 180, 120rand(frng) - 60)
            pts = [(c[1] + (5 + 3rand(frng)) * cosd(t), c[2] + (5 + 3rand(frng)) * sind(t)) for t in 0:30:330]
            push!(pts, pts[1])
            qs = [(c[1] + 20rand(frng) - 10, c[2] + 20rand(frng) - 10) for _ in 1:150]
            append!(qs, [(ll[1], c[2] + 12rand(frng) - 6) for ll in pts])   # share vertex longitudes
            check_sph_agreement(GI.Polygon([GI.LinearRing(pts)]), qs)
        end
    end

    #= An equatorial band whose meridian ends carry far more vertices than its
    parallels (as a buffer's round caps do). Vertices spanning more than a
    hemisphere admit no exterior anchor, so both locators must fall back to the
    winding bootstrap; membership is analytic: `0 < lon < L`, `|lat| < w`. =#
    function band(L; w = 20.0, step_along = 10.0, step_end = 0.25)
        pts = Tuple{Float64, Float64}[]
        for lon in 0.0:step_along:L; push!(pts, (lon, w)); end
        for lat in w-step_end:-step_end:-w+step_end; push!(pts, (L, lat)); end
        for lon in L:-step_along:0.0; push!(pts, (lon, -w)); end
        for lat in -w+step_end:step_end:w-step_end; push!(pts, (0.0, lat)); end
        push!(pts, pts[1])
        return GI.Polygon([GI.LinearRing(pts)])
    end
    function check_band(L; w = 20.0, N = 300, clearance = 3.0)
        poly = band(L; w)
        idx = GO.IndexedPointInAreaLocator(m, poly; exact = True())
        scan = GO.IndexedPointInAreaLocator(m, poly; exact = True(), indexed = false)
        brng = Xoshiro(7)
        n_wrong_idx = n_wrong_scan = 0
        n = 0
        while n < N
            lon = rand(brng) * 360 - 180
            lat = rand(brng) * 180 - 90
            # well clear of the band's parallels and meridian ends
            abs(abs(lat) - w) < clearance && continue
            min(mod(lon, 360), abs(mod(lon, 360) - L)) < clearance && continue
            n += 1
            truth = (0 < mod(lon, 360) < L) && abs(lat) < w
            q = kp((lon, lat))
            n_wrong_idx += (GO.locate(idx, q) == GO.LOC_INTERIOR) != truth
            n_wrong_scan += (GO.locate(scan, q) == GO.LOC_INTERIOR) != truth
        end
        @test n_wrong_idx == 0
        @test n_wrong_scan == 0
    end

    @testset "equatorial band spanning 240° of longitude" begin
        check_band(240.0)
    end

    @testset "equatorial band spanning 200° of longitude" begin
        check_band(200.0)
    end

    @testset "small ordinary ring: the vertex-mass anchor still answers" begin
        v = [GO.UnitSpherical.UnitSphereFromGeographic()(p)
             for p in GI.getpoint(GI.getexterior(band(60.0)))][1:end-1]
        mass = sum(normalize(p) for p in v)
        @test GO.UnitSpherical.spherical_exterior_anchor(v, length(v)) ==
              GO.UnitSpherical.UnitSphericalPoint(-normalize(mass))
        check_band(60.0)
    end

    @testset "empty spherical element" begin
        # as in the planar empty case: no rings, no index, everything exterior
        loc = GO.IndexedPointInAreaLocator(m, True(), nothing, GO._SphPolyRings[],
            GO._SPH_SOUTH_POLE, GO.LOC_EXTERIOR)
        @test GO.locate(loc, kp((0.0, 0.0))) == GO.LOC_EXTERIOR
        @test GO.locate(loc, kp((0.0, 90.0))) == GO.LOC_EXTERIOR
    end
end
