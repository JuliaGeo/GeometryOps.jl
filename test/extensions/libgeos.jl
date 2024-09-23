using Test
import GeometryOps as GO, GeoInterface as GI, LibGEOS as LG
using GeometryOps: GEOS

pt1 = GI.Point((0.0, 0.0))
mpt1 = GI.MultiPoint([pt1, pt1])
l1 = GI.Line([(0.0, 0.0), (0.0, 1.0)])

concave_coords = [(0.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 2.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)]
l2 = GI.LineString(concave_coords)
l3 = GI.LineString(concave_coords[1:(end - 1)])
r1 = GI.LinearRing(concave_coords)
r2 = GI.LinearRing(concave_coords[1:(end - 1)])
r3 = GI.LinearRing([(1.0, 1.0), (1.0, 1.5), (1.5, 1.5), (1.5, 1.0), (1.0, 1.0)])
concave_angles = [90.0, 270.0, 90.0, 90.0, 90.0, 90.0]

p1 = GI.Polygon([r3])
p2 = GI.Polygon([[(0.0, 0.0), (0.0, 4.0), (3.0, 0.0), (0.0, 0.0)]])
p3 = GI.Polygon([[(-3.0, -2.0), (0.0,0.0), (5.0, 0.0), (-3.0, -2.0)]])
p4 = GI.Polygon([r1])
p5 = GI.Polygon([r1, r3])

mp1 = GI.MultiPolygon([p2, p3])
c1 = GI.GeometryCollection([pt1, l2, p2])

@testset "GeometryOpsLibGEOSExt with GeoInterface Geometries" begin
    @testset "Functionality Tests" begin
        @testset "Buffer" begin
            @test GO.equals(GO.buffer(p1, 1.0), LG.buffer(p1, 1.0))
            @test_nowarn GO.buffer(GO.GEOS(), [p1 p2], 1.0) # test non-geometry inputs - if it works for this, it works for any `apply`.
            @test_nowarn GO.buffer(GO.GEOS(), c1, 1.0)      # test geometrycollection specifically - that's the reason for the awful hack in `GEOSExt/buffers.jl`.
        end
        # DE-9IM functions are tested in methods/geom_relations.jl, through the GO.GEOS wrapper.
        # Now, we test polygon intersection functions:
        @testset "Polygon clipping" begin
            # These tests would be nice but probably don't work, plus we'd have to munge LibGEOS output.
            #=
            @test GO.intersection(p1, p2) == GO.intersection(GO.GEOS(), p1, p2)
            @test GO.union(p1, p2) == GO.union(GO.GEOS(), p1, p2)
            @test GO.difference(p1, p2) == GO.difference(GO.GEOS(), p1, p2)
            =#
            # Instead, we just test that they work using `@test_nowarn`.
            @test_nowarn GO.intersection(GO.GEOS(), p1, p2)
            @test_nowarn GO.union(GO.GEOS(), p1, p2)
            @test_nowarn GO.difference(GO.GEOS(), p1, p2)
        end
        @testset "Segmentize" begin
            @test GI.npoint(GO.segmentize(GO.GEOS(; max_distance = 0.1), l2)) == GI.npoint(GO.segmentize(l2; max_distance = 0.1))
        end
        @testset "Simplify" begin
            @test GI.npoint(GO.simplify(GO.GEOS(; tol = 1.0, method = :TopologyPreserve), l2)) == GI.npoint(GO.simplify(l2; tol = 1.0))
            @test GI.npoint(GO.simplify(GO.GEOS(; tol = 1.0, method = :DouglasPeucker), l2)) == GI.npoint(GO.simplify(l2; tol = 1.0))
        end
    end
end
