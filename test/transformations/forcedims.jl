import GeometryOps as GO
import GeoInterface as GI
using GeometryOpsTestHelpers

using Test

geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])

@testset_implementations "force dimensions" begin

    # Test forcexy on 3D geometry
    geom3d = GO.transform($geom) do p
        (GI.x(p), GI.y(p), 1.0)
    end
    @test GI.is3d(geom3d)
    geom2d = GO.forcexy(geom3d)
    @test !GI.is3d(geom2d)
    @test GO.equals(geom2d, geom)

    # Test forcexyz with default z
    geom3d_default = GO.forcexyz($geom)
    @test GI.is3d(geom3d_default)
    points3d = collect(GO.flatten(GI.PointTrait, geom3d_default))
    @test all(p -> GI.z(p) == 0, points3d)

    # Test forcexyz with custom z
    geom3d_custom = GO.forcexyz($geom, 5.0)
    @test GI.is3d(geom3d_custom)
    points3d_custom = collect(GO.flatten(GI.PointTrait, geom3d_custom))
    @test all(p -> GI.z(p) == 5.0, points3d_custom)

    # Test forcexyz preserves existing z values
    @test GO.equals(GO.forcexyz(geom3d), geom3d)
end

@testset "forcexyz on mixed dimensions" begin
    line2d = GI.LineString([(1.0, 2.0), (3.0, 4.0)])
    line3d = GI.LineString([(5.0, 6.0, 7.0), (8.0, 9.0, 10.0)])
    filled = [(1.0, 2.0, 42.0), (3.0, 4.0, 42.0)]
    preserved = [(5.0, 6.0, 7.0), (8.0, 9.0, 10.0)]
    # Exercise both orders: the fill value must not depend on the previous point.
    for (geoms, expected) in (([line2d, line3d], vcat(filled, preserved)),
                              ([line3d, line2d], vcat(preserved, filled))),
        input in (geoms, GI.GeometryCollection(geoms))
        result = GO.forcexyz(input, 42)
        @test collect(GO.flatten(GI.PointTrait, result)) == expected
        @test all(GI.is3d, GO.flatten(GI.PointTrait, result))
    end
    @test (@inferred GO.forcexyz((1.0, 2.0), 42)) == (1.0, 2.0, 42.0)
    @test (@inferred GO.forcexyz((1.0, 2.0, 3.0), 42)) == (1.0, 2.0, 3.0)
end

@testset "forcexyz on a vector of geometries" begin
    geoms3d = GO.forcexyz([geom, geom])
    @test all(GI.is3d, geoms3d)
    @test all(p -> GI.z(p) == 0, GO.flatten(GI.PointTrait, geoms3d))
end
