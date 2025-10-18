import GeometryOps as GO
import GeoInterface as GI
using .TestHelpers

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
