using Test
 
import GeoInterface as GI, GeometryOps as GO

@testset "flip" begin
    geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                       GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])


    @test GO.flip(geom) == GI.Polygon([GI.LinearRing([(2, 1), (4, 3), (6, 5), (2, 1)]), 
                                       GI.LinearRing([(4, 3), (6, 5), (7, 6), (4, 3)])])
end
