using Test
 
import GeoInterface as GI, GeometryOps as GO
using CoordinateTransformations

@testset "transform" begin
    geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                       GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])
    translated = GI.Polygon([GI.LinearRing([[4.5, 3.5], [6.5, 5.5], [8.5, 7.5], [4.5, 3.5]]), 
                             GI.LinearRing([[6.5, 5.5], [8.5, 7.5], [9.5, 8.5], [6.5, 5.5]])])
    f = CoordinateTransformations.Translation(3.5, 1.5)
    @test GO.equals(GO.transform(f, geom), translated)
end
