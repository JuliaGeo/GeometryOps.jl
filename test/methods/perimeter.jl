using Test

import GeometryOps as GO, GeoInterface as GI, GeoFormatTypes as GFT
import Proj # make sure the GO extension on Proj is active

# First, test against R

@testset "R/SF readme" begin

    outer = GI.LinearRing([(0,0),(10,0),(10,10),(0,10),(0,0)])
    hole1 = GI.LinearRing([(1,1),(1,2),(2,2),(2,1),(1,1)])
    hole2 = GI.LinearRing([(5,5),(5,6),(6,6),(6,5),(5,5)])

    p = GI.Polygon([outer, hole1, hole2])
    mp = GI.MultiPolygon([
        p, 
        GO.transform(x -> x .+ 12, GI.Polygon([outer, hole1]))
    ])

    @test GO.perimeter(p) == 48
    @test GO.perimeter(mp) == 92
    @test GO.perimeter(GI.GeometryCollection([p, mp])) == 48+92
end
