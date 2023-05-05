using Test

import GeoInterface as GI
import GeometryOps as GO

pv1 = [(1, 2), (3, 4), (5, 6), (1, 2)]
pv2 = [(3, 4), (5, 6), (6, 7), (3, 4)]
lr1 = GI.LinearRing(pv1)
lr2 =  GI.LinearRing(pv2)
poly = GI.Polygon([lr1, lr2])

@testset "apply" begin

    flipped_poly = GO.apply(GI.PointTrait, poly) do p
        (GI.y(p), GI.x(p))
    end

    @test flipped_poly == GI.Polygon([GI.LinearRing([(2, 1), (4, 3), (6, 5), (2, 1)]), 
                                      GI.LinearRing([(4, 3), (6, 5), (7, 6), (4, 3)])])
end

@testset "unwrap" begin
    flipped_vectors = GO.unwrap(GI.PointTrait, poly) do p
        (GI.y(p), GI.x(p))
    end

    @test flipped_vectors == [[(2, 1), (4, 3), (6, 5), (2, 1)], [(4, 3), (6, 5), (7, 6), (4, 3)]]
end

@testset "flatten" begin
    very_wrapped = [[GI.FeatureCollection([GI.Feature(poly; properties=(;))])]]
    @test collect(GO.flatten(GI.PointTrait, very_wrapped)) == vcat(pv1, pv2)
    @test collect(GO.flatten(GI.LinearRingTrait, [poly])) == [lr1, lr2]
end

@testset "reconstruct" begin
    revlr1 =  GI.LinearRing(reverse(pv2))
    revlr2 = GI.LinearRing(reverse(pv1))
    revpoly = GI.Polygon([revlr1, revlr2])
    points = collect(GO.flatten(GI.PointTrait, poly))
    @test GO.reconstruct(poly, reverse(points)) == revpoly
end
