using Test

import GeoInterface as GI
import GeometryOps as GO
import GeometryBasics as GB

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
    reconstructed = GO.reconstruct(poly, reverse(points))
    @test reconstructed == revpoly
    @test reconstructed isa GI.Polygon


    gb_revlr1 = GB.LineString(GB.Point.(reverse(pv2)))
    gb_revlr2 = GB.LineString(GB.Point.(reverse(pv1)))
    gb_revpoly = GB.Polygon(gb_revlr1, [gb_revlr2])
    gb_lr1 = GB.LineString(GB.Point.(pv1))
    gb_lr2 = GB.LineString(GB.Point.(pv2))
    gb_poly = GB.Polygon(gb_lr1, [gb_lr2])
    gb_points = collect(GO.flatten(GI.PointTrait, gb_poly))
    gb_reconstructed = GO.reconstruct(gb_poly, reverse(gb_points))
    @test gb_reconstructed == gb_revpoly
    @test gb_reconstructed isa GB.Polygon
end
