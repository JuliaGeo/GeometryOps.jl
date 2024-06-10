using Test

import GeoInterface as GI
import GeometryOps as GO
import GeometryBasics as GB
import Proj
import Shapefile, DataFrames

pv1 = [(1, 2), (3, 4), (5, 6), (1, 2)]
pv2 = [(3, 4), (5, 6), (6, 7), (3, 4)]
lr1 = GI.LinearRing(pv1)
lr2 =  GI.LinearRing(pv2)
poly = GI.Polygon([lr1, lr2])

@testset "apply" begin

    @test_all_implementations "simple flip to tuple" poly begin
        flipped_poly = GO.apply(GI.PointTrait, poly) do p
            (GI.y(p), GI.x(p))
        end

        @test flipped_poly == GI.Polygon([GI.LinearRing([(2, 1), (4, 3), (6, 5), (2, 1)]), 
                                          GI.LinearRing([(4, 3), (6, 5), (7, 6), (4, 3)])])
    end

    @testset "Tables.jl support" begin
        mktempdir() do dir
        cd(dir) do

            download("https://rawcdn.githack.com/nvkelso/natural-earth-vector/v5.1.2/110m_cultural/ne_110m_admin_0_countries.shp", "countries.shp")
            download("https://rawcdn.githack.com/nvkelso/natural-earth-vector/v5.1.2/110m_cultural/ne_110m_admin_0_countries.shx", "countries.shx")
            download("https://rawcdn.githack.com/nvkelso/natural-earth-vector/v5.1.2/110m_cultural/ne_110m_admin_0_countries.dbf", "countries.dbf")
            download("https://rawcdn.githack.com/nvkelso/natural-earth-vector/v5.1.2/110m_cultural/ne_110m_admin_0_countries.prj", "countries.prj")
            countries_table = Shapefile.Table("countries.shp")

            @testset "Shapefile" begin
                centroid_table = GO.apply(GO.centroid, GO.TraitTarget(GI.PolygonTrait(), GI.MultiPolygonTrait()), countries_table);
                centroid_geometry = centroid_table.geometry
                # Test that the centroids are correct
                @test all(centroid_geometry .== GO.centroid.(countries_table.geometry))
                @testset "Columns are preserved" begin  
                    for column in Iterators.filter(!=(:geometry), GO.Tables.columnnames(countries_table))
                        @test all(GO.Tables.getcolumn(centroid_table, column) .== GO.Tables.getcolumn(countries_table, column))
                    end
                end
            end

            @testset "DataFrames" begin
                countries_df = DataFrames.DataFrame(countries_table)
                centroid_df = GO.apply(GO.centroid, GO.TraitTarget(GI.PolygonTrait(), GI.MultiPolygonTrait()), countries_df);
                @test centroid_df isa DataFrames.DataFrame
                centroid_geometry = centroid_df.geometry
                # Test that the centroids are correct
                @test all(centroid_geometry .== GO.centroid.(countries_df.geometry))
                @testset "Columns are preserved" begin  
                    for column in Iterators.filter(!=(:geometry), GO.Tables.columnnames(countries_df))
                        @test all(centroid_df[!, column] .== countries_df[!, column])
                    end
                end
            end
        end
        end
    end
end



@test_all_implementations "unwrap" poly begin
    flipped_vectors = GO.unwrap(GI.PointTrait, poly) do p
        (GI.y(p), GI.x(p))
    end

    @test flipped_vectors == [[(2, 1), (4, 3), (6, 5), (2, 1)], [(4, 3), (6, 5), (7, 6), (4, 3)]]
end

@test_all_implementations "flatten" (poly, lr1, lr2) begin
    very_wrapped = [[GI.FeatureCollection([GI.Feature(poly; properties=(;))])]]
    @test GO._tuple_point.(GO.flatten(GI.PointTrait, very_wrapped)) == vcat(pv1, pv2)
    @test collect(GO.flatten(GI.AbstractCurveTrait, [poly])) == [lr1, lr2]
    @test collect(GO.flatten(GI.x, GI.PointTrait, very_wrapped)) == first.(vcat(pv1, pv2))
end

# TODO test_all_implementations
@testset "reconstruct" begin
    revlr1 =  GI.LinearRing(reverse(pv2))
    revlr2 = GI.LinearRing(reverse(pv1))
    revpoly = GI.Polygon([revlr1, revlr2])
    points = collect(GO.flatten(GI.PointTrait, poly))
    reconstructed = GO.reconstruct(poly, reverse(points))
    @test reconstructed == revpoly
    @test reconstructed isa GI.Polygon


    revlr1 = GI.LineString(GB.Point.(reverse(pv2)))
    revlr2 = GI.LineString(GB.Point.(reverse(pv1)))
    revpoly = GI.Polygon([revlr1, revlr2])
    gb_lr1 = GB.LineString(GB.Point.(pv1))
    gb_lr2 = GB.LineString(GB.Point.(pv2))
    gb_poly = GB.Polygon(gb_lr1, [gb_lr2])
    gb_points = collect(GO.flatten(GI.PointTrait, gb_poly))
    gb_reconstructed = GO.reconstruct(gb_poly, reverse(gb_points))
    @test gb_reconstructed == revpoly
    @test gb_reconstructed isa GI.Polygon
end
