using Test

import ArchGDAL as AG
import GeometryBasics as GB 
import GeoFormatTypes as GFT
import GeometryOps as GO 
import GeoInterface as GI
import LibGEOS as LG
import Proj
import Shapefile
import DataFrames, Tables, DataAPI
using Downloads: download
using ..TestHelpers

pv1 = [(1, 2), (3, 4), (5, 6), (1, 2)]
pv2 = [(3, 4), (5, 6), (6, 7), (3, 4)]
lr1 = GI.LinearRing(pv1)
lr2 =  GI.LinearRing(pv2)
poly = GI.Polygon([lr1, lr2])

@testset "apply" begin

    @testset_implementations "apply flip" begin
        flipped_poly = GO.apply(GI.PointTrait, $poly) do p
            (GI.y(p), GI.x(p))
        end
        @test flipped_poly == GI.Polygon([GI.LinearRing([(2, 1), (4, 3), (6, 5), (2, 1)]), 
                                          GI.LinearRing([(4, 3), (6, 5), (7, 6), (4, 3)])])
    end

    @testset "Tables.jl support" begin
        # check to account for missing data
        missing_or_equal(x, y) = (ismissing(x) && ismissing(y)) || (x == y)
        # file setup
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
                        @test all(missing_or_equal.(GO.Tables.getcolumn(centroid_table, column), GO.Tables.getcolumn(countries_table, column)))
                    end
                end
            end

            @testset "DataFrames" begin
                countries_df = DataFrames.DataFrame(countries_table)
                GO.DataAPI.metadata!(countries_df, "note metadata", "note metadata value"; style = :note)
                GO.DataAPI.metadata!(countries_df, "default metadata", "default metadata value"; style = :default)
                centroid_df = GO.apply(GO.centroid, GO.TraitTarget(GI.PolygonTrait(), GI.MultiPolygonTrait()), countries_df; crs = GFT.EPSG(3031));
                # Test that the Tables.jl materializer is used
                @test centroid_df isa DataFrames.DataFrame
                # Test that the centroids are correct
                @test all(centroid_df.geometry .== GO.centroid.(countries_df.geometry))
                @testset "Columns are preserved" begin  
                    for column in filter(!=(:geometry), GO.Tables.columnnames(countries_df))
                        @test all(missing_or_equal.(centroid_df[!, column], countries_df[!, column]))
                    end
                end
                @testset "Metadata preservation (or not)" begin
                    @test DataAPI.metadata(centroid_df, "note metadata") == "note metadata value"
                    @test !("default metadata" in DataAPI.metadatakeys(centroid_df))
                    @test DataAPI.metadata(centroid_df, "GEOINTERFACE:geometrycolumns") == (:geometry,)
                    @test DataAPI.metadata(centroid_df, "GEOINTERFACE:crs") == GFT.EPSG(3031)
                end
                @testset "Multiple geometry columns in metadata" begin
                    # set up a dataframe with multiple geometry columns
                    countries_df2 = deepcopy(countries_df)
                    countries_df2.centroid = GO.centroid.(countries_df2.geometry)
                    GI.DataAPI.metadata!(countries_df2, "GEOINTERFACE:geometrycolumns", (:geometry, :centroid); style = :note)
                    transformed = GO.transform(p -> p .+ 3, countries_df2)
                    @test GI.DataAPI.metadata(transformed, "GEOINTERFACE:geometrycolumns") == (:geometry, :centroid)
                    @test GI.DataAPI.metadata(transformed, "GEOINTERFACE:crs") == GFT.EPSG(4326)
                    # Test that the transformation was actually applied to both geometry columns.
                    @test all(map(zip(countries_df2.geometry, transformed.geometry)) do (o, n)
                        GO.equals(GO.transform(p -> p .+ 3, o), n)
                    end)
                    @test all(map(zip(countries_df2.centroid, transformed.centroid)) do (o, n)
                        any(isnan, o) || GO.equals(GO.transform(p -> p .+ 3, o), n)
                    end)
                end
            end
        end
        end
        @testset "Wrong geometry column kwarg" begin
            tab = Tables.dictcolumntable((; geometry = [(1, 2), (3, 4), (5, 6)], other = [1, 2, 3]))
            @test_throws "got a Float64" GO.transform(identity, tab; geometrycolumn = 1000.0)
            @test_throws "but the table has columns" GO.transform(identity, tab; geometrycolumn = :somethingelse)
        end
    end
end

@testset "unwrap" begin
    flipped_vectors = GO.unwrap(GI.PointTrait, poly) do p
        (GI.y(p), GI.x(p))
    end

    @test flipped_vectors == [[(2, 1), (4, 3), (6, 5), (2, 1)], [(4, 3), (6, 5), (7, 6), (4, 3)]]
end

@testset "flatten" begin
    very_wrapped = [[GI.FeatureCollection([GI.Feature(poly; properties=(;))])]]
    @test GO._tuple_point.(GO.flatten(GI.PointTrait, very_wrapped)) == vcat(pv1, pv2)
    @test collect(GO.flatten(GI.AbstractCurveTrait, [poly])) == [lr1, lr2]
    @test collect(GO.flatten(GI.x, GI.PointTrait, very_wrapped)) == first.(vcat(pv1, pv2))
    @testset "flatten with tables" begin
        # Construct a simple table with a geometry column
        geom_column = [GI.Point(1.0,1.0), GI.Point(2.0,2.0), GI.Point(3.0,3.0)]
        table = (geometry = geom_column, id = [1, 2, 3])
        
        # Test flatten on the table
        flattened = collect(GO.flatten(GI.PointTrait, table))
        
        @test length(flattened) == 3
        @test all(p isa GI.Point for p in flattened)
        @test flattened == geom_column
        
        # Test flatten with a function
        flattened_coords = collect(GO.flatten(p -> (GI.x(p), GI.y(p)), GI.PointTrait, table))
        
        @test length(flattened_coords) == 3
        @test all(c isa Tuple{Float64,Float64} for c in flattened_coords)
        @test flattened_coords == [(1.0,1.0), (2.0,2.0), (3.0,3.0)]
    end
end

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
