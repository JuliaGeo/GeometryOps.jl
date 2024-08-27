using GeometryOps, GeoInterface, Test
import OffsetArrays, DimensionalData, Rasters
using ..TestHelpers

# Missing holes throw a warning, so testing there are
# no warnings in a range of randomisation is one way to test 
# that things are working, without testing explicit return values
for i in (100, 300), j in (100, 300)
    @testset "bool arrays without a function return MultiPolygon" begin
        A = rand(Bool, i, j)
        @test_nowarn multipolygon = polygonize(A);
        @test multipolygon isa GeoInterface.MultiPolygon
        @test GeoInterface.ngeom(multipolygon) > 0
    end

    A = rand(i, j)
    @testset "bool functions return MultiPolygon" begin
        multipolygon = @test_nowarn polygonize(>(0.5), A);
        @test multipolygon isa GeoInterface.MultiPolygon
        @test GeoInterface.ngeom(multipolygon) > 0
    end

    @testset "other functions return FeatureCollection" begin
        fc = @test_nowarn polygonize(x -> trunc(3x), A);
        @test fc isa GeoInterface.FeatureCollection
        @test GeoInterface.nfeature(fc) == 3
        @test map(GeoInterface.getfeature(fc)) do f
            GeoInterface.properties(f).value
        end == [0.0, 1.0, 2.0]
    end

    @testset "values are polygonized without a function" begin
        A = rand(1:3, i, j)
        fc = @test_nowarn polygonize(A)
        fc isa GeoInterface.FeatureCollection
        @test GeoInterface.nfeature(fc) == 3
        @test map(GeoInterface.getfeature(fc)) do f
            GeoInterface.properties(f).value
        end == [1, 2, 3]
    end
end


@testset "Polygonize with exotic arrays" begin
    @testset "OffsetArrays" begin
        data = rand(1:4, 100, 100) .== 1
        evil = OffsetArrays.Origin(-100, -100)(data)
        data_mp = polygonize(data)
        evil_mp = @test_nowarn polygonize(evil)
        evil_in_data_space_mp = GO.transform(evil_mp) do point
            point .- evil.offsets # undo the offset from the OffsetArray
        end
        @test GO.equals(data_mp, evil_in_data_space_mp)
    end
    @testset "DimensionalData" begin
        data = rand(1:4, 100, 100) .== 1
        evil = DimensionalData.DimArray(data, (DimensionalData.X(-100:-1), DimensionalData.Y(-100:-1)))
        data_mp = polygonize(data)
        evil_mp = @test_nowarn polygonize(evil)
        @test GO.equals(data_mp, evil_mp)
    end
    @testset "Rasters" begin
        data = rand(1:4, 100, 100) .== 1
        evil = Rasters.Raster(data; dims = (DimensionalData.X(-100:-1), DimensionalData.Y(-90:9)), crs = Rasters.GeoFormatTypes.EPSG(4326))
        data_mp = polygonize(data)
        evil_mp = @test_nowarn polygonize(evil)
        @test GO.equals(data_mp, evil_mp)
        @test GI.crs(evil_mp) == GI.crs(evil)
    end
end
