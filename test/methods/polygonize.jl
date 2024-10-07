using Test
import GeometryOps: polygonize
import DimensionalData as DD
import Rasters, DimensionalData
using GeometryOps, GeoInterface, Test
using ..TestHelpers

import GeometryOps as GO 
import GeoInterface as GI

@testset "Polygonize with xs and ys, without offsetarrays" begin
    @test !(@isdefined OffsetArrays) # to make sure this isn't loaded somewhere else
    data  = rand(1:4, 100, 100) .== 1
    unitrange = 1:100
    steprange = 1:1:100
    steprangelen = range(1, 100; length = 100)
    data_mp = polygonize(data)
    for range in (unitrange, steprange, steprangelen)
        data_mp_range = polygonize(range, range, data)
        @test GO.equals(data_mp, data_mp_range)
    end
    # ideally we'd have a better test to make sure this returns what we think it does
    data_mp_range200 = polygonize(2:2:200, 2:2:200, data)
    @test length(GI.coordinates(data_mp_range200)) == length(GI.coordinates(data_mp))

    # this is an example that could throw floating point error
    range_floats = -1.333333333333343:0.041666666666666664:0.374999999999986
    data2 = rand(1:4, length(range_floats), length(range_floats)) .== 1
    data_mp_range_floats = polygonize(range_floats, range_floats, data2)
end

import OffsetArrays

@testset "Polygonize with xs and ys, with offsetarrays" begin
    data  = rand(1:4, 100, 100) .== 1
    unitrange = 1:100
    steprange = 1:1:100
    steprangelen = range(1, 100; length = 100)
    data_mp = polygonize(data)
    for range in (unitrange, steprange, steprangelen)
        data_mp_range = polygonize(range, range, data)
        @test GO.equals(data_mp, data_mp_range)
    end
end

# Missing holes throw a warning, so testing there are
# no warnings in a range of randomisation is one way to test 
# that things are working, without testing explicit return values
for i in (100, 300), j in (100, 300)
    @testset "bool arrays without a function return MultiPolygon" begin
        A = rand(Bool, i, j)
        @test_nowarn multipolygon = polygonize(A);
        @test multipolygon isa GI.MultiPolygon
        @test GI.ngeom(multipolygon) > 0
    end

    A = rand(i, j)
    @testset "bool functions return MultiPolygon" begin
        multipolygon = @test_nowarn polygonize(>(0.5), A);
        @test multipolygon isa GI.MultiPolygon
        @test GI.ngeom(multipolygon) > 0
    end

    @testset "other functions return FeatureCollection" begin
        fc = @test_nowarn polygonize(x -> trunc(3x), A);
        @test fc isa GI.FeatureCollection
        @test GI.nfeature(fc) == 3
        @test map(GI.getfeature(fc)) do f
            GI.properties(f).value
        end == [0.0, 1.0, 2.0]
    end

    @testset "values are polygonized without a function" begin
        A = rand(1:3, i, j)
        fc = @test_nowarn polygonize(A)
        fc isa GI.FeatureCollection
        @test GI.nfeature(fc) == 3
        @test map(GI.getfeature(fc)) do f
            GI.properties(f).value
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
end

@testset "Polygonize with DimensionalData compatible arrays" begin
    data = rand(1:4, 100, 50) .== 1
    dd = DD.DimArray(data, (DD.X(51:150), DD.Y(151:200)))
    @testset "DimensionalData" begin
        data_mp = polygonize(51:150, 151:200, data);
        dd_mp = polygonize(dd);
        @test GO.equals(data_mp, dd_mp)
    end

    @testset "Rasters" begin
        data = rand(1:4, 100, 50) .== 1
        rast = Rasters.Raster(data; dims=(DD.X(51:150), DD.Y(151:200)), crs=Rasters.GeoFormatTypes.EPSG(4326))
        data_mp = polygonize(51:150, 151:200, data)
        rast_mp = @test_nowarn polygonize(rast)
        @test GO.equals(data_mp, rast_mp)
        @test GI.crs(rast_mp) == GI.crs(evil)
    end
end

