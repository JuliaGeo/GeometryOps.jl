using GeometryOps, GeoInterface, Test

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
