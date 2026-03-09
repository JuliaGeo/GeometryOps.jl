using Test
import GeoJSON, JLD2
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
using ..TestHelpers

datadir = realpath(joinpath(dirname(pathof(GO)), "../test/data"))
@testset "RadialDistance and VisvalingamWhyatt" begin
    fc = GeoJSON.read(joinpath(datadir, "simplify.geojson"))
    fc2 = GeoJSON.read(joinpath(datadir, "simplify2.geojson"))
    fcs = [fc for i in 1:100]

    # TODO: @testset_implementations doesn't handle feature collections yet
    for T in (GO.RadialDistance, GO.VisvalingamWhyatt)
        @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(T(number=10), fc)))) == 10
        @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(T(ratio=0.5), fc)))) == 39 # Half of 78
        GO.simplify(T(tol=0.001), fc; threaded=true, calc_extent=true)
        GO.simplify(T(tol=0.001), fcs; threaded=true, calc_extent=true)
    end
end

@testset "DouglasPeucker" begin
    # Test for issue #386: BoundsError when simplifying small geometries with low number/ratio
    @testset "small geometry simplification (issue #386)" begin
        # This would cause a BoundsError before the fix due to indexing bug in the while loop
        line = GI.LineString([(rand(), rand()) for _ in 1:4])
        @test_nowarn GO.simplify(line; ratio=0.1)
        @test_nowarn GO.simplify(line; tol=0.1)
        @test_nowarn GO.simplify(line; number=3)
        # Verify the output is valid
        result = GO.simplify(line; number=3)
        @test GI.npoint(result) == 3
    end

    poly_coords = JLD2.jldopen(joinpath(datadir, "complex_polygons.jld2"))["verts"][1:4]
    for c in poly_coords
        npoints = length(c[1])
        poly = LG.Polygon(c)
        @testset_implementations "Polygon coords match LibGEOS simplify" begin
            lg_vals = GI.coordinates(LG.simplify($poly, 100.0))[1]
            reduced_npoints = length(lg_vals)
            @test all(GI.coordinates(GO.simplify($poly; tol = 100.0))[1] .== lg_vals)
            @test all(GI.coordinates(GO.simplify($poly; number = reduced_npoints))[1] .== lg_vals)
            @test all(GI.coordinates(GO.simplify($poly; ratio = (reduced_npoints/npoints)))[1] .== lg_vals)
        end
    end
    # Ensure last point isn't removed with curve
    c = poly_coords[1]
    linestring = LG.LineString(c[1])
    @testset_implementations "LineString coords match LibGEOS simplify" begin
        lg_vals = GI.coordinates(LG.simplify($linestring, 100.0))
        @test all(GI.coordinates(GO.simplify($linestring; tol = 100.0)) .== lg_vals)
    end
end
