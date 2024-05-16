import GeoJSON
import JLD2
import LibGEOS as LG
import GeometryOps as GO
import GeoInterface as GI

datadir = realpath(joinpath(dirname(pathof(GO)), "../test/data"))
@testset "RadialDistance and VisvalingamWhyatt" begin
    fc = GeoJSON.read(joinpath(datadir, "simplify.geojson"))
    fc2 = GeoJSON.read(joinpath(datadir, "simplify2.geojson"))
    fcs = [fc for i in 1:100]

    # TODO: @test_all_implementations doesn't handle feature collections yet
    for T in (GO.RadialDistance, GO.VisvalingamWhyatt)
        @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(T(number=10), fc)))) == 10
        @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(T(ratio=0.5), fc)))) == 39 # Half of 78
        GO.simplify(T(tol=0.001), fc; threaded=true, calc_extent=true)
        GO.simplify(T(tol=0.001), fcs; threaded=true, calc_extent=true)
    end
end
@testset "DouglasPeucker" begin
    poly_coords = JLD2.jldopen(joinpath(datadir, "complex_polygons.jld2"))["verts"][1:4]
    for c in poly_coords
        npoints = length(c[1])
        poly = LG.Polygon(c)
        lg_vals = GI.coordinates(LG.simplify(poly, 100.0))[1]
        reduced_npoints = length(lg_vals)
        @test_all_implementations "Polygon coords match LibGEOS simplify" poly begin
            @test all(GI.coordinates(GO.simplify(poly; tol = 100.0))[1] .== lg_vals)
            @test all(GI.coordinates(GO.simplify(poly; number = reduced_npoints))[1] .== lg_vals)
            @test all(GI.coordinates(GO.simplify(poly; ratio = (reduced_npoints/npoints)))[1] .== lg_vals)
        end
    end
    # Ensure last point isn't removed with curve
    c = poly_coords[1]
    linestring = LG.LineString(c[1])
    lg_vals = GI.coordinates(LG.simplify(linestring, 100.0))
    @test_all_implementations "LineString coords match LibGEOS simplify" linestring begin
        @test all(GI.coordinates(GO.simplify(linestring; tol = 100.0)) .== lg_vals)
    end
end
