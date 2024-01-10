using Test
# import GeoInterface as GI
# import GeometryOps as GO
import GeoJSON
import JLD2

@testset "simplify" begin
    datadir = realpath(joinpath(dirname(pathof(GO)), "../test/data"))
    fc = GeoJSON.read(joinpath(datadir, "simplify.geojson"))
    fc2 = GeoJSON.read(joinpath(datadir, "simplify2.geojson"))
    fcs = [fc for i in 1:100]

    for T in (GO.RadialDistance, GO.VisvalingamWhyatt)
        @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(T(number=10), fc)))) == 10
        @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(T(ratio=0.5), fc)))) == 39 # Half of 78
        GO.simplify(T(tol=0.001), fc; threaded=true, calc_extent=true)
        GO.simplify(T(tol=0.001), fcs; threaded=true, calc_extent=true)
    end

    poly_coords = JLD2.jldopen(joinpath(datadir, "complex_polygons.jld2"))["verts"][1:3]
    for c in poly_coords
        npoints = length(c[1])
        third = (npoints == 50)  # only the third is failing
        third && @show c
        poly = LG.Polygon(c)
        lg_vals = GI.coordinates(LG.simplify(poly, 100.0))[1]
        reduced_npoints = length(lg_vals)
        third && @show reduced_npoints
        tol_vals = GI.coordinates(GO.simplify(poly; tol = 100.0))[1]
        @test all(tol_vals .== lg_vals)
        @test all(GI.coordinates(GO.simplify(poly; number = reduced_npoints))[1] .== lg_vals)
        @test all(GI.coordinates(GO.simplify(poly; ratio = (reduced_npoints/npoints)))[1] .== lg_vals)
    end

end
