using Test
 
import GeoInterface as GI
import GeometryOps as GO
import GeoJSON

@testset "simplify" begin
    datadir = realpath(joinpath(dirname(pathof(GO)), "../test/data"))
    fc = GeoJSON.read(joinpath(datadir, "simplify.geojson"))
    fc2 = GeoJSON.read(joinpath(datadir, "simplify2.geojson"))

    GO.simplify(GO.RadialDistance(tol=0.0001), fc) 
    GO.simplify(GO.DouglasPeucker(tol=0.0001), fc) 

    @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(GO.VisvalingamWhyatt(number=10), fc)))) == 10
    @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(GO.VisvalingamWhyatt(ratio=0.5), fc)))) == 39
    length(collect(GO.flatten(GI.PointTrait, GO.simplify(GO.VisvalingamWhyatt(tol=0.1), fc))))

    # DouglasPeucker is the default
	@test_broken collect(GO.flatten(GI.PointTrait, GO.simplify(fc; tol=0.0001))) [
        (26.14843, -28.297552),
        (26.150354, -28.302606),
        (26.135463, -28.304283),
        (26.14843, -28.297552),
    ]

	@test_throws ErrorException GO.simplify(fc2, tol=-1)
    @test_broken collect(GO.flatten(GI.PointTrait, GO.simplify(fc2; tol=0.01, prefilter=true))) == [
        (179.975281, -16.51477),
        (179.980431, -16.539127),
        (180.0103, -16.523328),
        (180.007553, -16.534848),
        (180.018196, -16.539127),
        (180.061455, -16.525632),
        (180.066605, -16.513124),
        (180.046349, -16.479547),
        (180.086861, -16.44761),
        (180.084114, -16.441354),
        (180.055618, -16.439707),
        (180.026093, -16.464732),
        (180.01442, -16.464073),
        (179.975281, -16.51477),
    ]
end
