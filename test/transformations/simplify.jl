using Test
 
import GeoInterface as GI
import GeometryOps as GO
import GeoJSON

# Unncomment when JSON3 bumps a patch version
# @testset "simplify" begin
#     datadir = realpath(joinpath(dirname(pathof(GO)), "../test/data"))
#     fc = GeoJSON.read(joinpath(datadir, "simplify.geojson"))
#     fc2 = GeoJSON.read(joinpath(datadir, "simplify2.geojson"))

#     for T in (GO.RadialDistance, GO.VisvalingamWhyatt, GO.DouglasPeucker)
#         @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(T(number=10), fc)))) == 10
#         @test length(collect(GO.flatten(GI.PointTrait, GO.simplify(T(ratio=0.5), fc)))) == 39 # Half of 78
#         GO.simplify(T(tol=0.001), fc)
#     end
# end
