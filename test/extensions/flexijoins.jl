using Test
using FlexiJoins
using DataFrames
import GeometryOps as GO
import GeoInterface as GI
using .TestHelpers

points = GI.MultiPoint(tuple.(rand(100), rand(100)))

pl = GI.Polygon([GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 0)])])
pu = GI.Polygon([GI.LinearRing([(0, 0), (0, 1), (1, 1), (0, 0)])])

@testset_implementations "Polygon DataDrame" begin
    points_df = DataFrame(geometry=collect(GI.getpoint($points)))
    poly_df = DataFrame(geometry=[$pl, $pu], color=[:red, :blue])
    # Test that the join happened correctly
    joined_df = FlexiJoins.innerjoin((poly_df, points_df), by_pred(:geometry, GO.contains, :geometry))
    @test all(GO.contains.(($pl,), joined_df.geometry_1[joined_df.color .== :red]))
    @test all(GO.contains.(($pu,), joined_df.geometry_1[joined_df.color .== :blue]))
    # Test that within also works
    @test_nowarn joined_df = FlexiJoins.innerjoin((points_df, poly_df), by_pred(:geometry, GO.within, :geometry))
end

