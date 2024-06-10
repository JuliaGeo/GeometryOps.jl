import GeometryOps as GO, GeoInterface as GI
using FlexiJoins, DataFrames

points = tuple.(rand(100), rand(100))
points_df = DataFrame(geometry = points)

pl = GI.Polygon([GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 0)])])
pu = GI.Polygon([GI.LinearRing([(0, 0), (0, 1), (1, 1), (0, 0)])])

@test_all_implementations "Polygon DataDrame" (pl, pu) begin
    poly_df = DataFrame(geometry = [pl, pu], color = [:red, :blue])
    # Test that the join happened correctly
    joined_df = FlexiJoins.innerjoin((poly_df, points_df), by_pred(:geometry, GO.contains, :geometry))
    @test all(GO.contains.((pl,), joined_df.geometry_1[joined_df.color .== :red]))
    @test all(GO.contains.((pu,), joined_df.geometry_1[joined_df.color .== :blue]))
    # Test that within also works
    @test_nowarn joined_df = FlexiJoins.innerjoin((points_df, poly_df), by_pred(:geometry, GO.within, :geometry))
end

