using Test
using FlexiJoins
using DataFrames
import GeometryOps as GO
import GeoInterface as GI
using GeometryOpsTestHelpers

points = GI.MultiPoint(tuple.(rand(100), rand(100)))

pl = GI.Polygon([GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 0)])])
pu = GI.Polygon([GI.LinearRing([(0, 0), (0, 1), (1, 1), (0, 0)])])

@testset_implementations "Polygon DataFrame" begin
    points_df = DataFrame(geometry=collect(GI.getpoint($points)))
    poly_df = DataFrame(geometry=[$pl, $pu], color=[:red, :blue])
    # Test that the join happened correctly
    joined_df = FlexiJoins.innerjoin((poly_df, points_df), by_pred(:geometry, GO.contains, :geometry))
    @test all(GO.contains.(($pl,), joined_df.geometry_1[joined_df.color .== :red]))
    @test all(GO.contains.(($pu,), joined_df.geometry_1[joined_df.color .== :blue]))
    # Test that within also works
    @test_nowarn joined_df = FlexiJoins.innerjoin((points_df, poly_df), by_pred(:geometry, GO.within, :geometry))

    points_by_location = DataFrame(location=points_df.geometry)
    zones_by_shape = DataFrame(shape=poly_df.geometry, color=poly_df.color)
    GI.DataAPI.metadata!(points_by_location, "GEOINTERFACE:geometrycolumns", (:location,); style=:note)
    GI.DataAPI.metadata!(zones_by_shape, "GEOINTERFACE:geometrycolumns", (:shape,); style=:note)
    @test GI.geometrycolumns(points_by_location) == (:location,)
    @test GI.geometrycolumns(zones_by_shape) == (:shape,)
    @test nrow(joined_df) == nrow(FlexiJoins.innerjoin((points_by_location, zones_by_shape), GO.within))

    GI.DataAPI.metadata!(points_by_location, "GEOINTERFACE:geometrycolumns", (:location, :alternate_location); style=:note)
    @test_logs (:warn, r"First input declares multiple geometry columns \(:location, :alternate_location\); using the first.") FlexiJoins.innerjoin((points_by_location, zones_by_shape), GO.within)
    @test nrow(joined_df) == nrow(FlexiJoins.innerjoin((points_by_location, zones_by_shape), GO.within))

    # GeoInterface defaults a DataFrame with no geometry metadata to `:geometry`.
    # This table only has `:location`, so the shorthand must reject it before delegation.
    no_metadata_location = DataFrame(location=points_df.geometry)
    @test GI.geometrycolumns(no_metadata_location) == (:geometry,)
    @test_throws "type NamedTuple has no field `geometry`" FlexiJoins.innerjoin((no_metadata_location, zones_by_shape), GO.within)
end
