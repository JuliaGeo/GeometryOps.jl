using Test
using FlexiJoins
using DataFrames
import GeometryOps as GO
import GeoInterface as GI
import GeoFormatTypes as GFT
import Proj
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
    @test_throws "type NamedTuple has no field" FlexiJoins.innerjoin((no_metadata_location, zones_by_shape), GO.within)
end

@testset "Join manifold from the CRS" begin
    # A box across the antimeridian: on the sphere it holds the points at ±179.5,
    # while its planar lon/lat polygon spans the other way round and holds (0, 0).
    box = GI.Polygon([[(170.0, -5.0), (-170.0, -5.0), (-170.0, 5.0), (170.0, 5.0), (170.0, -5.0)]])
    points = [GI.Point(179.5, 0.0), GI.Point(0.0, 0.0), GI.Point(-179.5, 0.0)]
    with_crs(df, crs) = (GI.DataAPI.metadata!(df, "GEOINTERFACE:crs", crs; style=:note); df)

    box_df = with_crs(DataFrame(geometry=[box]), GFT.EPSG(4326))
    points_df = with_crs(DataFrame(geometry=points, id=1:3), GFT.EPSG(4326))
    for mode in (FlexiJoins.Mode.Tree(), FlexiJoins.Mode.NestedLoopFast(), FlexiJoins.Mode.NestedLoop())
        joined = FlexiJoins.innerjoin((points_df, box_df), by_pred(:geometry, GO.within, :geometry); mode)
        @test sort(joined.id) == [1, 3]
        joined = FlexiJoins.innerjoin((box_df, points_df), by_pred(:geometry, GO.contains, :geometry); mode)
        @test sort(joined.id) == [1, 3]
    end
    @test sort(FlexiJoins.innerjoin((points_df, box_df), GO.within).id) == [1, 3]

    planar_box_df = DataFrame(geometry=[box])
    planar_points_df = DataFrame(geometry=points, id=1:3)
    @test FlexiJoins.innerjoin((planar_points_df, planar_box_df), GO.within).id == [2]

    @test_throws ArgumentError FlexiJoins.innerjoin((planar_points_df, box_df), GO.within)
    projected_points_df = with_crs(DataFrame(geometry=points, id=1:3), GFT.EPSG(3857))
    @test_throws ArgumentError FlexiJoins.innerjoin((projected_points_df, box_df), GO.within)
    @test FlexiJoins.innerjoin((projected_points_df, with_crs(DataFrame(geometry=[box]), GFT.EPSG(3857))), GO.within).id == [2]
end

