using Test

@testset "AutoManifold without Proj" begin
    child_test = """
        using Test
        import GeoFormatTypes
        import GeoInterface as GI
        import GeometryOps as GO

        @test Base.get_extension(GO, :GeometryOpsProjExt) === nothing

        unknown_crs_poly = GI.Polygon(
            [[
                (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0),
            ]];
            crs=GeoFormatTypes.EPSG(4326),
        )
        @test GI.crstrait(unknown_crs_poly) isa GI.UnknownTrait
        @test GO.area(unknown_crs_poly) == GO.area(GO.Planar(), unknown_crs_poly)
        @test GO.area(unknown_crs_poly, Float32) isa Float32
        @test GO.perimeter(unknown_crs_poly) == 4.0
        @test GO.distance((0.5, 2.0), unknown_crs_poly) == 1.0
        @test GI.npoint(GO.segmentize(unknown_crs_poly; max_distance = 0.5)) == 9

        planar_poly = GI.Polygon([[
            (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0),
        ]])
        @test GO.area(planar_poly) == GO.area(GO.Planar(), planar_poly) == 1.0

        # A geometry whose CRS trait is geographic is interpreted as lon/lat on a sphere.
        struct GeographicLineString
            points::Vector{Tuple{Float64, Float64}}
        end
        GI.isgeometry(::Type{GeographicLineString}) = true
        GI.geomtrait(::GeographicLineString) = GI.LineStringTrait()
        GI.ngeom(::GI.LineStringTrait, line::GeographicLineString) = length(line.points)
        GI.getgeom(::GI.LineStringTrait, line::GeographicLineString, i) = line.points[i]
        GI.crs(::GeographicLineString) = GeoFormatTypes.EPSG(4326)
        GI.crstrait(::GeographicLineString) = GI.GeographicTrait()

        points = [(0.0, 0.0), (90.0, 0.0)]
        geographic_line = GeographicLineString(points)
        @test GO.perimeter(geographic_line) == GO.perimeter(GO.Spherical(), GI.LineString(points))
        @test GO.distance((45.0, 30.0), geographic_line) ≈ GO.Spherical().radius * deg2rad(30.0)
        @test GI.npoint(GO.segmentize(geographic_line; max_distance = 5_100_000)) == 3
    """
    test_project = dirname(Base.active_project())
    cmd = `$(Base.julia_cmd()) --project=$test_project -e $child_test`
    @test success(cmd)
end
