using Test

@testset "Area without Proj" begin
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

        planar_poly = GI.Polygon([[
            (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0),
        ]])
        @test GO.area(planar_poly) == GO.area(GO.Planar(), planar_poly) == 1.0
    """
    test_project = dirname(Base.active_project())
    cmd = `$(Base.julia_cmd()) --project=$test_project -e $child_test`
    @test success(cmd)
end
