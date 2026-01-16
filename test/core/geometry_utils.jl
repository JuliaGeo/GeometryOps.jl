using Test
using GeometryOps
import GeometryOps as GO
import GeoInterface as GI
using Tables
using OffsetArrays

@testset "get_geometries" begin
    # Set up some test geometries
    pt1 = GI.Point(0.0, 0.0)
    pt2 = GI.Point(1.0, 1.0)
    pt3 = GI.Point(2.0, 2.0)
    line = GI.LineString([(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])
    poly = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])

    @testset "AbstractArray input" begin
        @testset "standard Vector" begin
            geoms = [pt1, pt2, pt3]
            result = GO.get_geometries(geoms)
            @test result === geoms  # Should return same object
            @test length(result) == 3
        end

        @testset "offset axes array" begin
            geoms = OffsetArray([pt1, pt2, pt3], -1:1)  # indices -1, 0, 1
            @test Base.has_offset_axes(geoms)
            result = GO.get_geometries(geoms)
            @test !Base.has_offset_axes(result)  # Should be 1-indexed now
            @test result isa Vector
            @test length(result) == 3
            @test result[1] == pt1
            @test result[2] == pt2
            @test result[3] == pt3
        end
    end

    @testset "Single geometry" begin
        @test GO.get_geometries(pt1) === pt1
        @test GO.get_geometries(line) === line
        @test GO.get_geometries(poly) === poly
    end

    @testset "GeometryCollection" begin
        @testset "simple collection" begin
            gc = GI.GeometryCollection([pt1, line, poly])
            result = GO.get_geometries(gc)
            @test result isa AbstractArray
            @test length(result) == 3
        end

        @testset "recursive processing" begin
            # GeometryCollection that when collected becomes an AbstractArray
            # The recursive call should handle this properly
            gc = GI.GeometryCollection([pt1, pt2])
            result = GO.get_geometries(gc)
            @test result isa AbstractArray
            @test length(result) == 2
        end
    end

    @testset "Feature input" begin
        feature = GI.Feature(poly; properties=(name="test",))
        result = GO.get_geometries(feature)
        # Should return the geometry from the feature
        @test GI.trait(result) isa GI.PolygonTrait
    end

    @testset "FeatureCollection input" begin
        f1 = GI.Feature(pt1; properties=(id=1,))
        f2 = GI.Feature(pt2; properties=(id=2,))
        f3 = GI.Feature(poly; properties=(id=3,))
        fc = GI.FeatureCollection([f1, f2, f3])

        result = GO.get_geometries(fc)
        @test result isa AbstractArray
        @test length(result) == 3
        @test GI.trait(result[1]) isa GI.PointTrait
        @test GI.trait(result[2]) isa GI.PointTrait
        @test GI.trait(result[3]) isa GI.PolygonTrait
    end

    @testset "Table input" begin
        @testset "default geometry column" begin
            # Create a simple table with a geometry column
            # Using NamedTuple which satisfies Tables.istable
            tbl = (
                geometry = [pt1, pt2, pt3],
                name = ["a", "b", "c"]
            )
            @test Tables.istable(tbl)

            # Define geometrycolumns for this table type
            result = GO.get_geometries(tbl; geometrycolumn=:geometry)
            @test result isa AbstractArray
            @test length(result) == 3
        end

        @testset "explicit geometrycolumn kwarg" begin
            tbl = (
                geom = [pt1, pt2],
                other_geom = [pt3, line],
                id = [1, 2]
            )

            result1 = GO.get_geometries(tbl; geometrycolumn=:geom)
            @test length(result1) == 2
            @test result1[1] == pt1

            result2 = GO.get_geometries(tbl; geometrycolumn=:other_geom)
            @test length(result2) == 2
            @test result2[1] == pt3
        end

        @testset "missing geometry column throws error" begin
            tbl = (
                name = ["a", "b"],
                value = [1, 2]
            )
            @test_throws Exception GO.get_geometries(tbl)
        end
    end

    @testset "Generic iterable" begin
        # Generator expression
        gen = (GI.Point(Float64(i), Float64(i)) for i in 1:3)
        result = GO.get_geometries(gen)
        @test result isa AbstractArray
        @test length(result) == 3
        @test GI.x(result[1]) == 1.0
        @test GI.x(result[2]) == 2.0
        @test GI.x(result[3]) == 3.0
    end

    @testset "Invalid input" begin
        @test_throws ArgumentError GO.get_geometries("not a geometry")
        @test_throws ArgumentError GO.get_geometries(42)
        @test_throws ArgumentError GO.get_geometries((a=1, b=2))  # NamedTuple without geometry
    end
end
