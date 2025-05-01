using Test
import GeometryOps as GO, GeometryOpsCore as GOCore, GeoInterface as GI

@testset "ApplyWithTrait" begin
    # Basic functionality test
    @testset "Basic functionality" begin
        # Create a simple function that uses the trait
        f = (trait, obj) -> (trait, obj)
        awt = GOCore.ApplyWithTrait(f)
        
        # Test with a point trait
        point = (1.0, 2.0)
        trait = GI.PointTrait()
        @test awt(trait, point) == (trait, point)
        
        # Test with a polygon trait
        poly = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)])])
        trait = GI.PolygonTrait()
        @test awt(trait, poly) == (trait, poly)
    end

    # Test rebuild method
    @testset "Rebuild method" begin
        f1 = (trait, obj) -> (trait, obj)
        f2 = (trait, obj) -> (trait, obj, "extra")
        
        awt = GOCore.ApplyWithTrait(f1)
        new_awt = GOCore.rebuild(awt, f2)
        
        point = (1.0, 2.0)
        trait = GI.PointTrait()
        @test new_awt(trait, point) == (trait, point, "extra")
    end

    # Test with apply
    @testset "Usage with apply" begin
        # Create a function that uses the trait to determine behavior
        f = (trait, obj) -> begin
            if trait isa GI.PointTrait
                (GI.x(obj) + 1, GI.y(obj) + 1)
            elseif trait isa GI.PolygonTrait
                GI.ngeom(obj)
            else
                error("Unexpected trait")
            end
        end
        
        awt = GOCore.ApplyWithTrait(f)
        
        # Test with a point
        point = (1.0, 2.0)
        result = GO.apply(awt, GI.PointTrait(), point)
        @test result == (2.0, 3.0)
        
        # Test with a polygon
        poly = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)])])
        result = GO.apply(awt, GI.PolygonTrait(), poly)
        @test result == 1  # One linear ring in the polygon
    end

    # Test with applyreduce
    @testset "Usage with applyreduce" begin
        # Create a function that uses the trait to determine behavior
        f = (trait, obj) -> begin
            if trait isa GI.PointTrait
                GI.x(obj) + GI.y(obj)
            elseif trait isa GI.PolygonTrait
                GI.ngeom(obj)
            else
                error("Unexpected trait")
            end
        end
        
        awt = GOCore.ApplyWithTrait(f)
        
        # Test with a point
        point = (1.0, 2.0)
        result = GO.applyreduce(awt, +, GI.PointTrait(), point)
        @test result == 3.0
        
        # Test with a polygon
        poly = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)])])
        result = GO.applyreduce(awt, +, GI.PolygonTrait(), poly)
        @test result == 1  # One linear ring in the polygon
    end

    # Test with keyword arguments
    @testset "Keyword arguments" begin
        f = (trait, obj; kw...) -> (trait, obj, kw)
        awt = GOCore.ApplyWithTrait(f)
        
        point = (1.0, 2.0)
        trait = GI.PointTrait()
        result = awt(trait, point; extra=1, more="test")
        @test result == (trait, point, pairs((;extra=1, more="test")))
    end
end

@testset "ApplyToPoint" begin
    @testset "Usage with apply" begin
        # Test WithXY with apply
        point = (1.0, 2.0)
        result = GO.apply(GOCore.WithXY((x, y) -> (x + 1, y + 1)), GI.PointTrait(), point)
        @test result == (2.0, 3.0)

        # Test WithXYZ with apply
        point = (1.0, 2.0, 3.0)
        result = GO.apply(GOCore.WithXYZ((x, y, z) -> (x + 1, y + 1, z + 1)), GI.PointTrait(), point)
        @test result == (2.0, 3.0, 4.0)

        # Test WithXYM with apply
        point = (1.0, 2.0, 3.0, 4.0)  # m value
        result = GO.apply(GOCore.WithXYM((x, y, m) -> (x + 1, y + 1, m + 1)), GI.PointTrait(), point)
        @test result == (2.0, 3.0, 5.0)

        # Test WithXYZM with apply
        point = (1.0, 2.0, 3.0, 4.0)  # x, y, z, m
        result = GO.apply(GOCore.WithXYZM((x, y, z, m) -> (x + 1, y + 1, z + 1, m + 1)), GI.PointTrait(), point)
        @test result == (2.0, 3.0, 4.0, 5.0)
    end

    @testset "Usage with applyreduce" begin
        # Test WithXY with applyreduce
        point = (1.0, 2.0)
        result = GO.applyreduce(GOCore.WithXY((x, y) -> x + y), +, GI.PointTrait(), point; init = 0.0)
        @test result == 3.0

        # Test WithXYZ with applyreduce
        point = (1.0, 2.0, 3.0)
        result = GO.applyreduce(GOCore.WithXYZ((x, y, z) -> x + y + z), +, GI.PointTrait(), point; init = 0.0)
        @test result == 6.0
    end
end

@testset "ApplyToArray" begin
    @testset "Usage with apply" begin
        # Test with array of geometries
        points = [(1.0, 2.0), (3.0, 4.0)]
        result = GO.apply(GOCore.WithXY((x, y) -> (x + 1, y + 1)), GI.PointTrait(), points)
        @test result == [(2.0, 3.0), (4.0, 5.0)]
    end

    @testset "Usage with applyreduce" begin
        # Test with array of geometries
        points = [(1.0, 2.0), (3.0, 4.0)]
        result = GO.applyreduce(GOCore.WithXY((x, y) -> x + y), +, GI.PointTrait(), points; init = 0.0)
        @test result == 10.0
    end
end

@testset "ApplyToFeatures" begin
    @testset "Usage with apply" begin
        # Test with feature collection
        features = [GI.Feature((1, 2)), GI.Feature((3, 4))]
        result = GO.apply(GOCore.WithXY((x, y) -> (x + 1, y + 1)), GI.PointTrait(), features)
        @test result == [GI.Feature((2, 3)), GI.Feature((4, 5))]
    end

    @testset "Usage with applyreduce" begin
        # Test with feature collection
        features = [GI.Feature((1, 2)), GI.Feature((3, 4))]
        result = GO.applyreduce(GOCore.WithXY((x, y) -> x + y), +, GI.PointTrait(), features; init = 0.0)
        @test result == 10.0
    end
end

@testset "ApplyToGeom" begin
    @testset "Usage with apply" begin
        # Test with polygon
        poly = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)])])
        result = GO.apply(GOCore.WithXY((x, y) -> (x + 1, y + 1)), GI.PointTrait(), poly)
        @test result == GI.Polygon([GI.LinearRing([(2, 3), (4, 5), (6, 7), (2, 3)])])
    end

    @testset "Usage with applyreduce" begin
        # Test with polygon
        poly = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)])])
        result = GO.applyreduce(GOCore.WithXY((x, y) -> x + y), +, GI.PointTrait(), poly; init = 0.0)
        @test result == 24.0  # Sum of all x+y values
    end
end

@testset "ApplyPointsToPolygon" begin
    @testset "Usage with apply" begin
        # Test with polygon
        poly = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)])])
        result = GO.apply(GOCore.WithXY((x, y) -> (x + 1, y + 1)), GI.PointTrait(), poly)
        @test result == GI.Polygon([GI.LinearRing([(2, 3), (4, 5), (6, 7), (2, 3)])])
    end

    @testset "Usage with applyreduce" begin
        # Test with polygon
        poly = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)])])
        result = GO.applyreduce(GOCore.WithXY((x, y) -> x + y), +, GI.PointTrait(), poly; init = 0.0)
        @test result == 24.0  # Sum of all x+y values
    end
end