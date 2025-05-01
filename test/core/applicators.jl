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