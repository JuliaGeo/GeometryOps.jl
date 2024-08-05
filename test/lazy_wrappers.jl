# Tests for lazy wrappers

# - Test that return type is inferred
# - Test that results are correct
# - Test that lazy evaluation is performed
# - Test for proper error handling
# - Test compatibility with different input types

using Test
using GeometryOps
import GeoInterface as GI, GeometryOps as GO

@testset "LazyClosedRing" begin
    # Helper function to create a simple LineString
    create_linestring(closed=false) = closed ? GI.LineString([
        (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0),
        
    ]) : GI.LineString([
        (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)
    ])

    @testset "Type inference" begin
        ls = create_linestring()
        wrapped = GO.LazyClosedRing(ls)
        @inferred GI.npoint(wrapped)
        @inferred GI.getpoint(wrapped, 1)
        @inferred collect(GI.getpoint(wrapped))
    end

    @testset "Correctness" begin
        ls = create_linestring()
        wrapped = GO.LazyClosedRing(ls)
        
        @test GI.npoint(wrapped) == 5
        @test GI.getpoint(wrapped, 1) == (0.0, 0.0)
        @test GI.getpoint(wrapped, 5) == (0.0, 0.0)
        @test collect(GI.getpoint(wrapped)) == [
            (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)
        ]
    end

    @testset "Lazy evaluation" begin
        ls = create_linestring()
        wrapped = GO.LazyClosedRing(ls)
        
        # Check that the original linestring is not modified
        @test GI.npoint(ls) == 4
        @test GI.npoint(wrapped) == 5
    end

    @testset "Error handling" begin
        ls = create_linestring()
        wrapped = GO.LazyClosedRing(ls)
        
        @test_throws BoundsError GI.getpoint(wrapped, 0)
        @test_throws BoundsError GI.getpoint(wrapped, 6)
    end

    @testset "Compatibility with different input types" begin
        # Test with a closed LineString
        closed_ls = create_linestring(true)
        wrapped_closed = GO.LazyClosedRing(closed_ls)
        @test GI.npoint(wrapped_closed) == 5
        @test collect(GI.getpoint(wrapped_closed)) == GI.getpoint(closed_ls)

        # Test with a 3D LineString
        ls_3d = GI.LineString([(x, y, 0.0) for (x, y) in GI.getpoint(create_linestring())])
        wrapped_3d = GO.LazyClosedRing(ls_3d)
        @test GI.is3d(wrapped_3d) == true
        @test GI.npoint(wrapped_3d) == 5
        @test GI.getpoint(wrapped_3d, 5) == (0.0, 0.0, 0.0)

        # Test with a measured LineString
        ls_measured = GI.LineString([GI.Point{false, true}(x, y, 0.0) for (x, y) in GI.getpoint(create_linestring())])
        wrapped_measured = GO.LazyClosedRing(ls_measured)
        @test GI.ismeasured(wrapped_measured) == true
        @test GI.npoint(wrapped_measured) == 5
    end
end