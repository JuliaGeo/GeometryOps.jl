using Test
import GeoInterface as GI
import GeometryOps as GO 
import LibGEOS as LG
using ..TestHelpers

@testset "Basic example" begin
    points = tuple.(rand(100), rand(100))
    hull = GO.convex_hull(points)
    @test !GO.isconcave(hull) skip=true # TODO: fix
    @test !GO.isclockwise(GI.getexterior(hull)) # exterior should be ccw
    @test all(x -> GO.covers(hull, x), points) # but the orientation is right
    # Test against LibGEOS, by testing that all the points are the same
    # This is robust against winding order and starting at a different point, etc.
    @test isempty(
        setdiff(
            collect(GO.flatten(GO.tuples, GI.PointTrait, hull)), 
            collect(GO.flatten(GO.tuples, GI.PointTrait, GO.convex_hull(GO.GEOS(), points)))
        )
    )
end

@testset "Duplicated points" begin
    points = tuple.(rand(100), rand(100))
    @test_nowarn hull = GO.convex_hull(vcat(points, points))
    single_hull = GO.convex_hull(points)
    double_hull = GO.convex_hull(vcat(points, points))

    @test GO.equals(GI.getexterior(single_hull), GI.getexterior(double_hull))
    @test !GO.isconcave(double_hull) skip=true # TODO: fix
    @test !GO.isclockwise(GI.getexterior(double_hull)) # exterior should be ccw
    @test all(x -> GO.covers(single_hull, x), points)
    @test all(x -> GO.covers(double_hull, x), points)
    # Test against LibGEOS, by testing that all the points are the same
    # This is robust against winding order and starting at a different point, etc.
    @test isempty(
        setdiff(
            collect(GO.flatten(GO.tuples, GI.PointTrait, double_hull)), 
            collect(GO.flatten(GO.tuples, GI.PointTrait, GO.convex_hull(GO.GEOS(), points)))
        )
    )
end

# The rest of the tests are handled in DelaunayTriangulation.jl, this is simply testing
# that the methods work as expected.
