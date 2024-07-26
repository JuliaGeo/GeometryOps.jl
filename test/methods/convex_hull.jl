using Test
import GeoInterface as GI, GeometryOps as GO, LibGEOS as LG

@testset "Hull winding order" begin
    points = tuple.(rand(100), rand(100))
    hull = GO.convex_hull(points)
    @test !GO.isconcave(hull) skip=true # TODO: fix
    @test !GO.isclockwise(GI.getexterior(hull)) # exterior should be ccw
end

@testset "Duplicated points" begin
    points = tuple.(rand(100), rand(100))
    @test_nowarn hull = GO.convex_hull(vcat(points, points))
    single_hull = GO.convex_hull(points)
    double_hull = GO.convex_hull(vcat(points, points))

    @test GO.equals(GI.getexterior(single_hull), GI.getexterior(double_hull))
    @test !GO.isconcave(double_hull) skip=true # TODO: fix
    @test !GO.isclockwise(GI.getexterior(double_hull)) # exterior should be ccw
end

# The rest of the tests are handled in DelaunayTriangulation.jl, this is simply testing
# that the methods work as expected.