using Test
import GeoInterface as GI, GeometryOps as GO, LibGEOS as LG
using GeometryOps, GeoInterface, GeometryBasics, barycentric_coordinates, MeanValue

t1, p1 = GI.LineString([(0.,0.), (1.,0.), (0.,1.)]), (0.5, 0.5)
t2, p2 = GI.LineString([(0.,0.), (1.,0.), (0.,1.)]), (1., 1.)
t3, p3 = GI.LineString(GeometryBasics.Point2d[(0,0), (1,0), (0,1)]), (1., 1.)

@testset "Triangle coordinates" begin
    # Test that the barycentric coordinates for (0,0), (1,0), (0,1), and (0.5,0.5) are (0.5,0.25,0.25)
    @test all(barycentric_coordinates(MeanValue(), t1, p1) .== (0.5,0.25,0.25))
    # Test that the barycentric coordinates for (0,0), (1,0), (0,1), and (1,1) are (-1,1,1)
    @test all(barycentric_coordinates(MeanValue(), t2, p2) .≈ (-1,1,1))
    # Test that calculations with different number types (in this case Float64) are also accurate
    @test all(barycentric_coordinates(MeanValue(), Point2{Float64}, Point2{Float64}(1,1)) .≈ (-1,1,1))
end

@testset "Preserving return type" begin
    @test barycentric_coordinates(MeanValue(), Point2{BigFloat}[(0,0), (1,0), (0,1)], Point2{BigFloat}(1,1)) isa Vector{BigFloat}
    @test barycentric_coordinates(MeanValue(), Point2{BigFloat}[(0,0), (1,0), (0,1)], Point2f(1,1)) isa Vector{BigFloat} # keep the most precise type
    @test barycentric_coordinates(MeanValue(), Point2{Float64}[(0,0), (1,0), (0,1)], Point2{Float64}(1,1)) isa Vector{Float64}
    @test barycentric_coordinates(MeanValue(), Point2{Float32}[(0,0), (1,0), (0,1)], Point2{Float32}(1,1)) isa Vector{Float32}
end

@testset "Tests for helper methods" begin
    @testset "`t_value`" begin
        @test GO.t_value(Point2f(0,0), Point2f(1,0), 1, 1) == 0
        @test GO.t_value(Point2f(0, 1), Point2f(1, 0), 1, 2) == -0.5f0
    end
    @testset "`_det`" begin
        @test GO._det((1,0), (0,1)) == 1f0
        @test GO._det(Point2f(1, 2), Point2f(3, 4)) == -2.0f0
    end
end
