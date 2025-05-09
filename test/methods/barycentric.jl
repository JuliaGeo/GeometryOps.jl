using Test
import GeoInterface as GI, GeometryOps as GO, LibGEOS as LG
using GeometryOps, GeoInterface, GeometryBasics
import GeometryOps: barycentric_coordinates, MeanValue

t1, p1 = GI.LineString([(0.,0.), (1.,0.), (0.,1.)]), (0.5, 0.5)
t2, p2 = GI.LineString([(0.,0.), (1.,0.), (0.,1.)]), (1., 1.)
t3, p3 = GI.LineString(GeometryBasics.Point2d[(0,0), (1,0), (0,1)]), (1., 1.)

@testset "Triangle coordinates" begin
    # Test that the barycentric coordinates for (0,0), (1,0), (0,1), and (0.5,0.5) are (0.5,0.25,0.25)
    @test all(barycentric_coordinates(MeanValue(), t1, p1) .== (0.5,0.25,0.25))
    # Test that the barycentric coordinates for (0,0), (1,0), (0,1), and (1,1) are (-1,1,1)
    @test all(barycentric_coordinates(MeanValue(), t2, p2) .≈ (-1,1,1))
end

@testset "Preserving return type" begin
    @test eltype(barycentric_coordinates(MeanValue(), GI.LinearRing(Point2{BigFloat}[(0,0), (1,0), (0,1)]), Point2{BigFloat}(1,1))) <: BigFloat
    @test_broken eltype(barycentric_coordinates(MeanValue(), GI.LinearRing(Point2{BigFloat}[(0,0), (1,0), (0,1)]), Point2f(1,1))) <: BigFloat # keep the most precise type
    @test eltype(barycentric_coordinates(MeanValue(), GI.LinearRing(Point2{Float64}[(0,0), (1,0), (0,1)]), Point2{Float64}(1,1))) <: Float64
    @test eltype(barycentric_coordinates(MeanValue(), GI.LinearRing(Point2{Float32}[(0,0), (1,0), (0,1)]), Point2{Float32}(1,1))) <: Float32
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
