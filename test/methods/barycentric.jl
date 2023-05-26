@testset "Mean value coordinates" begin
    @testset "Preserving return type" begin
        @test barycentric_coordinates(MeanValue(), Point2{BigFloat}[(0,0), (1,0), (0,1)], Point2{BigFloat}(1,1)) isa Vector{BigFloat}
        @test barycentric_coordinates(MeanValue(), Point2{BigFloat}[(0,0), (1,0), (0,1)], Point2f(1,1)) isa Vector{BigFloat} # keep the most precise type
        @test barycentric_coordinates(MeanValue(), Point2{Float64}[(0,0), (1,0), (0,1)], Point2{Float64}(1,1)) isa Vector{Float64}
        @test barycentric_coordinates(MeanValue(), Point2{Float32}[(0,0), (1,0), (0,1)], Point2{Float32}(1,1)) isa Vector{Float32}
    end
    @testset "Triangle coordinates" begin
        # Test that the barycentric coordinates for (0,0), (1,0), (0,1), and (0.5,0.5) are (0.5,0.25,0.25)
        @test all(barycentric_coordinates(MeanValue(), Point2f[(0,0), (1,0), (0,1)], Point2f(0.5,0.5)) .== (0.5,0.25,0.25))
        # Test that the barycentric coordinates for (0,0), (1,0), (0,1), and (1,1) are (-1,1,1)
        @test all(barycentric_coordinates(MeanValue(), Point2f[(0,0), (1,0), (0,1)], Point2f(1,1)) .≈ (-1,1,1))
        # Test that calculations with different number types (in this case Float64) are also accurate
        @test all(barycentric_coordinates(MeanValue(), Point2{Float64}[(0,0), (1,0), (0,1)], Point2{Float64}(1,1)) .≈ (-1,1,1))
    end
end