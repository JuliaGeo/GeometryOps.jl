using Test

import GeoInterface as GI
using GeometryOpsTestHelpers

@testset "implementation macros" begin
    point = GI.Point(1.0, 2.0)
    line = GI.LineString([(0.0, 0.0), (1.0, 1.0)])

    @test_implementations [GI] begin
        GI.x($point) == 1.0 && GI.y($point) == 2.0
    end

    @test_implementations [GI] begin
        GI.npoint($line) == 2 && GI.x(GI.getpoint($line, 2)) == 1.0
    end

    @testset_implementations "tiny geometry costume party" [GI] begin
        @test GI.trait($point) isa GI.PointTrait
        @test GI.trait($line) isa GI.LineStringTrait
    end
end
