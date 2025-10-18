using Test
import ArchGDAL as AG
import GeoInterface as GI
import GeometryBasics as GB
import GeometryOps as GO
using .TestHelpers

open_rectangle = GI.Polygon([collect.([(0, 0), (10, 0), (10, 10), (0, 10)])])

# LibGEOS fails in GI.convert if the ring is not closed
@testset_implementations "Closed Ring correction" [GI, AG, GB] begin
    closed_rectangle = GO.ClosedRing()($open_rectangle)
    @test GI.npoint(closed_rectangle) == GI.npoint($open_rectangle) + 1 # test that the rectangle is closed
    @test GI.getpoint(closed_rectangle.geom[1], 1) == GI.getpoint(closed_rectangle.geom[1], GI.npoint(closed_rectangle))
    @test all(GO.flatten(GI.PointTrait, closed_rectangle) .== GO.flatten(GI.PointTrait, GO.ClosedRing()(closed_rectangle)))
end
