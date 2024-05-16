using Test

import GeoInterface as GI
import GeometryOps as GO

open_rectangle = GI.Wrappers.Polygon([collect.([(0, 0), (10, 0), (10, 10), (0, 10)])])

@test_all_integrations "Closed Ring correction" open_rectangle begin
    closed_rectangle = GO.ClosedRing()(open_rectangle)
    @test GI.npoint(closed_rectangle) == GI.npoint(open_rectangle) + 1 # test that the rectangle is closed
    @test GI.getpoint(closed_rectangle.geom[1], 1) == GI.getpoint(closed_rectangle.geom[1], GI.npoint(closed_rectangle))
    @test all(GO.flatten(GI.PointTrait, closed_rectangle) .== GO.flatten(GI.PointTrait, GO.ClosedRing()(closed_rectangle)))
end
