using Test

import GeoInterface as GI, GeometryOps as GO

open_rectangle = GI.Polygon([collect.([(0, 0), (10, 0), (10, 10), (0, 10)])])

# LibGEOS fails in GI.convert if the ring is not closed
@test_all_implementations "fix open rectangle" open_rectangle [GeoInterface, ArchGDAL, GeometryBasics] begin
    closed_rectangle = GO.fix(open_rectangle)
    @test GI.npoint(closed_rectangle) == GI.npoint(open_rectangle) + 1 # test that the rectangle is closed
    @test GI.getpoint(closed_rectangle.geom[1], 1) == GI.getpoint(closed_rectangle.geom[1], GI.npoint(closed_rectangle))
    @test all(GO.flatten(GI.PointTrait, closed_rectangle) .== GO.flatten(GI.PointTrait, GO.fix(closed_rectangle)))
end
