using Test

import GeoInterface as GI
import GeometryOps as GO

p1 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
# (p1, p2) -> one polygon inside of the other, sharing an edge
p2 = GI.Polygon([[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 0.0)]])
# (p1, p3) -> polygons outside of one another, sharing an edge
p3 = GI.Polygon([[(1.0, 0.0), (2.0, 0.0), (2.0, -1.0), (1.0, -1.0), (1.0, 0.0)]])
# (p1, p4) -> polygons are completely disjoint (no holes)
p4 = GI.Polygon([[(1.0, -1.0), (2.0, -1.0), (2.0, -2.0), (1.0, -2.0), (1.0, -1.0)]])
# (p1, p5) -> polygons cut through one another
p5 = GI.Polygon([[(1.0, -1.0), (2.0, -1.0), (2.0, 4.0), (1.0, 4.0), (1.0, -1.0)]])

mp1 = GI.MultiPolygon([p1])
mp2 = GI.MultiPolygon([p1, p1])
mp3 = GI.MultiPolygon([p1, p2, p3])
mp4 = GI.MultiPolygon([p1, p4])
mp5 = GI.MultiPolygon([p1, p5])
mp6 = GI.MultiPolygon([  # four interlocking polygons forming a hole
    [[(-5.0, 10.0), (-5.0, 15.0), (15.0, 15.0), (15.0, 10.0), (-5.0, 10.0)]],
    [[(-5.0, -5.0), (-5.0, 0.0), (15.0, 0.0), (15.0, -5.0), (-5.0, -5.0)]],
    [[(10.0, -5.0), (10.0, 15.0), (15.0, 15.0), (15.0, -5.0), (10.0, -5.0)]],
    [[(-5.0, -5.0), (-5.0, 15.0), (0.0, 15.0), (0.0, -5.0), (-5.0, -5.0)]],
])

@testset_implementations begin
    union_fixed_mp1 = GO.UnionIntersectingPolygons()(mp1)
    @test GI.npolygon(union_fixed_mp1) == 1
    @test GO.equals(GI.getpolygon(union_fixed_mp1, 1), p1)

    diff_fixed_mp1 = GO.DiffIntersectingPolygons()(mp1)
    @test GO.equals(diff_fixed_mp1, union_fixed_mp1)

    union_fixed_mp2 = GO.UnionIntersectingPolygons()(mp2)
    @test GI.npolygon(union_fixed_mp2) == 1
    @test GO.equals(GI.getpolygon(union_fixed_mp2, 1), p1)

    diff_fixed_mp2 = GO.DiffIntersectingPolygons()(mp2)
    @test GO.equals(diff_fixed_mp2, union_fixed_mp2)

    union_fixed_mp3 = GO.UnionIntersectingPolygons()(mp3)
    @test GI.npolygon(union_fixed_mp3) == 1
    @test all((GO.coveredby(p, union_fixed_mp3) for p in GI.getpolygon(mp3)))
    diff_polys = GO.difference(union_fixed_mp3, mp3; target = GI.PolygonTrait(), fix_multipoly = nothing)
    @test sum(GO.area, diff_polys; init = 0.0) == 0

    diff_fixed_mp3 = GO.DiffIntersectingPolygons()(mp3)
    @test GI.npolygon(diff_fixed_mp3) == 3
    @test all((GO.coveredby(p, union_fixed_mp3) for p in GI.getpolygon(diff_fixed_mp3)))

    union_fixed_mp4 = GO.UnionIntersectingPolygons()(mp4)
    @test GI.npolygon(union_fixed_mp4) == 2
    @test (GO.equals(GI.getpolygon(union_fixed_mp4, 1), p1) && GO.equals(GI.getpolygon(union_fixed_mp4, 2), p4)) ||
        (GO.equals(GI.getpolygon(union_fixed_mp4, 2), p1) && GO.equals(GI.getpolygon(union_fixed_mp4, 1), p4))

    diff_fixed_mp4 = GO.DiffIntersectingPolygons()(mp4)
    @test GO.equals(diff_fixed_mp4, union_fixed_mp4)

    union_fixed_mp5 = GO.UnionIntersectingPolygons()(mp5)
    @test GI.npolygon(union_fixed_mp5) == 1
    @test all((GO.coveredby(p, union_fixed_mp5) for p in GI.getpolygon(mp5)))
    diff_polys = GO.difference(union_fixed_mp5, mp5; target = GI.PolygonTrait(), fix_multipoly = nothing)
    @test sum(GO.area, diff_polys; init = 0.0) == 0

    diff_fixed_mp5 = GO.DiffIntersectingPolygons()(mp5)
    @test GI.npolygon(diff_fixed_mp5) == 3
    @test all((GO.coveredby(p, union_fixed_mp5) for p in GI.getpolygon(diff_fixed_mp5)))

    union_fixed_mp6 = GO.UnionIntersectingPolygons()(mp6)
    @test GI.npolygon(union_fixed_mp6) == 1
    @test GI.nhole(GI.getpolygon(union_fixed_mp6, 1)) == 1
    @test all((GO.coveredby(p, union_fixed_mp6) for p in GI.getpolygon(mp6)))
    diff_polys = GO.difference(union_fixed_mp6, mp6; target = GI.PolygonTrait(), fix_multipoly = nothing)
    @test sum(GO.area, diff_polys; init = 0.0) == 0

    diff_fixed_mp6 = GO.DiffIntersectingPolygons()(mp6)
    @test GI.npolygon(diff_fixed_mp6) == 4
    @test all((GO.coveredby(p, union_fixed_mp6) for p in GI.getpolygon(diff_fixed_mp6)))
end
