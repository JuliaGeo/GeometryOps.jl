l1 = GI.Line([(5.0, -5.0), (5.0, 15.0)])
l2 = GI.Line([(-1.0, 8.0), (2.0, 11.0)])
l3 = GI.Line([(-1.0, 6.0), (11.0, 6.0)])
l4 = GI.Line([(-10.0, -10.0), (-1.0, -1.0)])
l5 = GI.Line([(-10.0, 5.0), (5.0, 5.0)])
r1 = GI.LinearRing([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)])
r2 = GI.LinearRing([(2.0, 2.0), (8.0, 2.0), (8.0, 4.0), (2.0, 4.0), (2.0, 2.0)])
p1 = GI.Polygon([r1])
p2 = GI.Polygon([r1, r2])
p3 = GI.Polygon([[(0.0, 0.0), (0.0, 5.0), (2.5, 7.5), (5.0, 5.0), (7.5, 7.5), (10.0, 5.0), (10.0, 0.0), (0.0, 0.0)]])

@test_all_implementations "Cut Polygon" (l1, l2, l3, l4, l5, r1, r2, p1, p2, p3), begin
    # Cut convex polygon
    cut_polys = GO.cut(p1, l1)
    @test all(GO.equals.(
        cut_polys,
        [
            GI.Polygon([[(0.0, 0.0), (5.0, 0.0), (5.0, 10.0), (0.0, 10.0), (0.0, 0.0)]]),
            GI.Polygon([[(5.0, 0.0), (10.0, 0.0), (10.0, 10.0), (5.0, 10.0), (5.0, 0.0)]]),
        ],
    ))

    cut_polys = GO.cut(p1, l2)
    @test all(GO.equals.(
        cut_polys,
        [
            GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (1.0, 10.0), (0.0, 9.0), (0.0, 0.0)]]),
            GI.Polygon([[(0.0, 9.0), (0.0, 10.0), (1.0, 10.0), (0.0, 9.0)]]),
        ],
    ))

    # Cut convex polygon with hole
    cut_polys = GO.cut(p2, l1)
    @test all(GO.equals.(
        cut_polys,
        [
            GI.Polygon([[(0.0, 0.0), (5.0, 0.0), (5.0, 2.0), (2.0, 2.0), (2.0, 4.0), (5.0, 4.0), (5.0, 10.0), (0.0, 10.0), (0.0, 0.0)]]),
            GI.Polygon([[(5.0, 0.0), (10.0, 0.0), (10.0, 10.0), (5.0, 10.0), (5.0, 4.0), (8.0, 4.0), (8.0, 2.0), (5.0, 2.0), (5.0, 0.0)]]),
        ],
    ))
    cut_polys = GO.cut(p2, l3)
    @test all(GO.equals.(
        cut_polys,
        [
            GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0), (0.0, 0.0)], r2]),
            GI.Polygon([[(0.0, 6.0), (10.0, 6.0), (10.0, 10.0), (0.0, 10.0), (0.0, 6.0)]]),
        ],
    ))

    # Cut concave polygon into three pieces
    cut_polys = GO.cut(p3, l3)
    @test all(GO.equals.(
        cut_polys,
        [
            GI.Polygon([[(0.0, 0.0), (0.0, 5.0), (1.0, 6.0), (4.0, 6.0), (5.0, 5.0), (6.0, 6.0), (9.0, 6.0), (10.0, 5.0), (10.0, 0.0), (0.0, 0.0)]]),
            GI.Polygon([[(1.0, 6.0), (2.5, 7.5), (4.0, 6.0), (1.0, 6.0)]]),
            GI.Polygon([[(6.0, 6.0), (7.5, 7.5), (9.0, 6.0), (6.0, 6.0)]]),
        ],
    ))

    # Line doesn't cut through polygon
    cut_polys = GO.cut(p1, l4)
    @test all(GO.equals.(cut_polys, [p1]))
    cut_polys = GO.cut(p1, l5)
    @test all(GO.equals.(cut_polys, [p1]))
end
