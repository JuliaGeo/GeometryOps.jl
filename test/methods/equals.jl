@testset "Points/MultiPoints" begin
    p1 = LG.Point([0.0, 0.0])
    p2 = LG.Point([0.0, 1.0])
    mp1 = LG.MultiPoint([[0.0, 1.0], [2.0, 2.0]])
    mp2 = LG.MultiPoint([[0.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    mp3 = LG.MultiPoint([p2])

    @test_all_integrations (p1, p2, mp1, mp2, mp3) begin
        # Same points
        @test GO.equals(p1, p1) == LG.equals(p1, p1)
        @test GO.equals(p2, p2) == LG.equals(p2, p2)
        # Different points
        @test GO.equals(p1, p2) == LG.equals(p1, p2)

        # Same points
        @test GO.equals(mp1, mp1) == LG.equals(mp1, mp1)
        @test GO.equals(mp2, mp2) == LG.equals(mp2, mp2)
        # Different points
        @test GO.equals(mp1, mp2) == LG.equals(mp1, mp2)
        @test GO.equals(mp1, p1) == LG.equals(mp1, p1)
        # Point and multipoint
        @test GO.equals(p2, mp3) == LG.equals(p2, mp3)
    end
end

@testset "Lines/Rings" begin
    l1 = LG.LineString([[0.0, 0.0], [0.0, 10.0]])
    l2 = LG.LineString([[0.0, -10.0], [0.0, 20.0]])
    l3 = LG.LineString([[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]])
    r1 = LG.LinearRing([[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]])
    r2 = LG.LinearRing([[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]])
    r3 = GI.LinearRing([[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0]])

    @test_all_integrations (l1, l2, l3, r1, r2, r3) begin
        # Equal lines
        @test GO.equals(l1, l1) == LG.equals(l1, l1)
        @test GO.equals(l2, l2) == LG.equals(l2, l2)
        # Different lines
        @test GO.equals(l1, l2) == GO.equals(l2, l1) == LG.equals(l1, l2)

        # Equal rings
        @test GO.equals(r1, r1) == LG.equals(r1, r1)
        @test GO.equals(r2, r2) == LG.equals(r2, r2)
        # Test equal rings without closing point
        @test GO.equals(r2, r3)
        @test GO.equals(r3, l3)
        # Different rings
        @test GO.equals(r1, r2) == GO.equals(r2, r1) == LG.equals(r1, r2)
        # Equal linear ring and line string
        @test GO.equals(r2, l3) == LG.equals(r2, l3)
        # Equal line string and line
        @test GO.equals(l1, GI.Line([(0.0, 0.0), (0.0, 10.0)]))
    end
end

@testset "Polygons/MultiPolygons" begin
    pt1 = LG.Point([0.0, 0.0])
    r1 = GI.LinearRing([(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)])
    p1 = GI.Polygon([[(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)]])
    p2 = GI.Polygon([[(1, 1), (1, 6), (6, 6), (6, 1), (1, 1)]])
    p3 = LG.Polygon(
        [
            [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
            [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]]
        ]
    )
    p4 = LG.Polygon(
        [
            [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
            [[16.0, 1.0], [16.0, 11.0], [25.0, 11.0], [25.0, 1.0], [16.0, 1.0]]
        ]
    )
    p5 = LG.Polygon(
        [
            [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
            [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]],
            [[11.0, 1.0], [11.0, 2.0], [12.0, 2.0], [12.0, 1.0], [11.0, 1.0]]
        ]
    )
    p6 = GI.Polygon([[(6, 6), (6, 1), (1, 1), (1, 6), (6, 6)]])
    p7 = GI.Polygon([[(6, 6), (1, 6), (1, 1), (6, 1), (6, 6)]])
    p8 = GI.Polygon([[(6, 6), (1, 6), (1, 1), (6, 1)]])
    p9 = LG.Polygon(
        [[
            [-53.57208251953125, 28.287451910503744],
            [-53.33038330078125, 28.29228897739706],
            [-53.34136962890625, 28.430052892335723],
            [-53.57208251953125, 28.287451910503744],
        ]]
    )
    m1 = LG.MultiPolygon([
        [[[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]]],
        [
            [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
            [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]]
        ]
    ])
    m2 = LG.MultiPolygon([
        [
            [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
            [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]]
        ],
        [[[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]]]
    ])
    m3 = LG.MultiPolygon([p3])

    @test_all_integrations "Polygons" (pt1, r1, p1, p2, p3, p4, p5, p6, p7, p8) begin
        # Point and polygon aren't equal
        GO.equals(pt1, p1) == LG.equals(pt1, p1)
        # Linear ring and polygon aren't equal
        @test GO.equals(r1, p1) == LG.equals(r1, p1)
        # Equal polygon
        @test GO.equals(p1, p1) == LG.equals(p1, p1)
        @test GO.equals(p2, p2) == LG.equals(p2, p2)
        # Equal but offset polygons
        @test GO.equals(p2, p6) == LG.equals(p2, p6)
        # Equal but opposite winding orders
        @test GO.equals(p2, p7) == LG.equals(p2, p7)
        # Equal but without closing point (implied)
        @test GO.equals(p7, p8) 
        # Different polygons
        @test GO.equals(p1, p2) == LG.equals(p1, p2)
        # Equal polygons with holes
        @test GO.equals(p3, p3) == LG.equals(p3, p3)
        # Same exterior, different hole
        @test GO.equals(p3, p4) == LG.equals(p3, p4)
        # Same exterior and first hole, has an extra hole
        @test GO.equals(p3, p5) == LG.equals(p3, p5)
        # Complex polygon
        @test GO.equals(p9, p9) == LG.equals(p9, p9)
    end

    @test_all_integrations "MultiPolygons" (p1, m1, m2, m3) begin
        # Equal multipolygon
        @test GO.equals(m1, m1) == LG.equals(m1, m1)
        # Equal multipolygon with different order
        @test GO.equals(m2, m2) == LG.equals(m2, m2)
        # Equal polygon to multipolygon
        @test GO.equals(p1, m3) == LG.equals(p1, m3)
    end
end
