@testset "Points/MultiPoints" begin
    p1 = LG.Point([0.0, 0.0])
    p2 = LG.Point([0.0, 1.0])
    # Same points
    @test GO.equals(p1, p1)
    @test GO.equals(p2, p2)
    # Different points
    @test !GO.equals(p1, p2)

    mp1 = LG.MultiPoint([[0.0, 1.0], [2.0, 2.0]])
    mp2 = LG.MultiPoint([[0.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    # Same points
    @test LG.equals(mp1, mp1)
    @test LG.equals(mp2, mp2)
    # Different points
    @test !LG.equals(mp1, mp2)
    @test !LG.equals(mp1, p1)
end

@testset "Lines/Rings" begin
    l1 = LG.LineString([[0.0, 0.0], [0.0, 10.0]])
    l2 = LG.LineString([[0.0, -10.0], [0.0, 20.0]])
    # Equal lines
    @test LG.equals(l1, l1)
    @test LG.equals(l2, l2)
    # Different lines
    @test !LG.equals(l1, l2) && !LG.equals(l2, l1)

    r1 = LG.LinearRing([[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]])
    r2 = LG.LinearRing([[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]])
    l3 = LG.LineString([[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]])
    # Equal rings
    @test GO.equals(r1, r1)
    @test GO.equals(r2, r2)
    # Different rings
    @test !GO.equals(r1, r2) && !GO.equals(r2, r1)
    # Equal linear ring and line string
    @test !GO.equals(r2, l3) # TODO: should these be equal?
end

@testset "Polygons/MultiPolygons" begin
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
    # Equal polygon
    @test GO.equals(p1, p1)
    @test GO.equals(p2, p2)
    # Different polygons
    @test !GO.equals(p1, p2)
    # Equal polygons with holes
    @test GO.equals(p3, p3)
    # Same exterior, different hole
    @test !GO.equals(p3, p4)
    # Same exterior and first hole, has an extra hole
    @test !GO.equals(p3, p5)

    p3 = GI.Polygon(
        [[
            [-53.57208251953125, 28.287451910503744],
            [-53.33038330078125, 28.29228897739706],
            [-53.34136962890625, 28.430052892335723],
            [-53.57208251953125, 28.287451910503744],
        ]]
    )
    # Complex polygon
    @test GO.equals(p3, p3)

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
    # Equal multipolygon
    @test GO.equals(m1, m1)
    # Equal multipolygon with different order
    @test GO.equals(m1, m2)
end