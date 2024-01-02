@testset "intersects" begin
    @testset "Lines/Rings" begin
        # Line test intersects -----------------------------------------------------

        # Test for parallel lines
        l1 = GI.Line([(0.0, 0.0), (2.5, 0.0)])
        l2 = GI.Line([(0.0, 1.0), (2.5, 1.0)])
        @test !GO.intersects(l1, l2)
        @test isnothing(GO.intersection(l1, l2))

        # Test for non-parallel lines that don't intersect
        l1 = GI.Line([(0.0, 0.0), (2.5, 0.0)])
        l2 = GI.Line([(2.0, -3.0), (3.0, 0.0)])
        @test !GO.intersects(l1, l2)
        @test isnothing(GO.intersection(l1, l2))

        # Test for lines only touching at endpoint
        l1 = GI.Line([(0.0, 0.0), (2.5, 0.0)])
        l2 = GI.Line([(2.0, -3.0), (2.5, 0.0)])
        @test GO.intersects(l1, l2)
        @test all(GO.intersection(l1, l2) .≈ (2.5, 0.0))

        # Test for lines that intersect in the middle
        l1 = GI.Line([(0.0, 0.0), (5.0, 5.0)])
        l2 = GI.Line([(0.0, 5.0), (5.0, 0.0)])
        @test GO.intersects(l1, l2)
        @test all(GO.intersection(l1, l2) .≈ (2.5, 2.5))

        # Line string test intersects ----------------------------------------------

        # Single element line strings crossing over each other
        l1 = LG.LineString([[5.5, 7.2], [11.2, 12.7]])
        l2 = LG.LineString([[4.3, 13.3], [9.6, 8.1]])
        @test GO.intersects(l1, l2)
        go_inter = GO.intersection(l1, l2)
        lg_inter = LG.intersection(l1, l2)
        @test go_inter[1][1] .≈ GI.x(lg_inter)
        @test go_inter[1][2] .≈ GI.y(lg_inter)

        # Multi-element line strings crossing over on vertex
        l1 = LG.LineString([[0.0, 0.0], [2.5, 0.0], [5.0, 0.0]])
        l2 = LG.LineString([[2.0, -3.0], [3.0, 0.0], [4.0, 3.0]])
        @test GO.intersects(l1, l2)
        go_inter = GO.intersection(l1, l2)
        @test length(go_inter) == 1
        lg_inter = LG.intersection(l1, l2)
        @test go_inter[1][1] .≈ GI.x(lg_inter)
        @test go_inter[1][2] .≈ GI.y(lg_inter)

        # Multi-element line strings crossing over with multiple intersections
        l1 = LG.LineString([[0.0, -1.0], [1.0, 1.0], [2.0, -1.0], [3.0, 1.0]])
        l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [3.0, 0.0]])
        @test GO.intersects(l1, l2)
        go_inter = GO.intersection(l1, l2)
        @test length(go_inter) == 3
        lg_inter = LG.intersection(l1, l2)
        @test issetequal(
            Set(go_inter),
            Set(GO._tuple_point.(GI.getpoint(lg_inter)))
        )

        # Line strings far apart so extents don't overlap
        l1 = LG.LineString([[100.0, 0.0], [101.0, 0.0], [103.0, 0.0]])
        l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [3.0, 0.0]])
        @test !GO.intersects(l1, l2)
        @test isnothing(GO.intersection(l1, l2))

        # Line strings close together that don't overlap
        l1 = LG.LineString([[3.0, 0.25], [5.0, 0.25], [7.0, 0.25]])
        l2 = LG.LineString([[0.0, 0.0], [5.0, 10.0], [10.0, 0.0]])
        @test !GO.intersects(l1, l2)
        @test isempty(GO.intersection(l1, l2))

        # Closed linear ring with open line string
        r1 = LG.LinearRing([[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]])
        l2 = LG.LineString([[0.0, -2.0], [12.0, 10.0],])
        @test GO.intersects(r1, l2)
        go_inter = GO.intersection(r1, l2)
        @test length(go_inter) == 2
        lg_inter = LG.intersection(r1, l2)
        @test issetequal(
            Set(go_inter),
            Set(GO._tuple_point.(GI.getpoint(lg_inter)))
        )

        # Closed linear ring with closed linear ring
        r1 = LG.LinearRing([[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]])
        r2 = LG.LineString([[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]])
        @test GO.intersects(r1, r2)
        go_inter = GO.intersection(r1, r2)
        @test length(go_inter) == 2
        lg_inter = LG.intersection(r1, r2)
        @test issetequal(
            Set(go_inter),
            Set(GO._tuple_point.(GI.getpoint(lg_inter)))
        )
    end

    @testset "Polygons" begin
        # Two polygons that intersect
        p1 = LG.Polygon([[[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]]])
        p2 = LG.Polygon([[[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]]])
        @test GO.intersects(p1, p2)
        @test all(GO.intersection_points(p1, p2) .== [(6.5, 3.5), (6.5, -3.5)])

        # Two polygons that don't intersect
        p1 = LG.Polygon([[[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]]])
        p2 = LG.Polygon([[[13.0, 0.0], [18.0, 5.0], [23.0, 0.0], [18.0, -5.0], [13.0, 0.0]]])
        @test !GO.intersects(p1, p2)
        @test isnothing(GO.intersection_points(p1, p2))

        # Polygon that intersects with linestring
        p1 = LG.Polygon([[[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]]])
        l2 = LG.LineString([[0.0, 0.0], [10.0, 0.0]])
        @test GO.intersects(p1, l2)
        GO.intersection_points(p1, l2)
        @test all(GO.intersection_points(p1, l2) .== [(0.0, 0.0), (10.0, 0.0)])

        # Polygon with a hole, line through polygon and hole
        p1 = LG.Polygon([
            [[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]],
            [[2.0, -1.0], [2.0, 1.0], [3.0, 1.0], [3.0, -1.0], [2.0, -1.0]]
        ])
        l2 = LG.LineString([[0.0, 0.0], [10.0, 0.0]])
        @test GO.intersects(p1, l2)
        @test all(GO.intersection_points(p1, l2) .== [(0.0, 0.0), (2.0, 0.0), (3.0, 0.0), (10.0, 0.0)])

        # Polygon with a hole, line only within the hole
        p1 = LG.Polygon([
            [[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]],
            [[2.0, -1.0], [2.0, 1.0], [3.0, 1.0], [3.0, -1.0], [2.0, -1.0]]
        ])
        l2 = LG.LineString([[2.25, 0.0], [2.75, 0.0]])
        @test !GO.intersects(p1, l2)
        @test isempty(GO.intersection_points(p1, l2))
    end

    @testset "MultiPolygons" begin
        # TODO: Add these tests
        # Multi-polygon and polygon that intersect

        # Multi-polygon and polygon that don't intersect

        # Multi-polygon that intersects with linestring
        
    end
end
@testset "overlaps" begin
    @testset "Points/MultiPoints" begin
        p1 = LG.Point([0.0, 0.0])
        p2 = LG.Point([0.0, 1.0])
        # Two points can't overlap
        @test GO.overlaps(p1, p1) == LG.overlaps(p1, p2)
    
        mp1 = LG.MultiPoint([[0.0, 1.0], [4.0, 4.0]])
        mp2 = LG.MultiPoint([[0.0, 1.0], [2.0, 2.0]])
        mp3 = LG.MultiPoint([[0.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
        # No shared points, doesn't overlap
        @test GO.overlaps(p1, mp1) == LG.overlaps(p1, mp1)
        # One shared point, does overlap
        @test GO.overlaps(p2, mp1) == LG.overlaps(p2, mp1)
        # All shared points, doesn't overlap
        @test GO.overlaps(mp1, mp1) == LG.overlaps(mp1, mp1)
        # Not all shared points, overlaps
        @test GO.overlaps(mp1, mp2) == LG.overlaps(mp1, mp2)
        # One set of points entirely inside other set, doesn't overlap
        @test GO.overlaps(mp2, mp3) == LG.overlaps(mp2, mp3)
        # Not all points shared, overlaps
        @test GO.overlaps(mp1, mp3) == LG.overlaps(mp1, mp3)
    
        mp1 = LG.MultiPoint([
            [-36.05712890625, 26.480407161007275],
            [-35.7220458984375, 27.137368359795584],
            [-35.13427734375, 26.83387451505858],
            [-35.4638671875, 27.254629577800063],
            [-35.5462646484375, 26.86328062676624],
            [-35.3924560546875, 26.504988828743404],
        ])
        mp2 = GI.MultiPoint([
            [-35.4638671875, 27.254629577800063],
            [-35.5462646484375, 26.86328062676624],
            [-35.3924560546875, 26.504988828743404],
            [-35.2001953125, 26.12091815959972],
            [-34.9969482421875, 26.455820238459893],
        ])
        # Some shared points, overlaps
        @test GO.overlaps(mp1, mp2) == LG.overlaps(mp1, mp2)
        @test GO.overlaps(mp1, mp2) == GO.overlaps(mp2, mp1)
    end
    
    @testset "Lines/Rings" begin
        l1 = LG.LineString([[0.0, 0.0], [0.0, 10.0]])
        l2 = LG.LineString([[0.0, -10.0], [0.0, 20.0]])
        l3 = LG.LineString([[0.0, -10.0], [0.0, 3.0]])
        l4 = LG.LineString([[5.0, -5.0], [5.0, 5.0]])
        # Line can't overlap with itself
        @test GO.overlaps(l1, l1) == LG.overlaps(l1, l1)
        # Line completely within other line doesn't overlap
        @test GO.overlaps(l1, l2) == GO.overlaps(l2, l1) == LG.overlaps(l1, l2)
        # Overlapping lines
        @test GO.overlaps(l1, l3) == GO.overlaps(l3, l1) == LG.overlaps(l1, l3)
        # Lines that don't touch
        @test GO.overlaps(l1, l4) == LG.overlaps(l1, l4)
        # Linear rings that intersect but don't overlap
        r1 = LG.LinearRing([[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]])
        r2 = LG.LinearRing([[1.0, 1.0], [1.0, 6.0], [6.0, 6.0], [6.0, 1.0], [1.0, 1.0]])
        @test LG.overlaps(r1, r2) == LG.overlaps(r1, r2)
    end
    
    @testset "Polygons/MultiPolygons" begin
        p1 = LG.Polygon([[[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]]])
        p2 = LG.Polygon([
            [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
            [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]]
        ])
        # Test basic polygons that don't overlap
        @test GO.overlaps(p1, p2) == LG.overlaps(p1, p2)
        @test !GO.overlaps(p1, (1, 1))
        @test !GO.overlaps((1, 1), p2)
    
        p3 = LG.Polygon([[[1.0, 1.0], [1.0, 6.0], [6.0, 6.0], [6.0, 1.0], [1.0, 1.0]]])
        # Test basic polygons that overlap
        @test GO.overlaps(p1, p3) == LG.overlaps(p1, p3)
    
        p4 = LG.Polygon([[[20.0, 5.0], [20.0, 10.0], [18.0, 10.0], [18.0, 5.0], [20.0, 5.0]]])
        # Test one polygon within the other
        @test GO.overlaps(p2, p4) == GO.overlaps(p4, p2) == LG.overlaps(p2, p4)
    
        p5 = LG.Polygon(
            [[
                [-53.57208251953125, 28.287451910503744],
                [-53.33038330078125, 28.29228897739706],
                [-53.34136352890625, 28.430052892335723],
                [-53.57208251953125, 28.287451910503744],
            ]]
        )
        # Test equal polygons
        @test GO.overlaps(p5, p5) == LG.overlaps(p5, p5)
    
        # Test multipolygons
        m1 = LG.MultiPolygon([
            [[[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]]],
            [
                [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
                [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]]
            ]
        ])
        # Test polygon that overlaps with multipolygon
        @test GO.overlaps(m1, p3) == LG.overlaps(m1, p3)
        # Test polygon in hole of multipolygon, doesn't overlap
        @test GO.overlaps(m1, p4) == LG.overlaps(m1, p4)
    end
end
@testset "Contains, within, disjoint" begin
    poly3 = GI.Polygon([[(1, 1), (1, 10), (10, 10), (10, 1), (1, 1)]])
	poly4 = GI.Polygon([[(1, 1), (2, 2), (3, 2), (1, 1)]])
	line5 = GI.LineString([(1.0, 1.0), (2.0, 3.0), (2.0, 3.5)])

	line6 = GI.LineString([(1.0, 1.0), (1.0, 2.0), (1.0, 3.0), (1.0, 4.0)])
	poly5 = GI.Polygon([[(1.0, 1.0), (1.0, 20.0), (1.0, 3.0), (1.0, 4.0), (1.0, 1.0)]])
	line7 = GI.LineString([(1.0, 2.0), (1.0, 3.0), (1.0, 3.5)])

	@test GO.contains(poly3, poly4) == true
	@test GO.contains(poly3, line5) == true
	@test GO.contains(line6, (1, 2)) == true
	@test GO.contains(poly3, poly5) == false
	@test GO.contains(poly3 , line7) == false

	@test GO.within(poly4, poly3) == true
	@test GO.within(line5, poly3) == true
	@test GO.within(poly5, poly3) == false
	@test GO.within((1, 2), line6) == true
	@test GO.within(line7, poly3) == false

	poly6 = GI.Polygon([[(-11, -12), (-13, -12), (-13, -13), (-11, -13), (-11, -12)]])
	poly7 = GI.Polygon([[(-1, 2), (3, 2), (3, 3), (-1, 3), (-1, 2)]])
	poly8 = GI.Polygon([[(-1, 2), (-13, -12), (-13, -13), (-11, -13), (-1, 2)]])

	@test GO.disjoint(poly7, poly6) == true
	@test GO.disjoint(poly7, (1, 1)) == true
	@test GO.disjoint(poly7, GI.LineString([(0, 0), (12, 2), (12, 3), (12, 4)])) == true
	@test GO.disjoint(poly8, poly7) == false

	line8 = GI.LineString([(124.584961, -12.768946), (126.738281, -17.224758)])
	line9 = GI.LineString([(123.354492, -15.961329), (127.22168, -14.008696)])

    @test all(GO.intersection(line8, line9)[1] .≈ (125.583754, -14.835723))

	line10 = GI.LineString([
        (142.03125, -11.695273),
        (138.691406, -16.804541),
        (136.40625, -14.604847),
        (135.966797, -12.039321),
        (131.308594, -11.436955),
        (128.232422, -15.36895),
        (125.947266, -13.581921),
        (121.816406, -18.729502),
        (117.421875, -20.632784),
        (113.378906, -23.402765),
        (114.169922, -26.667096),
    ])
	line11 = GI.LineString([
        (117.861328, -15.029686),
        (122.124023, -24.886436),
        (132.583008, -22.309426),
        (132.890625, -7.754537),
    ])

	points = GO.intersection(line10, line11)
    @test all(points[1] .≈ (119.832884, -19.58857))
    @test all(points[2] .≈ (132.808697, -11.6309378))

	@test GO.crosses(GI.LineString([(-2, 2), (4, 2)]), line6) == true
	@test GO.crosses(GI.LineString([(0.5, 2.5), (1.0, 1.0)]), poly7) == true
	@test GO.crosses(GI.MultiPoint([(1, 2), (12, 12)]), GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])) == true
	@test GO.crosses(GI.MultiPoint([(1, 0), (12, 12)]), GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])) == false
	@test GO.crosses(GI.LineString([(-2, 2), (-4, 2)]), poly7) == false
end