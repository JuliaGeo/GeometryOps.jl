using Test, GeometryOps
import GeoInterface as GI
import GeometryOps as GO

@testset "booleans" begin
    line1 = GI.LineString([[9.170356, 45.477985], [9.164434, 45.482551], [9.166644, 45.484003]])
    line2 = GI.LineString([[9.169356, 45.477985], [9.163434, 45.482551], [9.165644, 45.484003]])
    line3 = GI.LineString([
        (-111.544189453125, 24.186847428521244),
        (-110.687255859375, 24.966140159912975),
        (-110.4510498046875, 24.467150664739002),
        (-109.9951171875, 25.180087808990645)
    ])
	line4 = GI.LineString([
        (-111.4617919921875, 24.05148034322011),
        (-110.8795166015625, 24.681961205014595),
        (-110.841064453125, 24.14174098050432),
        (-109.97863769531249, 24.617057340809524)
    ])

    # @test isparallel(line1, line2) == true
	# @test isparallel(line3, line4) == false

	poly1 = GI.Polygon([[[0, 0], [1, 0], [1, 1], [0.5, 0.5], [0, 1], [0, 0]]])
	poly2 = GI.Polygon([[[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]]])

	@test isconcave(poly1) == true
	@test isconcave(poly2) == false

	l1 = GI.LineString([[0, 0], [1, 1], [1, 0], [0, 0]])
	l2 = GI.LineString([[0, 0], [1, 0], [1, 1], [0, 0]])

	@test isclockwise(l1) == true
	@test isclockwise(l2) == false

	l3 = GI.LineString([[0, 0], [3, 3], [4, 4]])
	p1 = GI.Point([1,1])

	l4 = GI.LineString([[0, 0], [3, 3]])
	p2 = GI.Point([0, 0])

	p3 = GI.Point([20, 20])
	l5 = GI.LineString([[0, 0], [3, 3], [38.32, 5.96]])

	@test GO.point_on_line(p2, l4; ignore_end_vertices=true) == false
	@test GO.point_on_line(p3, l5; ignore_end_vertices=true) == false
	@test GO.point_on_line(p1, l3) == true

	pt = (-77, 44)
	poly = GI.Polygon([[[-81, 41], [-81, 47], [-72, 47], [-72, 41], [-81, 41]]])

	@test point_in_polygon(pt, poly) == true

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
