# Same point --> intersects
@test GO.intersects(pt1, pt1) == LG.intersects(pt1, pt1) == true
# Different point --> don't intersect
@test GO.intersects(pt1, pt2) == LG.intersects(pt1, pt2) == false
# Point on line endpoint -> intersects
@test GO.intersects(pt1, l1) == GO.intersects(l1, pt1) == LG.intersects(l1, pt1) == true
# Point outside line -> don't intersect
@test GO.intersects(pt2, l1) == GO.intersects(l1, pt2) == LG.intersects(l1, pt2) == false
# Point on line segment -> intersects
@test GO.intersects(pt3, l1) == GO.intersects(l1, pt3) == LG.intersects(l1, pt3) == true
# Point on line vertex between segments -> intersects
@test GO.intersects(pt4, l1) == GO.intersects(l1, pt4) == LG.intersects(l1, pt4) == true
# Point on ring endpoint -> intersects
@test GO.intersects(pt1, r1) == GO.intersects(r1, pt1) == LG.intersects(r1, pt1) == true
# Point outside ring -> don't intersect
@test GO.intersects(pt2, r1) == GO.intersects(r1, pt2) == LG.intersects(r1, pt2) == false
# Point on ring segment -> intersects
@test GO.intersects(pt3, r1) == GO.intersects(r1, pt3) == LG.intersects(r1, pt3) == true
# Point on ring vertex between segments -> intersects
@test GO.intersects(pt4, r1) == GO.intersects(r1, pt4) == LG.intersects(r1, pt4) == true
# Point within hole formed by ring -> don't intersect
@test GO.intersects(pt5, r1) == GO.intersects(r1, pt5) == LG.intersects(r1, pt5) == false
# Point on vertex of polygon --> intersects
@test GO.intersects(pt1, p1) == GO.intersects(p1, pt1) == LG.intersects(p1, pt1) == true
# Point outside of polygon's external ring -> doesn't intersect
@test GO.intersects(pt2, p1) == GO.intersects(p1, pt2) == LG.intersects(p1, pt2) == false
# Point on polygon's edge -> intersects
@test GO.intersects(pt4, p1) == GO.intersects(p1, pt4) == LG.intersects(p1, pt4) == true
# Point inside of polygon -> intersects
@test GO.intersects(pt5, p1) == GO.intersects(p1, pt5) == LG.intersects(p1, pt5) == true
# Point on hole edge -> intersects
@test GO.intersects(pt6, p1) == GO.intersects(p1, pt6) == LG.intersects(p1, pt6) == true
# Point inside of polygon hole -> doesn't intersect
@test GO.intersects(pt7, p1)  == GO.intersects(p1, pt7) == LG.intersects(p1, pt7) == false

# # Line and Geometry

# Same line -> intersects
@test GO.intersects(l1, l1) == LG.intersects(l1, l1) == true
# Line overlaps line edge and endpoint -> intersects
@test GO.intersects(l1, l2) == LG.intersects(l1, l2) == true
# Line overlaps with one edge and is outside of other edge -> intersects
@test GO.intersects(l1, l3) == LG.intersects(l1, l3) == true
# Line segments both within other line segments -> intersects
@test GO.intersects(l1, l4) == LG.intersects(l1, l4) == true
# Line segments connect at endpoint -> intersects
@test GO.intersects(l1, l5) == LG.intersects(l1, l5) == true
# Line segments don't touch -> doesn't intersect
@test GO.intersects(l1, l6) == LG.intersects(l1, l6) == false
# Line segments cross -> intersects
@test GO.intersects(l1, l7) == LG.intersects(l1, l7) == true
# Line segments cross and go over and out -> intersects
@test GO.intersects(l1, l8) == LG.intersects(l1, l8) == true
# Line segments cross and overlap on endpoint -> intersects
@test GO.intersects(l1, l9) == LG.intersects(l1, l9) == true

# @testset "Lines/Rings" begin
#     # Line test intersects -----------------------------------------------------

#     # Test for parallel lines
    # l1 = GI.Line([(0.0, 0.0), (2.5, 0.0)])
    # l2 = GI.Line([(0.0, 1.0), (2.5, 1.0)])
    # @test !GO.intersects(l1, l2)
#     @test isnothing(GO.intersection(l1, l2))

#     # Test for non-parallel lines that don't intersect
#     l1 = GI.Line([(0.0, 0.0), (2.5, 0.0)])
#     l2 = GI.Line([(2.0, -3.0), (3.0, 0.0)])
#     @test !GO.intersects(l1, l2)
#     @test isnothing(GO.intersection(l1, l2))

#     # Test for lines only touching at endpoint
#     l1 = GI.Line([(0.0, 0.0), (2.5, 0.0)])
#     l2 = GI.Line([(2.0, -3.0), (2.5, 0.0)])
#     @test GO.intersects(l1, l2)
#     @test all(GO.intersection(l1, l2) .≈ (2.5, 0.0))

#     # Test for lines that intersect in the middle
#     l1 = GI.Line([(0.0, 0.0), (5.0, 5.0)])
#     l2 = GI.Line([(0.0, 5.0), (5.0, 0.0)])
#     @test GO.intersects(l1, l2)
#     @test all(GO.intersection(l1, l2) .≈ (2.5, 2.5))

#     # Line string test intersects ----------------------------------------------

#     # Single element line strings crossing over each other
#     l1 = LG.LineString([[5.5, 7.2], [11.2, 12.7]])
#     l2 = LG.LineString([[4.3, 13.3], [9.6, 8.1]])
#     @test GO.intersects(l1, l2)
#     go_inter = GO.intersection(l1, l2)
#     lg_inter = LG.intersection(l1, l2)
#     @test go_inter[1][1] .≈ GI.x(lg_inter)
#     @test go_inter[1][2] .≈ GI.y(lg_inter)

#     # Multi-element line strings crossing over on vertex
#     l1 = LG.LineString([[0.0, 0.0], [2.5, 0.0], [5.0, 0.0]])
#     l2 = LG.LineString([[2.0, -3.0], [3.0, 0.0], [4.0, 3.0]])
#     @test GO.intersects(l1, l2)
#     go_inter = GO.intersection(l1, l2)
#     @test length(go_inter) == 1
#     lg_inter = LG.intersection(l1, l2)
#     @test go_inter[1][1] .≈ GI.x(lg_inter)
#     @test go_inter[1][2] .≈ GI.y(lg_inter)

#     # Multi-element line strings crossing over with multiple intersections
#     l1 = LG.LineString([[0.0, -1.0], [1.0, 1.0], [2.0, -1.0], [3.0, 1.0]])
#     l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [3.0, 0.0]])
#     @test GO.intersects(l1, l2)
#     go_inter = GO.intersection(l1, l2)
#     @test length(go_inter) == 3
#     lg_inter = LG.intersection(l1, l2)
#     @test issetequal(
#         Set(go_inter),
#         Set(GO._tuple_point.(GI.getpoint(lg_inter)))
#     )

#     # Line strings far apart so extents don't overlap
#     l1 = LG.LineString([[100.0, 0.0], [101.0, 0.0], [103.0, 0.0]])
#     l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [3.0, 0.0]])
#     @test !GO.intersects(l1, l2)
#     @test isnothing(GO.intersection(l1, l2))

#     # Line strings close together that don't overlap
#     l1 = LG.LineString([[3.0, 0.25], [5.0, 0.25], [7.0, 0.25]])
#     l2 = LG.LineString([[0.0, 0.0], [5.0, 10.0], [10.0, 0.0]])
#     @test !GO.intersects(l1, l2)
#     @test isempty(GO.intersection(l1, l2))

#     # Closed linear ring with open line string
#     r1 = LG.LinearRing([[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]])
#     l2 = LG.LineString([[0.0, -2.0], [12.0, 10.0],])
#     @test GO.intersects(r1, l2)
#     go_inter = GO.intersection(r1, l2)
#     @test length(go_inter) == 2
#     lg_inter = LG.intersection(r1, l2)
#     @test issetequal(
#         Set(go_inter),
#         Set(GO._tuple_point.(GI.getpoint(lg_inter)))
#     )

#     # Closed linear ring with closed linear ring
#     r1 = LG.LinearRing([[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]])
#     r2 = LG.LineString([[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]])
#     @test GO.intersects(r1, r2)
#     go_inter = GO.intersection(r1, r2)
#     @test length(go_inter) == 2
#     lg_inter = LG.intersection(r1, r2)
#     @test issetequal(
#         Set(go_inter),
#         Set(GO._tuple_point.(GI.getpoint(lg_inter)))
#     )
# end

# @testset "Polygons" begin
#     # Two polygons that intersect
#     p1 = LG.Polygon([[[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]]])
#     p2 = LG.Polygon([[[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]]])
#     @test GO.intersects(p1, p2)
#     @test all(GO.intersection_points(p1, p2) .== [(6.5, 3.5), (6.5, -3.5)])

#     # Two polygons that don't intersect
#     p1 = LG.Polygon([[[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]]])
#     p2 = LG.Polygon([[[13.0, 0.0], [18.0, 5.0], [23.0, 0.0], [18.0, -5.0], [13.0, 0.0]]])
#     @test !GO.intersects(p1, p2)
#     @test isnothing(GO.intersection_points(p1, p2))

#     # Polygon that intersects with linestring
#     p1 = LG.Polygon([[[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]]])
#     l2 = LG.LineString([[0.0, 0.0], [10.0, 0.0]])
#     @test GO.intersects(p1, l2)
#     GO.intersection_points(p1, l2)
#     @test all(GO.intersection_points(p1, l2) .== [(0.0, 0.0), (10.0, 0.0)])

#     # Polygon with a hole, line through polygon and hole
#     p1 = LG.Polygon([
#         [[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]],
#         [[2.0, -1.0], [2.0, 1.0], [3.0, 1.0], [3.0, -1.0], [2.0, -1.0]]
#     ])
#     l2 = LG.LineString([[0.0, 0.0], [10.0, 0.0]])
#     @test GO.intersects(p1, l2)
#     @test all(GO.intersection_points(p1, l2) .== [(0.0, 0.0), (2.0, 0.0), (3.0, 0.0), (10.0, 0.0)])

#     # Polygon with a hole, line only within the hole
#     p1 = LG.Polygon([
#         [[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]],
#         [[2.0, -1.0], [2.0, 1.0], [3.0, 1.0], [3.0, -1.0], [2.0, -1.0]]
#     ])
#     l2 = LG.LineString([[2.25, 0.0], [2.75, 0.0]])
#     @test !GO.intersects(p1, l2)
#     @test isempty(GO.intersection_points(p1, l2))
# end

# @testset "MultiPolygons" begin
#     # TODO: Add these tests
#     # Multi-polygon and polygon that intersect

#     # Multi-polygon and polygon that don't intersect

#     # Multi-polygon that intersects with linestring
    
# end