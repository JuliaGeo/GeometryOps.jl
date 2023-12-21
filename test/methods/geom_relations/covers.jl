# # Point and Geometry

# Same point -> covers
@test GO.covers(pt1, pt1) == LG.covers(pt1, pt1) == true
# Different point -> doesn't cover
@test GO.covers(pt1, pt2) == LG.covers(pt1, pt2) == false
# Point on line endpoint -> covers
@test GO.covers(l1, pt1) == LG.covers(l1, pt1) == true
# Point outside line -> does not cover
@test GO.covers(l1, pt2) == LG.covers(l1, pt2) == false
# Point on line segment -> covers
@test GO.covers(l1, pt3) == LG.covers(l1, pt3) == true
# Point on line vertex between segments -> cover
@test GO.covers(l1, pt4) == LG.covers(l1, pt4) == true
# Point cannot cover a line -> doesn't cover
@test GO.covers(pt3, l1) == LG.covers(pt3, l1) == false
# Point on ring endpoint -> covers
@test GO.covers(r1, pt1) == LG.covers(r1, pt1) == true
# Point outside ring -> doesn't cover
@test GO.covers(r1, pt2) == LG.covers(r1, pt2) == false
# Point on ring segment -> covers
@test GO.covers(r1, pt3) == LG.covers(r1, pt3) == true
# Point on ring vertex between segments -> covers
@test GO.covers(r1, pt4) == LG.covers(r1, pt4) == true
# Point cannot cover a ring -> doesn't cover
@test GO.covers(pt3, r1) == LG.covers(pt3, r1) == false
# Point on vertex of polygon --> covers
@test GO.covers(p1, pt1) == LG.covers(p1, pt1) == true
# Point outside of polygon's external ring -> not covered
@test GO.covers(p1, pt2) == LG.covers(p1, pt2) == false
# Point on polygon's edge -> covers
@test GO.covers(p1, pt4) == LG.covers(p1, pt4) == true
# Point inside of polygon -> covers
@test GO.covers(p1, pt5) == LG.covers(p1, pt5) == true
# Point on hole edge -> covers
@test GO.covers(p1, pt6) == LG.covers(p1, pt6) == true
# Point inside of polygon hole -> not covered
@test GO.covers(p1, pt7) == LG.covers(p1, pt7) == false
# Point can't cover a polygon -> not covered
@test GO.covers(pt5, p1) == LG.covers(pt5, p1) == false

# # Line and Geometry

# Same line -> covers
@test GO.covers(l1, l1) == LG.covers(l1, l1) == true
# Line overlaps line edge and endpoint -> covers
@test GO.covers(l1, l2) == LG.covers(l1, l2) == true
# Line overlaps with one edge and is outside of other edge -> not covered
@test GO.covers(l1, l3) == LG.covers(l1, l3) == false
# Line segments both within other line segments -> covers
@test GO.covers(l1, l4) == LG.covers(l1, l4) == true
# Line segments connect at endpoint -> not covered
@test GO.covers(l1, l5) == LG.covers(l1, l5) == false
# Line segments don't touch -> not covered
@test GO.covers(l1, l6) == LG.covers(l1, l6) == false
# Line segments cross -> not covered
@test GO.covers(l1, l7) == LG.covers(l1, l7) == false
# Line segments cross and go over and out -> not covered
@test GO.covers(l1, l8) == LG.covers(l1, l8) == false
# Line segments cross and overlap on endpoint -> doesn't cover
@test GO.covers(l1, l9) == LG.covers(l1, l9) == false






# p1 = LG.Point([0.0, 0.0])
# p2 = LG.Point([0.0, 0.1])
# p3 = LG.Point([1.0, 0.0])

# l1 = LG.LineString([[0.0, 1.0]])
# l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1]])

# # Point and point
# @test GO.covers(p1, p1) == LG.covers(p1, p1)
# @test GO.covers(p1, p2) == LG.covers(p1, p2)

# # Point and line
# @test GO.covers(l1, p1) == LG.covers(l1, p1)
# @test GO.covers(l1, p2) == LG.covers(l1, p2)
# @test GO.covers(l1, p3) == LG.covers(l1, p3)
# @test GO.covers(l2, p1) == LG.covers(l2, p1)
# @test GO.covers(l2, p2) == LG.covers(l2, p2)
# @test GO.covers(l2, p3) == LG.covers(l2, p3)

# # Line and line
# @test GO.covers(l1, l1) == LG.covers(l1, l1)
# @test GO.covers(l2, l1) == LG.covers(l2, l1)
# @test GO.covers(l3, l1) == LG.covers(l3, l1)
# @test GO.covers(l4, l1) == LG.covers(l4, l1)
# @test GO.covers(l5, l1) == LG.covers(l5, l1)