# # Point and Geometry

# Same point -> contains
@test GO.contains(pt1, pt1) == LG.contains(pt1, pt1) == true
# Different point -> doesn't contain
@test GO.contains(pt1, pt2) == LG.contains(pt1, pt2) == false
# Point on line endpoint -> does not contain
@test GO.contains(l1, pt1) == LG.contains(l1, pt1) == false
# Point outside line -> does not contain
@test GO.contains(l1, pt2) == LG.contains(l1, pt2) == false
# Point on line segment -> contains
@test GO.contains(l1, pt3) == LG.contains(l1, pt3) == true
# Point on line vertex between segments -> contain
@test GO.contains(l1, pt4) == LG.contains(l1, pt4) == true
# Point cannot contain a line -> doesn't contain
@test GO.contains(pt3, l1) == LG.contains(pt3, l1) == false
# Point on ring endpoint -> contains
@test GO.contains(r1, pt1) == LG.contains(r1, pt1) == true
# Point outside ring -> does not contain
@test GO.contains(r1, pt2) == LG.contains(r1, pt2) == false
# Point on ring segment -> contains
@test GO.contains(r1, pt3) == LG.contains(r1, pt3) == true
# Point on ring vertex between segments -> contain
@test GO.contains(r1, pt4) == LG.contains(r1, pt4) == true
# Point cannot contain a ring -> doesn't contain
@test GO.contains(pt3, r1) == LG.contains(pt3, r1) == false
# Point on vertex of polygon --> doesn't contain
@test GO.contains(p1, pt1) == LG.contains(p1, pt1) == false
# Point outside of polygon's external ring -> doesn't contain
@test GO.contains(p1, pt2) == LG.contains(p1, pt2) == false
# Point on polygon's edge -> doesn't contain
@test GO.contains(p1, pt4) == LG.contains(p1, pt4) == false
# Point inside of polygon -> contains
@test GO.contains(p1, pt5) == LG.contains(p1, pt5) == true
# Point on hole edge -> doesn't contain
@test GO.contains(p1, pt6) == LG.contains(p1, pt6) == false
# Point inside of polygon hole -> doesn't contain
@test GO.contains(p1, pt7) == LG.contains(p1, pt7) == false
# Point cannot contain a polygon -> doesn't contain
@test GO.contains(pt5, p1) == LG.contains(pt5, p1) == false

# # Line and Geometry

# Same line -> contains
@test GO.contains(l1, l1) == LG.contains(l1, l1) == true
# Line overlaps line edge and endpoint -> contains
@test GO.contains(l1, l2) == LG.contains(l1, l2) == true
# Line overlaps with one edge and is outside of other edge -> doesn't contain
@test GO.contains(l1, l3) == LG.contains(l1, l3) == false
# Line segments both within other line segments -> contain
@test GO.contains(l1, l4) == LG.contains(l1, l4) == true
# Line segments connect at endpoint -> doesn't contain
@test GO.contains(l1, l5) == LG.contains(l1, l5) == false
# Line segments don't touch -> doesn't contain
@test GO.contains(l1, l6) == LG.contains(l1, l6) == false
# Line segments cross -> doesn't contain
@test GO.contains(l1, l7) == LG.contains(l1, l7) == false
# Line segments cross and go over and out -> doesn't contain
@test GO.contains(l1, l8) == LG.contains(l1, l8) == false
# Line segments cross and overlap on endpoint -> doesn't contain
@test GO.contains(l1, l9) == LG.contains(l1, l9) == false
# Line is within linear ring -> doesn't contain
@test GO.contains(l1, r1) == LG.contains(l1, r1) == false
# Line covers one edge of linera ring and has segment outside -> doesn't contain
@test GO.contains(l3, r1) == LG.contains(l3, r1) == false
# Line and linear ring are only connected at vertex -> doesn't contain
@test GO.contains(l5, r1) == LG.contains(l5, r1) == false
# Line and linear ring are disjoint -> doesn't contain
@test GO.contains(l6, r1) == LG.contains(l6, r1) == false
# Line crosses through two ring edges -> doesn't contain
@test GO.contains(l7, r1) == LG.contains(l7, r1) == false
# Line crosses through two ring edges and touches third edge -> doesn't contain
@test GO.contains(l8, r1) == LG.contains(l8, r1) == false
# Line is equal to linear ring -> contain
@test GO.contains(l10, r1) == LG.contains(l10, r1) == true
# Line covers linear ring and then has extra segment -> contain
@test GO.contains(l11, r1) == LG.contains(l11, r1) == true

# # Ring and Geometry

# Line is within linear ring -> contains
@test GO.contains(r1, l1) == LG.contains(r1, l1) == true
# Line covers one edge of linera ring and has segment outside -> doesn't contain
@test GO.contains(r1, l3) == LG.contains(r1, l3) == false
# Line and linear ring are only connected at vertex -> doesn't contain
@test GO.contains(r1, l5) == LG.contains(r1, l5) == false
# Line and linear ring are disjoint -> doesn't contain
@test GO.contains(r1, l6) == LG.contains(r1, l6) == false
# Line crosses through two ring edges -> doesn't contain
@test GO.contains(r1, l7) == LG.contains(r1, l7) == false
# Line crosses through two ring edges and touches third edge -> doesn't contain
@test GO.contains(r1, l8) == LG.contains(r1, l8) == false
# Line is equal to linear ring -> contain
@test GO.contains(r1, l10) == LG.contains(r1, l10) == true
# Line covers linear ring and then has extra segment -> doesn't contain
@test GO.contains(r1, l11) == LG.contains(r1, l11) == false



# p1 = LG.Point([0.0, 0.0])
# p2 = LG.Point([0.0, 0.1])
# p3 = LG.Point([1.0, 0.0])

# l1 = LG.LineString([[0.0, 1.0]])
# l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1]])

# # Point and point
# @test GO.contains(p1, p1) == LG.contains(p1, p1)
# @test GO.contains(p1, p2) == LG.contains(p1, p2)

# # Point and line
# @test GO.contains(l1, p1) == LG.contains(l1, p1)
# @test GO.contains(l1, p2) == LG.contains(l1, p2)
# @test GO.contains(l1, p3) == LG.contains(l1, p3)
# @test GO.contains(l2, p1) == LG.contains(l2, p1)
# @test GO.contains(l2, p2) == LG.contains(l2, p2)
# @test GO.contains(l2, p3) == LG.contains(l2, p3)

# # Line and line
# @test GO.contains(l1, l1) == LG.contains(l1, l1)
# @test GO.contains(l1, l2) == LG.contains(l1, l2)
# @test GO.contains(l1, l3) == LG.contains(l1, l3)
# @test GO.contains(l1, l4) == LG.contains(l1, l4)
# @test GO.contains(l1, l5) == LG.contains(l1, l5)

# @test GO.contais(l1, l1) == LG.contains(l1, l1)
# @test GO.contais(l2, l1) == LG.contains(l2, l1)
# @test GO.contais(l3, l1) == LG.contains(l3, l1)
# @test GO.contais(l4, l1) == LG.contains(l4, l1)
# @test GO.contais(l5, l1) == LG.contains(l5, l1)