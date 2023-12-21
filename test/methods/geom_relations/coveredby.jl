# # Point and Geometry

# Same point -> covered by
@test GO.coveredby(pt1, pt1) == LG.coveredby(pt1, pt1) == true
# Different point -> not covered by
@test GO.coveredby(pt1, pt2) == LG.coveredby(pt1, pt2) == false
# Point on line endpoint -> covered by
@test GO.coveredby(pt1, l1) == LG.coveredby(pt1, l1) == true
# Point outside line -> not covered by
@test GO.coveredby(pt2, l1) == LG.coveredby(pt2, l1) == false
# Point on line segment -> covered by
@test GO.coveredby(pt3, l1) == LG.coveredby(pt3, l1) == true
# Point on line vertex between segments -> covered by
@test GO.coveredby(pt4, l1) == LG.coveredby(pt4, l1) == true
# line cannot be covered by a point -> not covered by
@test GO.coveredby(l1, pt3) == LG.coveredby(l1, pt3) == false
# Point on ring endpoint -> covered by
@test GO.coveredby(pt1, r1) == LG.coveredby(pt1, r1) == true
# Point outside ring -> isn't covered by
@test GO.coveredby(pt2, r1) == LG.coveredby(pt2, r1) == false
# Point on ring segment -> covered by
@test GO.coveredby(pt3, r1) == LG.coveredby(pt3, r1) == true
# Point on ring vertex between segments -> covered by
@test GO.coveredby(pt4, r1) == LG.coveredby(pt4, r1) == true
# Ring cannot be covered by a point -> isn't covered by
@test GO.coveredby(r1, pt3) == LG.coveredby(r1, pt3) == false
# Point on vertex of polygon --> covered
@test GO.coveredby(pt1, p1) == LG.coveredby(pt1, p1) == true
# Point outside of polygon's external ring -> not covered by
@test GO.coveredby(pt2, p1) == LG.coveredby(pt2, p1) == false
# Point on polygon's edge -> covered by
@test GO.coveredby(pt4, p1) == LG.coveredby(pt4, p1) == true
# Point inside of polygon -> covered by
@test GO.coveredby(pt5, p1) == LG.coveredby(pt5, p1) == true
# Point on hole edge -> covered by
@test GO.coveredby(pt6, p1) == LG.coveredby(pt6, p1) == true
# Point inside of polygon hole -> not covered by
@test GO.coveredby(pt7, p1) == LG.coveredby(pt7, p1) == false
# Polygon can't be covered by a polygon -> not covered by
@test GO.coveredby(p1, pt5) == LG.coveredby(p1, pt5) == false

# # Line and Geometry

# Same line -> covered by
@test GO.coveredby(l1, l1) == LG.coveredby(l1, l1) == true
# Line overlaps line edge and endpoint -> covered by
@test GO.coveredby(l2, l1) == LG.coveredby(l2, l1) == true
# Line overlaps with one edge and is outside of other edge -> isn't covered by
@test GO.coveredby(l3, l1) == LG.coveredby(l3, l1) == false
# Line segments both within other line segments -> covered by
@test GO.coveredby(l4, l1) == LG.coveredby(l4, l1) == true
# Line segments connect at endpoint -> isn't covered by
@test GO.coveredby(l5, l1) == LG.coveredby(l5, l1) == false
# Line segments don't touch -> isn't covered by
@test GO.coveredby(l6, l1) == LG.coveredby(l6, l1) == false
# Line segments cross -> isn't covered by
@test GO.coveredby(l7, l1) == LG.coveredby(l7, l1) == false
# Line segments cross and go over and out -> isn't covered by
@test GO.coveredby(l8, l1) == LG.coveredby(l8, l1) == false
# Line segments cross and overlap on endpoint -> isn't covered by
@test GO.coveredby(l9, l1) == LG.coveredby(l9, l1) == false