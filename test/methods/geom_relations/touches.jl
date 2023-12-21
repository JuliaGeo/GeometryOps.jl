# # Point and Geometry

# Same point -> doesn't touch
@test GO.touches(pt1, pt1) == LG.touches(pt1, pt1) == false
# Different point -> doesn't touch
@test GO.touches(pt1, pt2) == LG.touches(pt1, pt2) == false
# Point on line endpoint -> touches
@test GO.touches(pt1, l1) == GO.touches(l1, pt1) == LG.touches(l1, pt1) == true
# Point outside line -> doesn't touch
@test GO.touches(pt2, l1) == GO.touches(l1, pt2) == LG.touches(l1, pt2) == false
# Point on line segment -> doesn't touch
@test GO.touches(pt3, l1) == GO.touches(l1, pt3) == LG.touches(l1, pt3) == false
# Point on line vertex between segments -> doesn't touch
@test GO.touches(pt4, l1) == GO.touches(l1, pt4) == LG.touches(l1, pt4) == false
# Point on ring endpoint -> doesn't touch
@test GO.touches(pt1, r1) == GO.touches(r1, pt1) == LG.touches(r1, pt1) == false
# Point outside ring -> doesn't touch
@test GO.touches(pt2, r1) == GO.touches(r1, pt2) == LG.touches(r1, pt2) == false
# Point on ring segment -> doesn't touch
@test GO.touches(pt3, r1) == GO.touches(r1, pt3) == LG.touches(r1, pt3) == false
# Point on ring vertex between segments -> doesn't touch
@test GO.touches(pt4, r1) == GO.touches(r1, pt4) == LG.touches(r1, pt4) == false
# Point within hole formed by ring -> doesn't touch
@test GO.touches(pt5, r1) == GO.touches(r1, pt5) == LG.touches(r1, pt5) == false
# Point on vertex of polygon --> touches
@test GO.touches(pt1, p1) == GO.touches(p1, pt1) == LG.touches(p1, pt1) == true
# Point outside of polygon's external ring -> doesn't touch
@test GO.touches(pt2, p1) == GO.touches(p1, pt2) == LG.touches(p1, pt2) == false
# Point on polygon's edge -> touches
@test GO.touches(pt4, p1) == GO.touches(p1, pt4) == LG.touches(p1, pt4) == true
# Point inside of polygon -> doesn't touch
@test GO.touches(pt5, p1) == GO.touches(p1, pt5) == LG.touches(p1, pt5) == false
# Point on hole edge -> touches
@test GO.touches(pt6, p1) == GO.touches(p1, pt6) == LG.touches(p1, pt6) == true
# Point inside of polygon hole -> doesn't touch
@test GO.touches(pt7, p1)  == GO.touches(p1, pt7) == LG.touches(p1, pt7) == false

# # Line and Geometry

# Same line -> doesn't touch
@test GO.touches(l1, l1) == LG.touches(l1, l1) == false
# Line overlaps line edge and endpoint (nothing exterior) -> doesn't touch
@test GO.touches(l1, l2) == LG.touches(l1, l2) == false
# Line overlaps with one edge and is outside of other edge -> doesn't touch
@test GO.touches(l1, l3) == LG.touches(l1, l3) == false
# Line segments both within other line segments -> doesn't touch
@test GO.touches(l1, l4) == LG.touches(l1, l4) == false
# Line segments connect at endpoint -> touches
@test GO.touches(l1, l5) == LG.touches(l1, l5) == true
# Line segments don't touch -> doesn't touch
@test GO.touches(l1, l6) == LG.touches(l1, l6) == false
# Line segments cross -> doesn't touch
@test GO.touches(l1, l7) == LG.touches(l1, l7) == false
# Line segments cross and go over and out -> doesn't touch
@test GO.touches(l1, l8) == LG.touches(l1, l8) == false
# Line segments cross and overlap on endpoint -> doesn't touch
@test GO.touches(l1, l9) == LG.touches(l1, l9) == false




