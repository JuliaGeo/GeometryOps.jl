# # Point and Geometry

# Same point -> doesn't cross
@test GO.crosses(pt1, pt1) == LG.crosses(pt1, pt1) == false
# Different point -> doesn't cross
@test GO.crosses(pt1, pt2) == LG.crosses(pt1, pt2) == false
# Point cannot cross line -> doesn't cross
@test GO.crosses(pt3, l1) == LG.crosses(pt3, l1) == false
# Line cannot cross point -> doesn't cross
@test GO.crosses(l1, pt3) == LG.crosses(l1, pt3) == false
# Point cannot cross ring -> doesn't cross
@test GO.crosses(pt3, r1) == LG.crosses(pt3, r1) == false
# Ring cannot cross point -> doesn't cross
@test GO.crosses(r1, pt3) == LG.crosses(r1, pt3) == false
# Point cannot cross polygon -> doesn't cross
@test GO.crosses(pt3, p1) == LG.crosses(pt3, p1) == false
# Polygon cannot cross point -> doesn't cross
@test GO.crosses(p1, pt3) == LG.crosses(p1, pt3) == false

# # Line and Geometry

# Same line -> doesn't cross
@test GO.crosses(l1, l1) == LG.crosses(l1, l1) == false
# Line overlaps line edge and endpoint -> doesn't cross
@test GO.crosses(l1, l2) == LG.crosses(l1, l2) == false
# Line overlaps with one edge and is outside of other edge -> not covered
@test GO.crosses(l1, l3) == LG.crosses(l1, l3) == false
# Line segments both within other line segments -> doesn't cross
@test GO.crosses(l1, l4) == LG.crosses(l1, l4) == false
# Line segments connect at endpoint -> doesn't cross
@test GO.crosses(l1, l5) == LG.crosses(l1, l5) == false
# Line segments don't touch -> doesn't cross
@test GO.crosses(l1, l6) == LG.crosses(l1, l6) == false
# Line segments cross -> crosses
@test GO.crosses(l1, l7) == LG.crosses(l1, l7) == true
# Line segments cross and go over and out -> doesn't cross
@test GO.crosses(l1, l8) == LG.crosses(l1, l8) == false
# Line segments cross and overlap on endpoint -> crosses
@test GO.crosses(l1, l9) == LG.crosses(l1, l9) == true