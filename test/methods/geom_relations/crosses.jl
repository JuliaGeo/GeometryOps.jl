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

#=
Extra tests since this uses different code due to curve-curve crosses definition
requiring meeting in a point and having lines cross one another
=#
h = LG.LineString([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
v1 = LG.LineString([[0.0, 1.0], [0.0, 0.0], [0.0, -1.0]])
v2 = LG.LineString([[2.0, 1.0], [2.0, 0.0], [2.0, -1.0]])
v3 = LG.LineString([[1.0, 1.0], [1.0, 0.0], [1.0, -1.0]])
v4 = LG.LineString([[1.5, 1.0], [1.5, 0.0], [1.5, -1.0]])
v5 = LG.LineString([[1.5, 1.0], [1.5, -1.0]])
d1 = LG.LineString([[-1.0, -1.0], [1.0, 1.0]])
d2 = LG.LineString([[-1.0, 1.0], [0.0, 0.0], [1.0, -1.0]])
b1 = LG.LineString([[0.0, 1.0], [1.5, 0.0], [0.0, -1.0]])
b2 = LG.LineString([[1.0, 1.0], [1.5, 0.0], [2.0, 1.0]])
b3 = LG.LineString([[0.0, 1.0], [1.0, 0.0], [2.0, 1.0]])
b4 = LG.LineString([[1.0, -0.5], [0.0, 0.0], [1.0, 0.0]])
b5 = LG.LineString([[-1.0, 0.0], [0.0, 0.0], [1.0, 1.0]])
b6 = LG.LineString([[-1.0, 0.0], [0.0, 0.0], [1.0, -1.0]])
b7 = LG.LineString([[1.0, 0.0], [0.0, 0.0], [-1.0, -1.0]])
# Crosses through line starting endpoint -> doesn't cross
@test GO.crosses(h, v1) == LG.crosses(h, v1) == false
# Crosses through line ending endpoint -> doesn't cross
@test GO.crosses(h, v2) == LG.crosses(h, v2) == false
# Crosses through line middle vertex -> crosses
@test GO.crosses(h, v3) == LG.crosses(h, v3) == true
# Crosses through line edge at vertex -> crosses
@test GO.crosses(h, v4) == LG.crosses(h, v4) == true
# Crosses through line edge -> crosses
@test GO.crosses(h, v5) == LG.crosses(h, v5) == true
# Crosses through line edge -> crosses
@test GO.crosses(v5, h) == LG.crosses(v5, h) == true
# Line bounces off of vertical curve on edge -> doesn't cross
@test GO.crosses(b1, v5) == GO.crosses(v5, b1) == LG.crosses(b1, v5) == false
# Line bounces off of horizontal curve on edge --> doesn't cross
@test GO.crosses(b2, h) == GO.crosses(h, b2) == LG.crosses(b2, h) == false
# Line bounces off of horizontal curve on vertex --> doesn't cross
@test GO.crosses(b3, h) == GO.crosses(h, b3) == LG.crosses(b3, h) == false
# Diagonal lines pass through one another --> crosses
@test GO.crosses(d1, d2) == GO.crosses(d2, d1) == LG.crosses(d1, d2) == true
# Curve bounces off of diagonal line -> doesn't cross
@test GO.crosses(d1, b4) == GO.crosses(b4, d1) == LG.crosses(d1, b4) == false
# Lines with parallel segments cross -> cross
@test GO.crosses(b5, b7) == GO.crosses(b7, b5) == LG.crosses(b7, b5) == true
# Lines with parallel segments bounce -> crosses
@test GO.crosses(b6, b7) == GO.crosses(b7, b6) == LG.crosses(b7, b6) == true

