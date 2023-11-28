pt1 = LG.Point([0.0, 0.0])
pt2 = LG.Point([0.0, 0.1])
pt3 = LG.Point([1.0, 0.0])
pt4 = LG.Point([0.5, 1.0])
pt5 = LG.Point([0.2, 0.5])
pt6 = LG.Point([0.3, 0.55])
pt7 = LG.Point([0.6, 0.49])

l1 = LG.LineString([[0.0, 0.0], [0.0, 1.0]])
l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1]])
l3 = LG.LineString([[0.0, -1.0], [0.0, 0.5]])
l4 = LG.LineString([[0.0, -1.0], [0.0, 1.5]])
l5 =  LG.LineString([[0.0, -1.0], [0.0, 0.0]])
l6 = LG.LineString([[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0]])
l7 = LG.LineString([[0.6, 0.6], [0.8, 0.6]])
l8 = LG.LineString([[0.3, 0.55], [0.3, 0.65]])
l9 = LG.LineString([[0.2, 0.5], [0.3, 0.7]])
l10 = LG.LineString([[1.0, 0.0], [1.0, 1.0], [0.8, 0.4]])
l11 = LG.LineString([[1.0, 0.0], [1.0, 1.0], [0.7, 0.21]])
l12 = LG.LineString([[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0], [-1.0, 0.0]])

r1 = LG.LinearRing([[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0]])
r2 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1], [0.0, 0.0]])
r3 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [0.0, 0.2], [0.0, 0.0]])

p1 = LG.Polygon([[[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0]]])
p2 = LG.Polygon([
    [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]],
    [[0.2, 0.5], [0.3, 0.7], [0.4, 0.5], [0.2, 0.5]],
    [[0.5, 0.5], [0.8, 0.4], [0.7, 0.2], [0.5, 0.5]]
])

# Point and point
@test GO.within(pt1, pt1) == LG.within(pt1, pt1)
@test GO.within(pt1, pt2) == LG.within(pt1, pt2)

# Point and line
@test GO.within(pt1, l1) == LG.within(pt1, l1)
@test GO.within(pt2, l2) == LG.within(pt2, l2)
@test GO.within(pt2, l1) == LG.within(pt2, l1)
@test GO.within(pt3, l1) == LG.within(pt3, l1)
@test GO.within(pt1, l2) == LG.within(pt1, l2)
@test GO.within(pt2, l2) == LG.within(pt2, l2)
@test GO.within(pt3, l2) == LG.within(pt3, l2)
@test GO.within(pt1, l6) == GO.within(pt1, l6)

# Point and Ring
@test GO.within(pt1, r1) == LG.within(pt1, r1)
@test GO.within(pt2, r1) == LG.within(pt2, r1)
@test GO.within(pt3, r1) == LG.within(pt3, r1)
@test GO.within(pt4, r1) == LG.within(pt4, r1)

# Point and polygon
@test GO.within(pt1, p1) == LG.within(pt1, p1)
@test GO.within(pt2, p1) == LG.within(pt2, p1)
@test GO.within(pt3, p1) == LG.within(pt3, p1)
@test GO.within(pt4, p1) == LG.within(pt4, p1)

@test GO.within(pt1, p2) == LG.within(pt1, p2)
@test GO.within(pt2, p2) == LG.within(pt2, p2)
@test GO.within(pt3, p2) == LG.within(pt3, p2)
@test GO.within(pt4, p2) == LG.within(pt4, p2)
@test GO.within(pt5, p2) == LG.within(pt5, p2)
@test GO.within(pt6, p2) == LG.within(pt6, p2)
@test GO.within(pt7, p2) == LG.within(pt7, p2)

# Line and line
@test GO.within(l1, l1) == LG.within(l1, l1)
@test GO.within(l1, l2) == LG.within(l1, l2)
@test GO.within(l1, l3) == LG.within(l1, l3)
@test GO.within(l1, l4) == LG.within(l1, l4)
@test GO.within(l1, l5) == LG.within(l1, l5)
@test GO.within(l5, l1) == LG.within(l5, l1)

# Line and ring
@test GO.within(l6, r1) == LG.within(l6, r1)
@test GO.within(l2, r2) == LG.within(l2, r2)
@test GO.within(l2, r3) == LG.within(l2, r3)
@test GO.within(l12, r1) == LG.within(l12, r1)

# Line and polygon
@test GO.within(l6, p1) == LG.within(l6, p1)
@test GO.within(l1, p2) == LG.within(l1, p2)
@test GO.within(l2, p2) == LG.within(l2, p2)
@test GO.within(l3, p2) == LG.within(l3, p2)
@test GO.within(l7, p2) == LG.within(l7, p2)
@test GO.within(l8, p2) == LG.within(l8, p2)
@test GO.within(l9, p2) == LG.within(l9, p2)
@test GO.within(l10, p2) == LG.within(l10, p2)
@test GO.within(l11, p2) == LG.within(l11, p2)

# Ring and line
@test GO.within(r1, l6) == LG.within(r1, l6)
@test GO.within(r2, l2) == LG.within(r2, l2)
@test GO.within(r3, l2) == LG.within(r3, l2)
@test GO.within(r1, l12) == LG.within(r1, l12)

# Ring and Ring
@test GO.within(r1, r1) == LG.within(r1, r1)
@test GO.within(r1, r2) == LG.within(r1, r2)
@test GO.within(r1, r3) == LG.within(r1, r3)

# Ring and polygon

# Polygon in polygon
