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

r1 = LG.LinearRing([[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0]])

p1 = LG.Polygon([[[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0]]])
p2 = LG.Polygon([
    [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]],
    [[0.2, 0.5], [0.3, 0.7], [0.4, 0.5], [0.2, 0.5]],
    [[0.5, 0.5], [0.8, 0.4], [0.7, 0.2], [0.5, 0.5]]
])

# Point and point
@test GO.disjoint(pt1, pt1) == LG.disjoint(pt1, pt1)
@test GO.disjoint(pt1, pt2) == LG.disjoint(pt1, pt2)

# Point and line
@test GO.disjoint(pt1, l1) == LG.disjoint(pt1, l1)
@test GO.disjoint(pt2, l1) == LG.disjoint(pt2, l1)
@test GO.disjoint(pt3, l1) == LG.disjoint(pt3, l1)
@test GO.disjoint(pt1, l2) == LG.disjoint(pt1, l2)
@test GO.disjoint(pt2, l2) == LG.disjoint(pt2, l2)
@test GO.disjoint(pt3, l2) == LG.disjoint(pt3, l2)
@test GO.within(pt1, l6) == GO.within(pt1, l6)

# Point and Ring
@test GO.disjoint(pt1, r1) == LG.disjoint(pt1, r1)
@test GO.disjoint(pt2, r1) == LG.disjoint(pt2, r1)
@test GO.disjoint(pt3, r1) == LG.disjoint(pt3, r1)
@test GO.disjoint(pt4, r1) == LG.disjoint(pt4, r1)

# Point and polygon
@test GO.disjoint(pt1, p1) == LG.disjoint(pt1, p1)
@test GO.disjoint(pt2, p1) == LG.disjoint(pt2, p1)
@test GO.disjoint(pt3, p1) == LG.disjoint(pt3, p1)
@test GO.disjoint(pt4, p1) == LG.disjoint(pt4, p1)

@test GO.disjoint(pt1, p2) == LG.disjoint(pt1, p2)
@test GO.disjoint(pt2, p2) == LG.disjoint(pt2, p2)
@test GO.disjoint(pt3, p2) == LG.disjoint(pt3, p2)
@test GO.disjoint(pt4, p2) == LG.disjoint(pt4, p2)
@test GO.disjoint(pt5, p2) == LG.disjoint(pt5, p2)
@test GO.disjoint(pt6, p2) == LG.disjoint(pt6, p2)
@test GO.disjoint(pt7, p2) == LG.disjoint(pt7, p2)