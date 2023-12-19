pt1 = LG.Point([0.0, 0.0])
pt2 = LG.Point([0.0, 0.1])
pt3 = LG.Point([1.0, 0.0])

l1 = LG.LineString([[0.0, 0.0], [0.0, 1.0]])
l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1]])
l3 = LG.LineString([[0.0, -1.0], [0.0, 0.5]])
l4 = LG.LineString([[0.0, -1.0], [0.0, 1.5]])
l5 =  LG.LineString([[0.0, -1.0], [0.0, 0.0]])

# Point and point
@test GO.touches(pt1, pt1) == LG.touches(pt1, pt1)
@test GO.touches(pt1, pt2) == LG.touches(pt1, pt2)

# Point and line
@test GO.touches(pt1, l1) == LG.touches(pt1, l1)
@test GO.touches(pt2, l1) == LG.touches(pt2, l1)
@test GO.touches(pt3, l1) == LG.touches(pt3, l1)
@test GO.touches(pt1, l2) == LG.touches(pt1, l2)
@test GO.touches(pt2, l2) == LG.touches(pt2, l2)
@test GO.touches(pt3, l2) == LG.touches(pt3, l2)

# Line and line
@test GO.touches(l1, l1) == LG.touches(l1, l1)
@test GO.touches(l1, l2) == LG.touches(l1, l2)
@test GO.touches(l1, l3) == LG.touches(l1, l3)
@test GO.touches(l1, l4) == LG.touches(l1, l4)
@test GO.touches(l1, l5) == LG.touches(l1, l5)