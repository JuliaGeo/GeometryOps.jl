p1 = LG.Point([0.0, 0.0])
p2 = LG.Point([0.0, 0.1])
p3 = LG.Point([1.0, 0.0])

l1 = LG.LineString([[0.0, 0.0], [0.0, 1.0]])
l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1]])

# Point and point
@test GO.touches(p1, p1) == LG.touches(p1, p1)
@test GO.touches(p1, p2) == LG.touches(p1, p2)

# Point and line
@test GO.touches(p1, l1) == LG.touches(p1, l1)
@test GO.touches(p2, l1) == LG.touches(p2, l1)
@test GO.touches(p3, l1) == LG.touches(p3, l1)
@test GO.touches(p1, l2) == LG.touches(p1, l2)
@test GO.touches(p2, l2) == LG.touches(p2, l2)
@test GO.touches(p3, l2) == LG.touches(p3, l2)

# Line and line
LG.touches(l1, l1)
LG.touches(l1, l2)
LG.touches(l1, l3)
LG.touches(l1, l4)
LG.touches(l1, l5)