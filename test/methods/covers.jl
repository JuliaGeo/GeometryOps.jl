p1 = LG.Point([0.0, 0.0])
p2 = LG.Point([0.0, 0.1])
p3 = LG.Point([1.0, 0.0])

l1 = LG.LineString([[0.0, 1.0]])
l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1]])

# Point and point
@test GO.covers(p1, p1) == LG.covers(p1, p1)
@test GO.covers(p1, p2) == LG.covers(p1, p2)

# Point and line
@test GO.covers(l1, p1) == LG.covers(l1, p1)
@test GO.covers(l1, p2) == LG.covers(l1, p2)
@test GO.covers(l1, p3) == LG.covers(l1, p3)
@test GO.covers(l2, p1) == LG.covers(l2, p1)
@test GO.covers(l2, p2) == LG.covers(l2, p2)
@test GO.covers(l2, p3) == LG.covers(l2, p3)

# Line and line
@test GO.covers(l1, l1) == LG.covers(l1, l1)
@test GO.covers(l2, l1) == LG.covers(l2, l1)
@test GO.covers(l3, l1) == LG.covers(l3, l1)
@test GO.covers(l4, l1) == LG.covers(l4, l1)
@test GO.covers(l5, l1) == LG.covers(l5, l1)