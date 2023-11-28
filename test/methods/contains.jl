p1 = LG.Point([0.0, 0.0])
p2 = LG.Point([0.0, 0.1])
p3 = LG.Point([1.0, 0.0])

l1 = LG.LineString([[0.0, 1.0]])
l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1]])

# Point and point
@test GO.contains(p1, p1) == LG.contains(p1, p1)
@test GO.contains(p1, p2) == LG.contains(p1, p2)

# Point and line
@test GO.contains(l1, p1) == LG.contains(l1, p1)
@test GO.contains(l1, p2) == LG.contains(l1, p2)
@test GO.contains(l1, p3) == LG.contains(l1, p3)
@test GO.contains(l2, p1) == LG.contains(l2, p1)
@test GO.contains(l2, p2) == LG.contains(l2, p2)
@test GO.contains(l2, p3) == LG.contains(l2, p3)

# Line and line
@test GO.contains(l1, l1) == LG.contains(l1, l1)
@test GO.contains(l1, l2) == LG.contains(l1, l2)
@test GO.contains(l1, l3) == LG.contains(l1, l3)
@test GO.contains(l1, l4) == LG.contains(l1, l4)
@test GO.contains(l1, l5) == LG.contains(l1, l5)

@test GO.contais(l1, l1) == LG.contains(l1, l1)
@test GO.contais(l2, l1) == LG.contains(l2, l1)
@test GO.contais(l3, l1) == LG.contains(l3, l1)
@test GO.contais(l4, l1) == LG.contains(l4, l1)
@test GO.contais(l5, l1) == LG.contains(l5, l1)