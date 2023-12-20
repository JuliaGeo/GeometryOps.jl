pt = LG.Point([0.0, 0.0])
l1 = LG.LineString([[0.0, 0.0], [0.5, 0.5], [1.0, 0.5]])
r1 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 0.0]])

@test GO.area(pt) == GO.signed_area(pt) == LG.area(pt)
@test GO.area(l1) == GO.signed_area(l1) == LG.area(l1)
@test GO.area(r1) == GO.signed_area(r1) == LG.area(r1)
p1 = LG.Polygon([[[10.0, 0.0], [30.0, 0.0], [30.0, 20.0], [10.0, 20.0], [10.0, 0.0]]])
p2 = LG.Polygon([
    [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
    [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]]
])
@test GO.area(p1) == GO.signed_area(p1) == LG.area(p1)
@test GO.area(p2) == LG.area(p2)
@test GO.signed_area(p2) == -LG.area(p2)
