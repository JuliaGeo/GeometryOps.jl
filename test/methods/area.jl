pt = LG.Point([0.0, 0.0])
l1 = LG.LineString([[0.0, 0.0], [0.5, 0.5], [1.0, 0.5]])
r1 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 0.0]])
p1 = LG.Polygon([
    [[10.0, 0.0], [30.0, 0.0], [30.0, 20.0], [10.0, 20.0], [10.0, 0.0]],
])
p2 = LG.Polygon([
    [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
    [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]],
])
p3 = LG.Polygon([
    [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
    [[15.0, 1.0], [25.0, 1.0], [25.0, 11.0], [15.0, 11.0], [15.0, 1.0]],
])
p4 = LG.Polygon([
    [
        [0.0, 5.0], [2.0, 2.0], [5.0, 2.0], [2.0, -2.0], [5.0, -5.0],
        [0.0, -2.0], [-5.0, -5.0], [-2.0, -2.0], [-5.0, 2.0], [-2.0, 2.0],
        [0.0, 5.0],
    ],
])
mp1 = LG.MultiPolygon([p2, p4])


# Points, lines, and rings have zero area
@test GO.area(pt) == GO.signed_area(pt) == LG.area(pt) == 0
@test GO.area(l1) == GO.signed_area(l1) == LG.area(l1) == 0
@test GO.area(r1) == GO.signed_area(r1) == LG.area(r1) == 0

# Polygons have non-zero area
# CCW polygons have positive signed area
@test GO.area(p1) == GO.signed_area(p1) == LG.area(p1)
@test GO.signed_area(p1) > 0
# CW polygons have negative signed area
a2 = LG.area(p2)
@test GO.area(p2) == a2
@test GO.signed_area(p2) == -a2
# Winding order of holes doesn't affect sign of signed area
@test GO.signed_area(p3) == -a2
# Concave polygon correctly calculates area
a4 = LG.area(p4)
@test GO.area(p4) == a4
@test GO.signed_area(p4) == -a4
# Multipolygon calculations work
@test GO.area(mp1) == a2 + a4
