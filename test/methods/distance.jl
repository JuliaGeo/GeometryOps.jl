pt1 = LG.Point([0.0, 0.0])
pt2 = LG.Point([0.0, 1.0])
pt3 = LG.Point([2.5, 2.5])
pt4 = LG.Point([3.0, 3.0])
pt5 = LG.Point([5.1, 5.0])
pt6 = LG.Point([3.0, 1.0])
pt7 = LG.Point([0.1, 4.9])
pt8 = LG.Point([2.0, 1.1])
pt9 = LG.Point([3.5, 3.1])
pt10 = LG.Point([10.0, 10.0])
pt11 = LG.Point([2.5, 7.0])

mpt1 = LG.MultiPoint([pt1, pt2, pt3])

l1 = LG.LineString([[0.0, 0.0], [0.0, 5.0], [5.0, 5.0]])

r1 = LG.LinearRing([[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [0.0, 0.0]])
r2 = LG.LinearRing([[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]])
r3 = LG.LinearRing([[1.0, 1.0], [3.0, 2.0], [4.0, 1.0], [1.0, 1.0]])
r4 = LG.LinearRing([[4.0, 3.0], [3.0, 3.0], [4.0, 4.0], [4.0, 3.0]])
r5 = LG.LinearRing([[0.0, 6.0], [2.5, 8.0], [5.0, 6.0], [0.0, 6.0]])

p1 = LG.Polygon(r2, [r3, r4])
p2 = LG.Polygon(r5)

mp1 = LG.MultiPolygon([p1, p2])

c1 = LG.GeometryCollection([pt1, r1, p1])

# Point and Point

# Distance from point to same point
@test GO.distance(pt1, pt1) == LG.distance(pt1, pt1)
# Distance from point to different point
@test GO.distance(pt1, pt2) ≈ GO.distance(pt2, pt1) ≈ LG.distance(pt1, pt2)
# Return types
@test GO.distance(pt1, pt1) isa Float64
@test GO.distance(pt1, pt1, Float32) isa Float32

# Point and Line

#Point on line vertex
@test GO.distance(pt1, l1) ==  GO.distance(l1, pt1) == LG.distance(pt1, l1)
# Point on line edge
@test GO.distance(pt2, l1) ==  GO.distance(l1, pt2) == LG.distance(pt2, l1)
# Point equidistant from both segments
@test GO.distance(pt3, l1) ≈  GO.distance(l1, pt3) ≈ LG.distance(pt3, l1)
# Point closer to one segment than another
@test GO.distance(pt4, l1) ≈  GO.distance(l1, pt4) ≈ LG.distance(pt4, l1)
# Return types
@test GO.distance(pt1, l1) isa Float64
@test GO.distance(pt1, l1, Float32) isa Float32

# Point and Ring

# Point on linear ring
@test GO.distance(pt1, r1) == LG.distance(pt1, r1)
@test GO.distance(pt3, r1) == LG.distance(pt3, r1)
# Point outside of linear ring
@test GO.distance(pt5, r1) ≈ LG.distance(pt5, r1)
# Point inside of hole created by linear ring
@test GO.distance(pt3, r1) ≈ LG.distance(pt3, r1)
@test GO.distance(pt4, r1) ≈ LG.distance(pt4, r1)

# Point and Polygon
# Point on polygon exterior edge
@test GO.distance(pt1, p1) == LG.distance(pt1, p1)
@test GO.signed_distance(pt1, p1) == 0
@test GO.distance(pt2, p1) == LG.distance(pt2, p1)
# Point on polygon hole edge
@test GO.distance(pt4, p1) == LG.distance(pt4, p1)
@test GO.signed_distance(pt4, p1) == 0
@test GO.distance(pt6, p1) == LG.distance(pt6, p1)
# Point inside of polygon
@test GO.distance(pt3, p1) == LG.distance(pt3, p1)
@test GO.signed_distance(pt3, p1) ≈
    -(min(LG.distance(pt3, r2), LG.distance(pt3, r3), LG.distance(pt3, r4)))
@test GO.distance(pt7, p1) == LG.distance(pt7, p1)
@test GO.signed_distance(pt7, p1) ≈
    -(min(LG.distance(pt7, r2), LG.distance(pt7, r3), LG.distance(pt7, r4)))
# Point outside of polyon exterior
@test GO.distance(pt5, p1) ≈ LG.distance(pt5, p1)
@test GO.signed_distance(pt5, p1) ≈ LG.distance(pt5, p1)
# Point inside of polygon hole
@test GO.distance(pt8, p1) ≈ LG.distance(pt8, p1)
@test GO.signed_distance(pt8, p1) ≈ LG.distance(pt8, p1)
@test GO.distance(pt9, p1) ≈ LG.distance(pt9, p1)
# Return types
@test GO.distance(pt1, p1) isa Float64
@test GO.distance(pt1, p1, Float32) isa Float32

# Point and MultiPoint
@test GO.distance(pt4, mpt1) == LG.distance(pt4, mpt1)
@test GO.distance(pt4, mpt1) isa Float64
@test GO.distance(pt4, mpt1, Float32) isa Float32

# Point and MultiPolygon
# Point outside of either polygon
@test GO.distance(pt5, mp1) ≈ LG.distance(pt5, mp1)
@test GO.distance(pt10, mp1) ≈ LG.distance(pt10, mp1)
# Point within one polygon
@test GO.distance(pt3, mp1) == LG.distance(pt3, mp1)
@test GO.signed_distance(pt3, mp1) ≈
    -(min(LG.distance(pt3, r2), LG.distance(pt3, r3), LG.distance(pt3, r4), LG.distance(pt3, r5)))
@test GO.distance(pt11, mp1) == LG.distance(pt11, mp1)
@test GO.signed_distance(pt11, mp1) ≈
    -(min(LG.distance(pt11, r2), LG.distance(pt11, r3), LG.distance(pt11, r4), LG.distance(pt11, r5)))

# Point and Geometry Collection
@test GO.distance(pt1, c1) == LG.distance(pt1, c1)

