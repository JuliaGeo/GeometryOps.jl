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
r4 = LG.LinearRing([[0.2, 0.5], [0.3, 0.7], [0.4, 0.5], [0.2, 0.5]])
r5 = LG.LinearRing([[0.6, 0.9], [0.7, 0.8], [0.6, 0.8], [0.6, 0.9]])
r6 = LG.LinearRing([[0.25, 0.55], [0.3, 0.65], [0.35, 0.55], [0.25, 0.55]])

p1 = LG.Polygon([[[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0]]])
p2 = LG.Polygon([
    [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]],
    [[0.2, 0.5], [0.3, 0.7], [0.4, 0.5], [0.2, 0.5]],
    [[0.5, 0.5], [0.8, 0.4], [0.7, 0.2], [0.5, 0.5]]
])
p3 = LG.Polygon([[[0.0, 0.0], [1.0, 0.0], [0.0, 0.2], [0.0, 0.0]]])
p4 = LG.Polygon([[[0.6, 0.9], [0.7, 0.8], [0.6, 0.8], [0.6, 0.9]]])
p5 = LG.Polygon([[[0.25, 0.55], [0.3, 0.65], [0.35, 0.55], [0.25, 0.55]]])
p6 = LG.Polygon([[[0.1, 0.4], [0.1, 0.8], [0.3, 0.8], [0.3, 0.4], [0.1, 0.4]]])
p7 = LG.Polygon([
    [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]],
    [[0.4, 0.4], [0.4, 0.6], [0.6, 0.6], [0.6, 0.4], [0.4, 0.4]],
])
p8 = LG.Polygon([
    [[0.1, 0.1], [0.1, 0.9], [0.9, 0.9], [0.9, 0.1], [0.1, 0.1]],
    [[0.3, 0.3], [0.3, 0.7], [0.7, 0.7], [0.7, 0.3], [0.3, 0.3]]
])
p9 = LG.Polygon([
    [[0.1, 0.1], [0.1, 0.9], [0.9, 0.9], [0.9, 0.1], [0.1, 0.1]],
    [[0.45, 0.45], [0.45, 0.55], [0.55, 0.55], [0.55, 0.45], [0.45, 0.45]]
])
p10 = LG.Polygon([
    [[0.1, 0.1], [0.1, 0.9], [0.9, 0.9], [0.9, 0.1], [0.1, 0.1]],
    [[0.4, 0.4], [0.4, 0.6], [0.6, 0.6], [0.6, 0.4], [0.4, 0.4]]
])
p11 = LG.Polygon([
    [[0.45, 0.5], [0.45, 0.75], [0.55, 0.75], [0.55, 0.5], [0.45, 0.5]]
])
p12 = LG.Polygon([
    [[0.4, 0.4], [0.4, 0.6], [0.6, 0.6], [0.6, 0.4], [0.4, 0.4]]
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
# Ring is equal to polygon's external ring, no holes
@test GO.within(r1, p1) == LG.within(r1, p1)
# Ring goes outside of polygon's external ring
@test GO.within(r1, p2) == LG.within(r1, p1)
# Ring is within polygon, but also on edges
@test GO.within(r2, p2) == LG.within(r2, p2)
# Ring is within polygon, but also on edges
@test GO.within(r3, p2) == LG.within(r3, p2)
# Ring is one of polygon's holes
@test GO.within(r4, p2) == LG.within(r4, p2)
# Ring is fully within polygon that has holes
@test GO.within(r5, p2) == LG.within(r5, p2)
# Ring is fully within polygon's hole
@test GO.within(r6, p2) == LG.within(r6, p2)

# Polygon in polygon

# Same polygon
@test GO.within(p1, p1) == LG.within(p1, p1)
@test GO.within(p2, p2) == LG.within(p2, p2)
# Polygon not in polygon
@test GO.within(p1, p2) == LG.within(p1, p2)
@test GO.within(p2, p1) == LG.within(p2, p1)
# Ring is within polygon, but also on edges
@test GO.within(p3, p2) == LG.within(p3, p2)
# Polygon within polygon with holes
@test GO.within(p4, p2) == LG.within(p4, p2)
# Polygon within polygon hole --> not within
@test GO.within(p5, p2) == LG.within(p5, p2)
# Polygon overlapping with other polygon's hole
@test GO.within(p6, p2) == LG.within(p6, p2)
# Polygon with hole nested with other polygon's hole --> within
@test GO.within(p8, p7) == LG.within(p8, p7)
# Nested holes but not within
@test GO.within(p9, p7) == LG.within(p9, p7)
# Nested with same hole
@test GO.within(p10, p7) == LG.within(p10, p7)
# within external ring but intersects with hole
@test GO.within(p11, p7) == LG.within(p11, p7)

@test GO.within(p12, p7) == LG.within(p12, p7)


