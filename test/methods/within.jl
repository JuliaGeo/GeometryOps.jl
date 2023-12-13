# Test points
pt1 = LG.Point([0.0, 0.0])
pt2 = LG.Point([0.0, 0.1])
pt3 = LG.Point([1.0, 0.0])
pt4 = LG.Point([0.5, 1.0])
pt5 = LG.Point([0.2, 0.5])
pt6 = LG.Point([0.3, 0.55])
pt7 = LG.Point([0.6, 0.49])
pt8 = LG.Point([0.25, 0.75])
pt9 = LG.Point([0.5, 0.1])
# Test lines
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
l13 = LG.LineString([[0.5, 0.01], [0.5, 0.09]])
# Test rings
r1 = LG.LinearRing([[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0]])
r2 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1], [0.0, 0.0]])
r3 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [0.0, 0.2], [0.0, 0.0]])
r4 = LG.LinearRing([[0.2, 0.5], [0.3, 0.7], [0.4, 0.5], [0.2, 0.5]])
r5 = LG.LinearRing([[0.6, 0.9], [0.7, 0.8], [0.6, 0.8], [0.6, 0.9]])
r6 = LG.LinearRing([[0.25, 0.55], [0.3, 0.65], [0.35, 0.55], [0.25, 0.55]])
r7 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1], [0.5, 0.3], [0.0, 0.3],[0.0, 0.0]])
# Test polygons
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
# Test multipolygons
m1 = LG.MultiPolygon([p3, p6])
m2 = LG.MultiPolygon([p3, p4])

# # Point and point
# Equal points -> within
@test GO.within(pt1, pt1) == LG.within(pt1, pt1)
# Different points -> not within
@test GO.within(pt1, pt2) == LG.within(pt1, pt2)

# # Point and line
# Line endpoint (1 segment) -> not within
@test GO.within(pt1, l1) == LG.within(pt1, l1)
# Line endpoint (2 segments) -> not within
@test GO.within(pt2, l2) == LG.within(pt2, l2)
# Middle of line (1 segment) -> within
@test GO.within(pt2, l1) == LG.within(pt2, l1)
# Not on line (1 segment) -> not within
@test GO.within(pt3, l1) == LG.within(pt3, l1)
# Middle of line on joint (2 segments) -> within
@test GO.within(pt3, l2) == LG.within(pt3, l2)
# Endpoint on closed line -> within
@test GO.within(pt1, l6) == LG.within(pt1, l6)

# # Point and Ring
# On ring corner -> within
@test GO.within(pt1, r1) == LG.within(pt1, r1)
# Outside of ring -> not within
@test GO.within(pt2, r1) == LG.within(pt2, r1)
# Inside of ring center (not on line) -> not within
@test GO.within(pt3, r1) == LG.within(pt3, r1)
# On ring edge -> within
@test GO.within(pt8, r1) == LG.within(pt8, r1)

# # Point and polygon
# On polygon vertex -> not within
@test GO.within(pt1, p1) == LG.within(pt1, p1)
# Outside of polygon -> not within
@test GO.within(pt2, p1) == LG.within(pt2, p1)
# Inside of polygon -> within
@test GO.within(pt3, p1) == LG.within(pt3, p1)
# On polygon vertex (with holes) -> not within
@test GO.within(pt1, p2) == LG.within(pt1, p2)
# On polygon edge (with holes) -> not within
@test GO.within(pt2, p2) == LG.within(pt2, p2)
# On hole vertex -> not within
@test GO.within(pt5, p2) == LG.within(pt5, p2)
# Within hole -> not within
@test GO.within(pt6, p2) == LG.within(pt6, p2)
# Inside of polygon (with holes) -> within
@test GO.within(pt7, p2) == LG.within(pt7, p2)

# # Geometry and point

# # Line and line
# Equal lines -> within
@test GO.within(l1, l1) == LG.within(l1, l1)
# Lines share 2 endpoints, but don't overlap -> not within
@test GO.within(l1, l2) == LG.within(l1, l2)
# Lines overlap, but neither is within other -> not within
@test GO.within(l1, l3) == LG.within(l1, l3)
# Within line (no shared endpoints) -> within
@test GO.within(l1, l4) == LG.within(l1, l4)
# Line shares just one endpoint -> not within
@test GO.within(l1, l5) == LG.within(l1, l5)

# # Line and ring
# Shares all endpoints -> within
@test GO.within(l6, r1) == LG.within(l6, r1)
# Shares all endpoints, but ring has extra edge -> within
@test GO.within(l2, r2) == LG.within(l2, r2)
# Doesn't share all edges -> not within
@test GO.within(l2, r3) == LG.within(l2, r3)
# Shares all endpoints, but adds one extra segment -> not within
@test GO.within(l12, r1) == LG.within(l12, r1)

# Line and polygon
# Line traces entire outline of polygon edges -> not within
@test GO.within(l6, p1) == LG.within(l6, p1)
# Line is edge of polygon -> not within
@test GO.within(l1, p2) == LG.within(l1, p2)
# Line is on edge + inside of polygon -> within
@test GO.within(l2, p2) == LG.within(l2, p2)
# Line goes outside of polygon -> not within
@test GO.within(l3, p2) == LG.within(l3, p2)
# Line is fully within polygon -> within
@test GO.within(l7, p2) == LG.within(l7, p2)
# Line is fully within hole -> not within
@test GO.within(l8, p2) == LG.within(l8, p2)
# Line is on hole edge -> not within
@test GO.within(l9, p2) == LG.within(l9, p2)
# Line on polygon edge and then enters polygon to end on hole vertex -> within
@test GO.within(l10, p2) == LG.within(l10, p2)
# Line is on polygon edge and then cuts through hole -> not within
@test GO.within(l11, p2) == LG.within(l11, p2)

# # Ring and line
# Shares all endpoints -> within
@test GO.within(r1, l6) == LG.within(r1, l6)
# Shares all endpoints but ring has closing edge -> not within
@test GO.within(r2, l2) == LG.within(r2, l2)
# Doesn't share all edges -> not within
@test GO.within(r3, l2) == LG.within(r3, l2)
# Shares all endpoints, but line has one extra segment -> within
@test GO.within(r1, l12) == LG.within(r1, l12)

# # Ring and ring
# Equal ring -> within
@test GO.within(r1, r1) == LG.within(r1, r1)
# Not equal ring -> not within
@test GO.within(r1, r2) == LG.within(r1, r2)
# Not equal ring -> not within
@test GO.within(r1, r3) == LG.within(r1, r3)
# Rings share all edges, but second ring has extra edges -> within
@test GO.within(r2, r7) == LG.within(r2, r7)

# # Ring and polygon
# Ring is equal to polygon's external ring, no holes -> not within
@test GO.within(r1, p1) == LG.within(r1, p1)
# Ring goes outside of polygon's external ring -> not within
@test GO.within(r1, p2) == LG.within(r1, p1)
# Ring is within polygon, but also on edges -> within
@test GO.within(r2, p2) == LG.within(r2, p2)
# Ring is within polygon, but also on edges -> within
@test GO.within(r3, p2) == LG.within(r3, p2)
# Ring is one of polygon's holes -> not within
@test GO.within(r4, p2) == LG.within(r4, p2)
# Ring is fully within polygon that has holes -> within
@test GO.within(r5, p2) == LG.within(r5, p2)
# Ring is fully within polygon's hole -> not within
@test GO.within(r6, p2) == LG.within(r6, p2)

# # Polygon in polygon
# Same polygon -> within
@test GO.within(p1, p1) == LG.within(p1, p1)
@test GO.within(p2, p2) == LG.within(p2, p2)
# Polygon not in polygon -> not within
@test GO.within(p1, p2) == LG.within(p1, p2)
@test GO.within(p2, p1) == LG.within(p2, p1)
# Polygon is within polygon, but also on edges -> within
@test GO.within(p3, p2) == LG.within(p3, p2)
# Polygon within polygon with holes -> within
@test GO.within(p4, p2) == LG.within(p4, p2)
# Polygon within polygon hole --> not within
@test GO.within(p5, p2) == LG.within(p5, p2)
# Polygon overlapping with other polygon's hole -> not within
@test GO.within(p6, p2) == LG.within(p6, p2)
# Polygon with hole nested with other polygon's hole --> within
@test GO.within(p8, p7) == LG.within(p8, p7)
# Nested holes but not within -> not within
@test GO.within(p9, p7) == LG.within(p9, p7)
# Nested with same hole -> within
@test GO.within(p10, p7) == LG.within(p10, p7)
# within external ring but intersects with hole -> not within
@test GO.within(p11, p7) == LG.within(p11, p7)
# polygon extactly overlaps with other polygon's hole -> not within
@test GO.within(p12, p7) == LG.within(p12, p7)

# # Multipolygon tests
# Point in multipolygon
@test GO.within(pt5, m1) == LG.within(pt5, m1)
@test GO.within(pt9, m1) == LG.within(pt9, m1)
# Point outside of multipolygon
@test GO.within(pt4, m1) == LG.within(pt4, m1)
# Line in multipolygon
@test GO.within(l13, m1) == LG.within(l13, m1)
@test GO.within(l9, m1) == LG.within(l9, m1)
# Line outside of multipolygon
@test GO.within(l1, m1) == LG.within(l1, m1)
# Ring in multipolygon
@test GO.within(r1, m2) == LG.within(r1, m2)
# Ring outside of multipolygon
@test GO.within(r1, m1) == LG.within(r1, m1)
# Polygon in multipolygon
@test GO.within(p3, m1) == LG.within(p3, m1)
@test GO.within(p6, m1) == LG.within(p6, m1)
# Polygon outside of multipolygon
@test GO.within(p1, m1) == LG.within(p1, m1)
# Multipolygon in multipolygon
@test GO.within(m1, m1) == LG.within(m1, m1)
# Multipolygon outside of multipolygon
@test GO.within(m2, m1) == LG.within(m2, m1)
