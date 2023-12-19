# Test points
pt1 = LG.Point([0.0, 0.0])
pt2 = LG.Point([0.0, 0.1])
pt3 = LG.Point([1.0, 0.0])
pt4 = LG.Point([0.5, 1.0])
pt5 = LG.Point([0.2, 0.5])
pt6 = LG.Point([0.3, 0.55])
pt7 = LG.Point([0.6, 0.49])
pt8 = LG.Point([0.25, 0.75])
pt9 = LG.Point([-1.0, 0.0])
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
r5 = LG.LinearRing([[5.0, 5.0], [6.0, 6.0], [7.0, 5.0], [5.0, 5.0]])
r6 = LG.LinearRing([[0.25, 0.55], [0.3, 0.65], [0.35, 0.55], [0.25, 0.55]])
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
p7 = LG.Polygon([[[-2.0, 0.0], [-1.0, 0.0], [-1.5, 1.5], [-2.0, 0.0]]])
p8 = LG.Polygon([
    [[0.4, 0.4], [0.4, 0.6], [0.6, 0.6], [0.6, 0.4], [0.4, 0.4]]
])
# Test multipolygons
m1 = LG.MultiPolygon([p3, p6])
m2 = LG.MultiPolygon([p3, p4])
m3 = LG.MultiPolygon([p2, p7])
m4 = LG.MultiPolygon([p7])

# # Point and point
# Equal points -> not disjoint
@test GO.disjoint(pt1, pt1) == LG.disjoint(pt1, pt1)
# Non-equal points -> disjoint
@test GO.disjoint(pt1, pt2) == LG.disjoint(pt1, pt2)

# # Point and line
# Line endpoint (1 segment) -> not disjoint
@test GO.disjoint(pt1, l1) == LG.disjoint(pt1, l1)
# Middle of line (1 segment) -> not disjoint
@test GO.disjoint(pt2, l1) == LG.disjoint(pt2, l1)
# Not on line (1 segment) -> disjoint
@test GO.disjoint(pt3, l1) == LG.disjoint(pt3, l1)
# Line endpoint (2 segments) -> not disjoing
@test GO.disjoint(pt2, l2) == LG.disjoint(pt2, l2)
# Middle of line on joint (2 segments) -> not disjoint
@test GO.disjoint(pt3, l2) == LG.disjoint(pt3, l2)
# Endpoint on closed line -> not disjoint
@test GO.disjoint(pt1, l6) == LG.disjoint(pt1, l6)

# # Point and ring
# On ring corner -> not disjoint
@test GO.disjoint(pt1, r1) == LG.disjoint(pt1, r1)
# Outside of ring -> disjoint
@test GO.disjoint(pt2, r1) == LG.disjoint(pt2, r1)
# Inside of ring center (not on line) -> disjoint
@test GO.disjoint(pt3, r1) == LG.disjoint(pt3, r1)
# On ring edge -> not disjoint
@test GO.disjoint(pt8, r1) == LG.disjoint(pt8, r1)

# # Point and polygon
# Point on polygon vertex -> not disjoint
@test GO.disjoint(pt1, p2) == LG.disjoint(pt1, p2)
# Point on polygon edge -> not disjoint
@test GO.disjoint(pt2, p2) == LG.disjoint(pt2, p2)
# Point on edge of hold --> not disjoint
@test GO.disjoint(pt5, p2) == LG.disjoint(pt5, p2)
# Point in hole -> disjoint
@test GO.disjoint(pt6, p2) == LG.disjoint(pt6, p2)
# Point inside of polygon -> not disjoint
@test GO.disjoint(pt7, p2) == LG.disjoint(pt7, p2)
# Point outside of polygon -> disjoint
@test GO.disjoint(pt9, p2) == LG.disjoint(pt9, p2)

# # Geometry and point (switched direction)
@test GO.disjoint(pt1, l1) == GO.disjoint(l1, pt1)
@test GO.disjoint(pt1, r1) == GO.disjoint(r1, pt1)
@test GO.disjoint(pt1, p2) == GO.disjoint(p2, pt1)

# # Line and line
# Equal lines -> not disjoint
@test GO.disjoint(l1, l1) == LG.disjoint(l1, l1)
# Lines share 2 endpoints, but don't overlap -> not disjoint
@test GO.disjoint(l1, l2) == LG.disjoint(l1, l2)
# Lines overlap, but neither is within other -> not disjoint
@test GO.disjoint(l1, l3) == LG.disjoint(l1, l3)
# Within line (no shared endpoints) -> not disjoint
@test GO.disjoint(l1, l4) == LG.disjoint(l1, l4)
# Line shares just 1 endpoint -> not disjoint
@test GO.disjoint(l1, l5) == LG.disjoint(l1, l5)
# Lines don't touch at all -> disjoint
@test GO.disjoint(l7, l1) == LG.disjoint(l7, l1)

# # Line and ring
# Shares all endpoints -> not disjoint
@test GO.disjoint(l6, r1) == LG.disjoint(l6, r1)
# Shares only some edges -> not disjoint
@test GO.disjoint(l2, r3) == LG.disjoint(l2, r3)
# line inside of ring -> disjoint
@test GO.disjoint(l7, r1) == LG.disjoint(l7, r1)
# line outside of ring -> disjoint
@test GO.disjoint(l7, r2) == LG.disjoint(l7, r2)

# # # Line and polygon
# # Line traces entire outline of polygon edges -> not disjoint
# @test GO.disjoint(l6, p1) == LG.disjoint(l6, p1)
# # Line is on edge + inside of polygon -> not disjoint
# @test GO.disjoint(l2, p2) == LG.disjoint(l2, p2)
# # Line goes outside of polygon -> not disjoint
# @test GO.disjoint(l3, p2) == LG.disjoint(l3, p2)
# # Line is fully within hole -> disjoint
# @test GO.disjoint(l8, p2) == LG.disjoint(l8, p2)
# # Line is on polygon edge and then cuts through hole -> not disjoint
# @test GO.disjoint(l11, p2) == LG.disjoint(l11, p2)

# # Geometry and line (switched direction)
@test GO.disjoint(l7, r1) == GO.disjoint(r1, l7)
# @test GO.disjoint(l8, p2) == GO.disjoint(p2, l8)

# # Ring and line
# Shares all endpoints -> not disjoint
@test GO.disjoint(r1, l6) == LG.disjoint(r1, l6)
# Shares some edges -> not disjoint
@test GO.disjoint(r3, l2) == LG.disjoint(r3, l2)
# Doesn't share any edges -> disjoint
@test GO.disjoint(r4, l2) == LG.disjoint(r4, l2)

# # Ring and ring
# Equal ring -> not disjoint
@test GO.disjoint(r1, r1) == LG.disjoint(r1, r1)
# Not equal ring but share a vertex -> not disjoint
@test GO.disjoint(r1, r2) == LG.disjoint(r1, r2)
# Rings not touching -> not disjoint
@test GO.disjoint(r3, r4) == LG.disjoint(r3, r4)
# Ring inside of ring -> disjoint
@test GO.disjoint(r4, r2) == LG.disjoint(r4, r2)
# Ring outside of other ring -> disjoint
@test GO.disjoint(r2, r4) == LG.disjoint(r2, r4)

# Ring and polygon
# Ring goes outside of polygon's external ring -> not disjoint
@test GO.within(r1, p2) == LG.within(r1, p1)
# Ring is one of polygon's holes -> not disjoint
@test GO.disjoint(r4, p2) == LG.disjoint(r4, p2)
# Ring is fully within polygon -> not disjoint
@test GO.disjoint(r5, p2) == LG.disjoint(r5, p2)
# Ring is fully within polygon's hole -> disjoint
@test GO.disjoint(r6, p2) == LG.disjoint(r6, p2)
# Ring is fully outside of the polygon -> disjoint
@test GO.disjoint(r5, p2) == LG.disjoint(r5, p2)

# # Geometry and ring (switched direction)
@test GO.disjoint(r4, r2) == GO.disjoint(r2, r4)
@test GO.disjoint(p2, r6) == GO.disjoint(r6, p2)

# # Polygon and polygon
# Overlapping polygons -> not disjoint
@test GO.disjoint(p1, p2) == LG.disjoint(p1, p2)
# Polygon is within polygon, but also on edges -> not disjoint
@test GO.disjoint(p3, p2) == LG.disjoint(p3, p2)
# Polygon within polygon hole --> disjoint
@test GO.disjoint(p5, p2) == LG.disjoint(p5, p2)
# polygon extactly overlaps with other polygon's hole -> not disjoint
@test GO.disjoint(p8, p7) == LG.disjoint(p8, p7)

# # Multipolygon tests
# Point in multipolygon -> not disjoint
@test GO.disjoint(pt5, m1) == LG.disjoint(pt5, m1)
# Point outside of multipolygon -> disjoint
@test GO.disjoint(pt4, m1) == LG.disjoint(pt4, m1)
# Line in multipolygon -> not disjoint
@test GO.disjoint(l13, m1) == LG.disjoint(l13, m1)
# Line outside of multipolygon -> disjoint
@test GO.disjoint(l8, m3) == LG.disjoint(l8, m3)
# Ring in multipolygon -> not disjoint
@test GO.disjoint(r1, m2) == LG.disjoint(r1, m2)
# Ring outside of multipolygon
@test GO.disjoint(r6, m3) == LG.disjoint(r6, m3)
# Polygon in multipolygon -> not disjoint
@test GO.disjoint(p3, m1) == LG.disjoint(p3, m1)
# Polygon outside of multipolygon -> disjoint
@test GO.disjoint(p5, m3) == LG.disjoint(p5, m3)
# Multipolygon in multipolygon -> not disjoint
@test GO.disjoint(m1, m1) == LG.disjoint(m1, m1)
# Multipolygon outside of multipolygon -> disjoint
@test GO.disjoint(m1, m4) == LG.disjoint(m1, m4)