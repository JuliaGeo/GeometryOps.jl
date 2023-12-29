# # Point and Geometry

# Same point -> contains
@test GO.contains(pt1, pt1) == LG.contains(pt1, pt1) == true
# Different point -> doesn't contain
@test GO.contains(pt1, pt2) == LG.contains(pt1, pt2) == false

# Point on line endpoint -> does not contain
@test GO.contains(l1, pt1) == LG.contains(l1, pt1) == false
# Point outside line -> does not contain
@test GO.contains(l1, pt2) == LG.contains(l1, pt2) == false
# Point on line segment -> contains
@test GO.contains(l1, pt3) == LG.contains(l1, pt3) == true
# Point on line vertex between segments -> contain
@test GO.contains(l1, pt4) == LG.contains(l1, pt4) == true
# Point cannot contain a line -> doesn't contain
@test GO.contains(pt3, l1) == LG.contains(pt3, l1) == false

# Point on ring endpoint -> contains
@test GO.contains(r1, pt1) == LG.contains(r1, pt1) == true
# Point outside ring -> does not contain
@test GO.contains(r1, pt2) == LG.contains(r1, pt2) == false
# Point on ring segment -> contains
@test GO.contains(r1, pt3) == LG.contains(r1, pt3) == true
# Point on ring vertex between segments -> contain
@test GO.contains(r1, pt4) == LG.contains(r1, pt4) == true
# Point cannot contain a ring -> doesn't contain
@test GO.contains(pt3, r1) == LG.contains(pt3, r1) == false

# Point on vertex of polygon --> doesn't contain
@test GO.contains(p1, pt1) == LG.contains(p1, pt1) == false
# Point outside of polygon's external ring -> doesn't contain
@test GO.contains(p1, pt2) == LG.contains(p1, pt2) == false
# Point on polygon's edge -> doesn't contain
@test GO.contains(p1, pt4) == LG.contains(p1, pt4) == false
# Point inside of polygon -> contains
@test GO.contains(p1, pt5) == LG.contains(p1, pt5) == true
# Point on hole edge -> doesn't contain
@test GO.contains(p1, pt6) == LG.contains(p1, pt6) == false
# Point inside of polygon hole -> doesn't contain
@test GO.contains(p1, pt7) == LG.contains(p1, pt7) == false
# Point cannot contain a polygon -> doesn't contain
@test GO.contains(pt5, p1) == LG.contains(pt5, p1) == false

# # Line and Geometry

# Same line -> contains
@test GO.contains(l1, l1) == LG.contains(l1, l1) == true
# Line overlaps line edge and endpoint -> contains
@test GO.contains(l1, l2) == LG.contains(l1, l2) == true
# Line overlaps with one edge and is outside of other edge -> doesn't contain
@test GO.contains(l1, l3) == LG.contains(l1, l3) == false
# Line segments both within other line segments -> contain
@test GO.contains(l1, l4) == LG.contains(l1, l4) == true
# Line segments connect at endpoint -> doesn't contain
@test GO.contains(l1, l5) == LG.contains(l1, l5) == false
# Line segments don't touch -> doesn't contain
@test GO.contains(l1, l6) == LG.contains(l1, l6) == false
# Line segments cross -> doesn't contain
@test GO.contains(l1, l7) == LG.contains(l1, l7) == false
# Line segments cross and go over and out -> doesn't contain
@test GO.contains(l1, l8) == LG.contains(l1, l8) == false
# Line segments cross and overlap on endpoint -> doesn't contain
@test GO.contains(l1, l9) == LG.contains(l1, l9) == false

# Line is within linear ring -> doesn't contain
@test GO.contains(l1, r1) == LG.contains(l1, r1) == false
# Line covers one edge of linera ring and has segment outside -> doesn't contain
@test GO.contains(l3, r1) == LG.contains(l3, r1) == false
# Line and linear ring are only connected at vertex -> doesn't contain
@test GO.contains(l5, r1) == LG.contains(l5, r1) == false
# Line and linear ring are disjoint -> doesn't contain
@test GO.contains(l6, r1) == LG.contains(l6, r1) == false
# Line crosses through two ring edges -> doesn't contain
@test GO.contains(l7, r1) == LG.contains(l7, r1) == false
# Line crosses through two ring edges and touches third edge -> doesn't contain
@test GO.contains(l8, r1) == LG.contains(l8, r1) == false
# Line is equal to linear ring -> contain
@test GO.contains(l10, r1) == LG.contains(l10, r1) == true
# Line covers linear ring and then has extra segment -> contain
@test GO.contains(l11, r1) == LG.contains(l11, r1) == true

# Line on polygon edge -> doesn't contain
@test GO.contains(p1, l1) == LG.contains(p1, l1) == false
# Line on polygon edge and extending beyond polygon edge -> doesn't contain
@test GO.contains(p1, l3) == LG.contains(p1, l3) == false
# Line outside polygon connected by an vertex -> doesn't contain
@test GO.contains(p1, l5) == LG.contains(p1, l5) == false
# Line through polygon cutting to outside -> doesn't contain
@test GO.contains(p1, l7) == LG.contains(p1, l7) == false
# Line inside of polygon -> contains
@test GO.contains(p1, l12) == LG.contains(p1, l12) == true
# Line outside of polygon -> doesn't contain
@test GO.contains(p1, l13) == LG.contains(p1, l13) == false
# Line in polygon hole -> doesn't contain
@test GO.contains(p1, l14) == LG.contains(p1, l14) == false
# Line outside crown polygon but touching edges -> doesn't contains
@test GO.contains(p8, l15) == LG.contains(p8, l15) == false
# Line within crown polygon but touching edges -> contains
@test GO.contains(p9, l15) == LG.contains(p9, l15) == true

# # Ring and Geometry

# Line is within linear ring -> contains
@test GO.contains(r1, l1) == LG.contains(r1, l1) == true
# Line covers one edge of linera ring and has segment outside -> doesn't contain
@test GO.contains(r1, l3) == LG.contains(r1, l3) == false
# Line and linear ring are only connected at vertex -> doesn't contain
@test GO.contains(r1, l5) == LG.contains(r1, l5) == false
# Line and linear ring are disjoint -> doesn't contain
@test GO.contains(r1, l6) == LG.contains(r1, l6) == false
# Line crosses through two ring edges -> doesn't contain
@test GO.contains(r1, l7) == LG.contains(r1, l7) == false
# Line crosses through two ring edges and touches third edge -> doesn't contain
@test GO.contains(r1, l8) == LG.contains(r1, l8) == false
# Line is equal to linear ring -> contain
@test GO.contains(r1, l10) == LG.contains(r1, l10) == true
# Line covers linear ring and then has extra segment -> doesn't contain
@test GO.contains(r1, l11) == LG.contains(r1, l11) == false

# Same ring -> contains
@test GO.contains(r1, r1) == LG.contains(r1, r1) == true
# Disjoint ring with one "inside" of hole created => doesn't contain
@test GO.contains(r1, r2) == LG.contains(r1, r2) == false
# Disjoint ring with one "outside" of hole created => doesn't contain
@test GO.contains(r1, r3) == LG.contains(r1, r3) == false
# Rings share two sides and rest of sides don't touch -> doesn't contain
@test GO.contains(r1, r4) == LG.contains(r1, r4) == false
# Ring shares all edges with other ring, plus an extra loop -> contains
@test GO.contains(r5, r1) == LG.contains(r5, r1) == true
# Rings share just one vertex  -> doesn't contain
@test GO.contains(r6, r1) == LG.contains(r6, r1) == false
# Rings cross over one another -> doesn't contain
@test GO.contains(r7, r1) == LG.contains(r7, r1) == false

# Ring on bounday of polygon -> doesn't contain
@test GO.contains(p1, r4) == LG.contains(p1, r4) == false
# Ring on boundary and cutting through polygon -> contains
@test GO.contains(p1, r1) == LG.contains(p1, r1) == true
# Ring on hole boundary -> doesn't contain
@test GO.contains(p1, r2) == LG.contains(p1, r2) == false
# Ring touches polygon at one vertex -> doesn't contain
@test GO.contains(p1, r6) == LG.contains(p1, r6) == false
# Ring crosses through polygon -> doesn't contain
@test GO.contains(p1, r7) == LG.contains(p1, r7) == false
# Ring inside polygon -> contains
@test GO.contains(p1, r8) == LG.contains(p1, r8) == true
# Ring outside -> doesn't contain
@test GO.contains(p1, r9) == LG.contains(p1, r9) == false
# Ring inside polygon and shares holes edge -> contains
@test GO.contains(p1, r10) == LG.contains(p1, r10) == true
# Ring inside of polygon hole -> doesn't contain
@test GO.contains(p1, r11) == LG.contains(p1, r11) == false

# # Polygon and Geometry

# Polygon with holes in polygon without holes -> contains
@test GO.contains(p2, p1) == LG.contains(p2, p1) == true
# Polygon without holes in polygon with holes -> doesn't contain
@test GO.contains(p1, p2) == LG.contains(p1, p2) == false
# Polygon is the same as other poylgon hole -> doesn't contain
@test GO.contains(p1, p3) == LG.contains(p1, p3) == false
# Polygon touches other polygon by vertex -> doesn't contain
@test GO.contains(p1, p4) == LG.contains(p1, p4) == false
# Polygon outside of other polygon -> doesn't contain
@test GO.contains(p1, p5) == LG.contains(p1, p5) == false
# Polygon inside of hole -> doesn't contain
@test GO.contains(p1, p6) == LG.contains(p1, p6) == false
# Polygon overlaps other polygon -> doesn't contain
@test GO.contains(p1, p7) == LG.contains(p1, p7) == false
# Polygon with hole inside polygon with hole (holes nested) -> contains
@test GO.contains(p1, p10) == LG.contains(p1, p10) == true

# # Multi-geometries and collections

mp1 = LG.MultiPolygon([p1, p3, p11])
mp2 = LG.MultiPolygon([p2, p5])
c1 = LG.GeometryCollection([r1, r8, pt5])
# Multipolygon plugs all holes -> doesn't contain
@test GO.contains(mp1, p2) == LG.contains(mp1, p2) == false
# Polygon is one of the multipolygons -> contains
@test GO.contains(mp2, p2) == LG.contains(mp2, p2) == true
# Polygon touches one of the multipolygons -> doesn't contain
@test GO.contains(mp2, p4) == LG.contains(mp2, p4) == false
# Polygon contains all multipolygon elements
@test GO.contains(p9, mp1) == LG.contains(p9, mp1) == true
# Polygon contains all collection elements
@test GO.contains(p1, c1) == LG.contains(p1, c1) == true
# Collection doesn't contain all multipolygon elements
@test GO.contains(c1, mp1) == LG.contains(c1, mp1) == false