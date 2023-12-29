# # Point and Geometry

# Same point -> covered by
@test GO.coveredby(pt1, pt1) == LG.coveredby(pt1, pt1) == true
# Different point -> not covered by
@test GO.coveredby(pt1, pt2) == LG.coveredby(pt1, pt2) == false

# Point on line endpoint -> covered by
@test GO.coveredby(pt1, l1) == LG.coveredby(pt1, l1) == true
# Point outside line -> not covered by
@test GO.coveredby(pt2, l1) == LG.coveredby(pt2, l1) == false
# Point on line segment -> covered by
@test GO.coveredby(pt3, l1) == LG.coveredby(pt3, l1) == true
# Point on line vertex between segments -> covered by
@test GO.coveredby(pt4, l1) == LG.coveredby(pt4, l1) == true
# line cannot be covered by a point -> not covered by
@test GO.coveredby(l1, pt3) == LG.coveredby(l1, pt3) == false

# Point on ring endpoint -> covered by
@test GO.coveredby(pt1, r1) == LG.coveredby(pt1, r1) == true
# Point outside ring -> isn't covered by
@test GO.coveredby(pt2, r1) == LG.coveredby(pt2, r1) == false
# Point on ring segment -> covered by
@test GO.coveredby(pt3, r1) == LG.coveredby(pt3, r1) == true
# Point on ring vertex between segments -> covered by
@test GO.coveredby(pt4, r1) == LG.coveredby(pt4, r1) == true
# Ring cannot be covered by a point -> isn't covered by
@test GO.coveredby(r1, pt3) == LG.coveredby(r1, pt3) == false

# Point on vertex of polygon --> covered
@test GO.coveredby(pt1, p1) == LG.coveredby(pt1, p1) == true
# Point outside of polygon's external ring -> not covered by
@test GO.coveredby(pt2, p1) == LG.coveredby(pt2, p1) == false
# Point on polygon's edge -> covered by
@test GO.coveredby(pt4, p1) == LG.coveredby(pt4, p1) == true
# Point inside of polygon -> covered by
@test GO.coveredby(pt5, p1) == LG.coveredby(pt5, p1) == true
# Point on hole edge -> covered by
@test GO.coveredby(pt6, p1) == LG.coveredby(pt6, p1) == true
# Point inside of polygon hole -> not covered by
@test GO.coveredby(pt7, p1) == LG.coveredby(pt7, p1) == false
# Polygon can't be covered by a polygon -> not covered by
@test GO.coveredby(p1, pt5) == LG.coveredby(p1, pt5) == false

# # Line and Geometry

# Same line -> covered by
@test GO.coveredby(l1, l1) == LG.coveredby(l1, l1) == true
# Line overlaps line edge and endpoint -> covered by
@test GO.coveredby(l2, l1) == LG.coveredby(l2, l1) == true
# Line overlaps with one edge and is outside of other edge -> isn't covered by
@test GO.coveredby(l3, l1) == LG.coveredby(l3, l1) == false
# Line segments both within other line segments -> covered by
@test GO.coveredby(l4, l1) == LG.coveredby(l4, l1) == true
# Line segments connect at endpoint -> isn't covered by
@test GO.coveredby(l5, l1) == LG.coveredby(l5, l1) == false
# Line segments don't touch -> isn't covered by
@test GO.coveredby(l6, l1) == LG.coveredby(l6, l1) == false
# Line segments cross -> isn't covered by
@test GO.coveredby(l7, l1) == LG.coveredby(l7, l1) == false
# Line segments cross and go over and out -> isn't covered by
@test GO.coveredby(l8, l1) == LG.coveredby(l8, l1) == false
# Line segments cross and overlap on endpoint -> isn't covered by
@test GO.coveredby(l9, l1) == LG.coveredby(l9, l1) == false

# Line is within linear ring -> covered by
@test GO.coveredby(l1, r1) == LG.coveredby(l1, r1) == true
# Line covers one edge of linera ring and has segment outside -> isn't covered by
@test GO.coveredby(l3, r1) == LG.coveredby(l3, r1) == false
# Line and linear ring are only connected at vertex -> isn't covered by
@test GO.coveredby(l5, r1) == LG.coveredby(l5, r1) == false
# Line and linear ring are disjoint -> isn't covered by
@test GO.coveredby(l6, r1) == LG.coveredby(l6, r1) == false
# Line crosses through two ring edges -> isn't covered by
@test GO.coveredby(l7, r1) == LG.coveredby(l7, r1) == false
# Line crosses through two ring edges and touches third edge -> isn't covered by
@test GO.coveredby(l8, r1) == LG.coveredby(l8, r1) == false
# Line is equal to linear ring -> covered by
@test GO.coveredby(l10, r1) == LG.coveredby(l10, r1) == true
# Line covers linear ring and then has extra segment -> isn't covered by
@test GO.coveredby(l11, r1) == LG.coveredby(l11, r1) == false

# Line on polygon edge -> coveredby
@test GO.coveredby(l1, p1) == LG.coveredby(l1, p1) == true
# Line on polygon edge and extending beyond polygon edge -> not coveredby
@test GO.coveredby(l3, p1) == LG.coveredby(l3, p1) == false
# Line outside polygon connected by an vertex -> not coveredby
@test GO.coveredby(l5, p1) == LG.coveredby(l5, p1) == false
# Line through polygon cutting to outside -> not coveredby
@test GO.coveredby(l7, p1) == LG.coveredby(l7, p1) == false
# Line inside of polygon -> coveredby
@test GO.coveredby(l12, p1) == LG.coveredby(l12, p1) == true
# Line outside of polygon -> not coveredby
@test GO.coveredby(l13, p1) == LG.coveredby(l13, p1) == false
# Line in polygon hole -> not coveredby
@test GO.coveredby(l14, p1) == LG.coveredby(l14, p1) == false
# Line outside crown polygon but touching edges -> not coveredby
@test GO.coveredby(l15, p8) == LG.coveredby(l15, p8) == false
# Line within crown polygon but touching edges -> not coveredby
@test GO.coveredby(l15, p9) == LG.coveredby(l15, p9) == true

# # Ring and Geometry

# Line is within linear ring -> not coveredby
@test GO.coveredby(r1, l1) == LG.coveredby(r1, l1) == false
# Line covers one edge of linera ring and has segment outside -> not coveredby
@test GO.coveredby(r1, l3) == LG.coveredby(r1, l3) == false
# Line and linear ring are only connected at vertex -> not coveredby
@test GO.coveredby(r1, l5) == LG.coveredby(r1, l5) == false
# Line and linear ring are disjoint -> not coveredby
@test GO.coveredby(r1, l6) == LG.coveredby(r1, l6) == false
# Line crosses through two ring edges -> not coveredby
@test GO.coveredby(r1, l7) == LG.coveredby(r1, l7) == false
# Line crosses through two ring edges and touches third edge -> not coveredby
@test GO.coveredby(r1, l8) == LG.coveredby(r1, l8) == false
# Line is equal to linear ring -> coveredby
@test GO.coveredby(r1, l10) == LG.coveredby(r1, l10) == true
# Line covers linear ring and then has extra segment -> coveredby
@test GO.coveredby(r1, l11) == LG.coveredby(r1, l11) == true

# Same ring -> coveredby
@test GO.coveredby(r1, r1) == LG.coveredby(r1, r1) == true
# Disjoint ring with one "inside" of hole created => not coveredby
@test GO.coveredby(r2, r1) == LG.coveredby(r2, r1) == false
# Disjoint ring with one "outside" of hole created => not coveredby
@test GO.coveredby(r3, r1) == LG.coveredby(r3, r1) == false
# Rings share two sides and rest of sides don't touch -> not coveredby
@test GO.coveredby(r4, r1) == LG.coveredby(r4, r1) == false
# Ring shares all edges with other ring, plus an extra loop -> coveredby
@test GO.coveredby(r1, r5) == LG.coveredby(r1, r5) == true
# Rings share just one vertex  -> not coveredby
@test GO.coveredby(r1, r6) == LG.coveredby(r1, r6) == false
# Rings cross over one another -> not coveredby
@test GO.coveredby(r1, r7) == LG.coveredby(r1, r7) == false

# Ring on bounday of polygon -> coveredby
@test GO.coveredby(r4, p1) == LG.coveredby(r4, p1) == true
# Ring on boundary and cutting through polygon -> coveredby
@test GO.coveredby(r1, p1) == LG.coveredby(r1, p1) == true
# Ring on hole boundary -> coveredby
@test GO.coveredby(r2, p1) == LG.coveredby(r2, p1) == true
# Ring touches polygon at one vertex -> not coveredby
@test GO.coveredby(r6, p1) == LG.coveredby(r6, p1) == false
# Ring crosses through polygon -> not coveredby
@test GO.coveredby(r7, p1) == LG.coveredby(r7, p1) == false
# Ring inside polygon -> coveredby
@test GO.coveredby(r8, p1) == LG.coveredby(r8, p1) == true
# Ring outside -> not coveredby
@test GO.coveredby(r9, p1) == LG.coveredby(r9, p1) == false
# Ring inside polygon and shares holes edge -> coveredby
@test GO.coveredby(r10, p1) == LG.coveredby(r10, p1) == true
# Ring inside of polygon hole -> not coveredby
@test GO.coveredby(r11, p1) == LG.coveredby(r11, p1) == false

# # Polygon and Geometry

# Polygon with holes in polygon without holes -> coveredby
@test GO.coveredby(p1, p2) == LG.coveredby(p1, p2) == true
# Polygon without holes in polygon with holes -> not coveredby
@test GO.coveredby(p2, p1) == LG.coveredby(p2, p1) == false
# Polygon is the same as other poylgon hole -> not coveredby
@test GO.coveredby(p3, p1) == LG.coveredby(p3, p1) == false
# Polygon touches other polygon by vertex -> not coveredby
@test GO.coveredby(p4, p1) == LG.coveredby(p4, p1) == false
# Polygon outside of other polygon -> not coveredby
@test GO.coveredby(p5, p1) == LG.coveredby(p5, p1) == false
# Polygon inside of hole -> not coveredby
@test GO.coveredby(p6, p1) == LG.coveredby(p6, p1) == false
# Polygon overlaps other polygon -> not coveredby
@test GO.coveredby(p7, p1) == LG.coveredby(p7, p1) == false
# Polygon with hole inside polygon with hole (holes nested) -> coveredby
@test GO.coveredby(p10, p1) == LG.coveredby(p10, p1) == true

# # Multi-geometries and collections
mpt1 = LG.MultiPoint([pt1, pt4, pt5])
ml1 = LG.MultiLineString([l1, l2, l3])
c1 = LG.GeometryCollection([l1, pt1, ])
# Three points all in polygon
@test GO.coveredby(mpt1, p1) == LG.coveredby(mpt1, p1) == true
# Polygon can't be covered by multipoints
@test GO.coveredby(p1, mpt1) == LG.coveredby(p1, mpt1)
# Three lines not all covered by line
@test GO.coveredby(ml1, l1) == LG.coveredby(ml1, l1) == false
