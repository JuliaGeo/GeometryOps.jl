pt1 = GI.Point((0.0, 0.0))
l1 = GI.Line([(0.0, 0.0), (0.0, 1.0)])

concave_coords = [(0.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 2.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)]
l2 = GI.LineString(concave_coords)
l3 = GI.LineString(concave_coords[1:(end - 1)])
r1 = GI.LinearRing(concave_coords)
r2 = GI.LinearRing(concave_coords[1:(end - 1)])
r3 = GI.LinearRing([(1.0, 1.0), (1.0, 1.5), (1.5, 1.5), (1.5, 1.0), (1.0, 1.0)])
concave_angles = [90.0, 270.0, 90.0, 90.0, 90.0, 90.0]

p1 = GI.Polygon([r3])
p2 = GI.Polygon([[(0.0, 0.0), (0.0, 4.0), (3.0, 0.0), (0.0, 0.0)]])
p3 = GI.Polygon([[(-3.0, -2.0), (0.0,0.0), (5.0, 0.0), (-3.0, -2.0)]])
p4 = GI.Polygon([r1])
p5 = GI.Polygon([r1, r3])

mp1 = GI.MultiPolygon([p2, p3])
c1 = GI.GeometryCollection([pt1, l2, p2])

# Points and lines
@test isempty(GO.angles(pt1))
@test isempty(GO.angles(l1))

# LineStrings and Linear Rings
@test all(isapprox.(GO.angles(l2), concave_angles, atol = 1e-3))
@test all(isapprox.(GO.angles(l3), concave_angles[2:(end - 1)], atol = 1e-3))
@test all(isapprox.(GO.angles(r1), concave_angles, atol = 1e-3))
@test all(isapprox.(GO.angles(r2), concave_angles, atol = 1e-3))

# Polygons
p2_angles = [90.0, 36.8699, 53.1301]
p3_angles = [19.6538, 146.3099, 14.0362]
@test all(isapprox.(GO.angles(p1), [90.0 for _ in 1:4], atol = 1e-3))
@test all(isapprox.(GO.angles(p2), p2_angles, atol = 1e-3))
@test all(isapprox.(GO.angles(p3), p3_angles, atol = 1e-3))
@test all(isapprox.(GO.angles(p4), concave_angles, atol = 1e-3))
@test all(isapprox.(GO.angles(p5), vcat(concave_angles, [270.0 for _ in 1:4]), atol = 1e-3))

# Multi-geometries
@test all(isapprox.(GO.angles(mp1), [p2_angles; p3_angles], atol = 1e-3))
@test all(isapprox.(GO.angles(c1), [concave_angles; p2_angles], atol = 1e-3))