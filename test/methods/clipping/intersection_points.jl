using Test
import GeoInterface as GI, GeometryOps as GO, LibGEOS as LG

l1 = GI.LineString([(90000.0, 1000.0), (90000.0, 22500.0), (95000.0, 22500.0), (95000.0, 1000.0), (90000.0, 1000.0)])
l2 = GI.LineString([(90000.0, 7500.0), (107500.0, 27500.0), (112500.0, 27500.0), (95000.0, 7500.0), (90000.0, 7500.0)])
l3 = GI.LineString([(90000.0, 90000.0), (90000.0, 105000.0), (105000.0, 105000.0), (105000.0, 90000.0), (90000.0, 90000.0)])
l4 = GI.LineString([(-98000.0, 90000.0), (-98000.0, 105000.0), (98000.0, 105000.0), (98000.0, 90000.0), (-98000.0, 90000.0)])
l5 = GI.LineString([(19999.999, 25000.0), (19999.999, 29000.0), (39999.998999999996, 29000.0), (39999.998999999996, 25000.0), (19999.999, 25000.0)])
l6 = GI.LineString([(0.0, 25000.0), (0.0, 29000.0), (20000.0, 29000.0), (20000.0, 25000.0), (0.0, 25000.0)])

p1, p2 = GI.Polygon([l1]), GI.Polygon([l2])

# Three intersection points
LG_l1_l2_mp = GI.MultiPoint(collect(GI.getpoint(LG.intersection(l1, l2))))
@test_implementations GO.equals(GI.MultiPoint(GO.intersection_points($l1, $l2)), LG_l1_l2_mp)

# Four intersection points with large intersection
LG_l3_l4_mp = GI.MultiPoint(collect(GI.getpoint(LG.intersection(l3, l4))))
@test_implementations GO.equals(GI.MultiPoint(GO.intersection_points($l3, $l4)), LG_l3_l4_mp)

# Four intersection points with very small intersection
LG_l5_l6_mp = GI.MultiPoint(collect(GI.getpoint(LG.intersection(l5, l6))))
@test_implementations GO.equals(GI.MultiPoint(GO.intersection_points($l5, $l6)), LG_l5_l6_mp)

# Test that intersection points between lines and polygons is equivalent
@test_implementations GO.equals(GI.MultiPoint(GO.intersection_points($p1, $p2)), GI.MultiPoint(GO.intersection_points($l1, $l2)))

# No intersection points between polygon and line
@test_implementations isempty(GO.intersection_points($p1, $l6))
