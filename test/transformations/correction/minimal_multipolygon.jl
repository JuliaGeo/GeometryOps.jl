using Test

import GeoInterface as GI
import GeometryOps as GO

p1 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
# (p1, p2) -> one polygon inside of the other, sharing an edge
p2 = GI.Polygon([[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 0.0)]])
# (p1, p3) -> polygons outside of one another, sharing an edge
p3 = GI.Polygon([[(1.0, 0.0), (2.0, 0.0), (2.0, -1.0), (1.0, -1.0), (1.0, 0.0)]])
# (p1, p4) -> polygons are completly disjoint (no holes)
p4 = GI.Polygon([[(1.0, -1.0), (2.0, -1.0), (2.0, -2.0), (1.0, -2.0), (1.0, -1.0)]])

mp1 = GI.MultiPolygon([p1])
mp2 = GI.MultiPolygon([p1, p1])
mp3 = GI.MultiPolygon([p1, p2, p3])
mp4 = GI.MultiPolygon([p1, p4])

fixed_mp1 = GO.MinimalMultiPolygon()(mp1)
@test GI.npolygon(fixed_mp1) == 1
@test GO.equals(GI.getpolygon(fixed_mp1, 1), p1)

fixed_mp2 = GO.MinimalMultiPolygon()(mp2)
@test GI.npolygon(fixed_mp2) == 1
@test GO.equals(GI.getpolygon(fixed_mp2, 1), p1)

fixed_mp3 = GO.MinimalMultiPolygon()(mp3)
@test GI.npolygon(fixed_mp3) == 1
@test GO.coveredby(p1, fixed_mp3) && GO.coveredby(p2, fixed_mp3) && GO.coveredby(p3, fixed_mp3)

fixed_mp4 = GO.MinimalMultiPolygon()(mp4)
@test GI.npolygon(fixed_mp4) == 2
@test (GO.equals(GI.getpolygon(fixed_mp4, 1), p1) && GO.equals(GI.getpolygon(fixed_mp4, 2), p4)) ||
    (GO.equals(GI.getpolygon(fixed_mp4, 2), p1) && GO.equals(GI.getpolygon(fixed_mp4, 1), p4))