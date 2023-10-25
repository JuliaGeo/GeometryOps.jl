import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

r1 = LG.LinearRing([[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]])
r2 = LG.LinearRing([[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]])
GO.equals(r1, r1)
