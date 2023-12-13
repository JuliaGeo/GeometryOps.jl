import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

r1 = LG.LinearRing([[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0]])
p3 = LG.Polygon([[[0.0, 0.0], [1.0, 0.0], [0.0, 0.2], [0.0, 0.0]]])

GO.disjoint(r1, p3)
    