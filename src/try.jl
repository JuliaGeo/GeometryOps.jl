import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG


r3 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [0.0, 0.2], [0.0, 0.0]])
l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [0.0, 0.1]])

GO.within(l2, r3)
    