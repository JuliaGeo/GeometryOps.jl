import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

p1 = LG.Polygon([[[0.0, 0.0], [0.5, 1.5], [2.5, -0.5], [0.0, 0.0]]])
a = GO.within(p1, p1)
    