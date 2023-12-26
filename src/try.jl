import GeometryOps as GO
#import GeoInterface as GI
import LibGEOS as LG

l1 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])
l9 = LG.LineString([[0.0, 1.0], [0.0, -1.0], [1.0, 1.0]])

GO.crosses(l1, l9)
    