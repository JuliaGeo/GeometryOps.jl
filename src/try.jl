import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

p2 = LG.Point([0.0, 1.0])
mp3 = LG.MultiPoint([p2])
GO.equals(p2, mp3)
