import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

diamond = LG.Polygon([[
    [0.0, 0.0], [-5.0, 5.0], [0.0, 10.0], [5.0, 5.0], [0.0, 0.0],
]])
GO.point_in_polygon((-2.5, 2.5), diamond) == on_geom