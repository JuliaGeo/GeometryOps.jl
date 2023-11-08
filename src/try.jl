import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

closed_string = LG.LineString([
    [0.0, 0.0], [1.0, 1.0], [3.0, 1.25], [2.0, 3.0], [-1.0, 2.75], [0.0, 0.0]
])
GO.line_in_geom(
    LG.LineString([[0.0, 0.0], [3.0, 1.25]]),
    closed_string
) == not_in_on_geom
    