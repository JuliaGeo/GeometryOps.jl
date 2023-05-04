using Test

import GeoInterface as GI
import GeometryOps as GO

geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                   GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])

flipped_geom = GO.map(GI.PointTrait, geom) do p
    (GI.y(p), GI.x(p))
end

@test flipped_geom == GI.Polygon([GI.LinearRing([(2, 1), (4, 3), (6, 5), (2, 1)]), 
                                  GI.LinearRing([(4, 3), (6, 5), (7, 6), (4, 3)])])

multipoints = GeometryOps.map(GI.LinearRingTrait, geom) do poly
    GI.MultiPoint(GI.getpoint(poly))
end

@test multipoints == [GI.MultiPoint([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                      GI.MultiPoint([(3, 4), (5, 6), (6, 7), (3, 4)])]

lengths = GeometryOps.map(GI.LinearRingTrait, geom) do poly
    GI.npoint(poly)
end

@test lengths == [4, 4]
