using Test
import GeoInterface as GI, GeometryOps as GO, LibGEOS as LG
using GeoInterface
using GeoInterface.Extents

cell_extremes = (0.0, 20.0, 0.0, 20.0)
cell_extent = Extents.Extent(X = (0.0, 20.0), Y = (0.0, 20.0))
cell_poly = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)]])
cell_area = 400.0

# # Points, lines, curves
pt1 = (1.0, 0.0)
@test GO.coverage(pt1, cell_extremes...) == 0.0
mpt1 = GI.MultiPoint([pt1, (25.0, 0.0), (100.0, 50.0)])
@test GO.coverage(mpt1, cell_extremes...) == 0.0
l1 = GI.LineString([(-1.0, 5.0), (5.0, 10.0), (10.0, -2.0), (22, 5.0)])
@test GO.coverage(l1, cell_extremes...) == 0.0

# # Polygons

# polygon is the same as the cell (input is extent)
p1 = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)]])
@test GO.coverage(p1, cell_extent) == cell_area
# polygon is the same as the cell (input is min/max values)
@test GO.coverage(p1, cell_extremes...) == cell_area
# polygon is bigger than the cell
p2 = GI.Polygon([[(-10, -10.0), (-10.0, 30.0), (30.0, 30.0), (300.0, -10.0), (-10.0, -10.0)]])
@test GO.coverage(p2, cell_extremes...) == cell_area
# polygon is completly inside of cell
p3 = GI.Polygon([[(5.0, 5.0), (5.0, 15.0), (15.0, 15.0), (15.0, 5.0), (5.0, 5.0)]])
@test GO.coverage(p3, cell_extremes...) ≈ LG.area(LG.intersection(cell_poly, p3))
# polygon exits cell through one edge
p4 = GI.Polygon([[(5.0, 5.0), (5.0, 25.0), (15.0, 25.0), (15.0, 5.0), (5.0, 5.0)]])
@test GO.coverage(p4, cell_extremes...) ≈ LG.area(LG.intersection(cell_poly, p4))
p5 = GI.Polygon([[(5.0, 5.0), (5.0, 25.0), (25.0, 25.0), (25.0, 5.0), (5.0, 5.0)]])
@test GO.coverage(p5, cell_extremes...) ≈ LG.area(LG.intersection(cell_poly, p5))
# polygon exits cell through multiple edges (north and east)
p6 = GI.Polygon([[(20.8826, 6.4239), (15.9663, 2.3014), (8.6078, 2.0995), (2.6849, 6.4088),
    (0.8449, 12.7452), (3.0813, 19.1654), (9.1906, 23.2520), (15.5835, 22.9101),
    (20.9143, 18.5933), (20.8826, 6.4239)]])
@test GO.coverage(p6, cell_extremes...) ≈ LG.area(LG.intersection(p6, cell_poly))
# polygon exits cell through multiple edges (west and east)
p7 = GI.Polygon([[(-5.0, 10.0), (-5.0, 25.0), (25.0, 25.0), (25.0, 10.0), (-5.0, 10.0)]])
@test GO.coverage(p7, cell_extremes...) ≈ LG.area(LG.intersection(p7, cell_poly))
# non-convex polygon split into two pieces (same wall)
p8 =  GI.Polygon([[(-10.0, 15.0), (10.0, 15.0), (10.0, 12.0), (-5.0, 12.0), (-5.0, 9.0), (10.0, 9.0), (10.0, 6.0), (-10.0, 6.0), (-10.0, 15.0)]])
@test GO.coverage(p8, cell_extremes...) ≈ LG.area(LG.intersection(p8, cell_poly))
# counter-clockwise polygon
p9 =  GI.Polygon([[(-10.0, 15.0), (-10.0, 6.0), (10.0, 6.0), (10.0, 9.0), (-5.0, 9.0), (-5.0, 12.0), (10.0, 12.0), (10.0, 15.0), (-10.0, 15.0),]])
@test GO.coverage(p9, cell_extremes...) ≈ LG.area(LG.intersection(p9, cell_poly))
# polygon with a hole
p10 = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)],
    [(10.0, 10.0), (10.0, 15.0), (15.0, 15.0), (15.0, 10.0), (10.0, 10.0)],
    [(7.0, 7.0), (5.0, 7.0), (5.0, 5.0), (7.0, 7.0)],
])
@test GO.coverage(p10, cell_extremes...) ≈ LG.area(LG.intersection(p10, cell_poly))
# non-convex polygon split into two pieces (differnet walls)
p11 = GI.Polygon([[(-10.0, 15.0), (-10.0, -10.0), (15.0, -10.0), (15.0, 5.0), (10.0, 5.0),
    (10.0, -5.0), (-5.0, -5.0), (-5.0, 10.0), (5.0, 10.0), (5.0, 15.0), (-10.0, 15.0)]])
@test GO.coverage(p11, cell_extremes...) ≈ LG.area(LG.intersection(p11, cell_poly))

# # MultiPolygons
mp1 = GI.MultiPolygon([p1, p2]) == 2cell_area
gc1 = GI.GeometryCollection([pt1, l1, p1]) == cell_area