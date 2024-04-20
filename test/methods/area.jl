using Test
import GeoInterface as GI, GeometryOps as GO, LibGEOS as LG

pt = LG.Point([0.0, 0.0])
empty_pt = LG.readgeom("POINT EMPTY")
mpt = LG.MultiPoint([[0.0, 0.0], [1.0, 0.0]])
empty_mpt = LG.readgeom("MULTIPOINT EMPTY")
l1 = LG.LineString([[0.0, 0.0], [0.5, 0.5], [1.0, 0.5]])
empty_l = LG.readgeom("LINESTRING EMPTY")
ml1 = LG.MultiLineString([[[0.0, 0.0], [0.5, 0.5], [1.0, 0.5]], [[0.0, 0.0], [0.0, 0.1]]])
empty_ml = LG.readgeom("MULTILINESTRING EMPTY")
empty_l = LG.readgeom("LINESTRING EMPTY")
r1 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 0.0]])
empty_r = LG.readgeom("LINEARRING EMPTY")
p1 = LG.Polygon([
    [[10.0, 0.0], [30.0, 0.0], [30.0, 20.0], [10.0, 20.0], [10.0, 0.0]],
])
p2 = LG.Polygon([
    [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
    [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]],
])
p3 = LG.Polygon([
    [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
    [[15.0, 1.0], [25.0, 1.0], [25.0, 11.0], [15.0, 11.0], [15.0, 1.0]],
])
p4 = LG.Polygon([
    [
        [0.0, 5.0], [2.0, 2.0], [5.0, 2.0], [2.0, -2.0], [5.0, -5.0],
        [0.0, -2.0], [-5.0, -5.0], [-2.0, -2.0], [-5.0, 2.0], [-2.0, 2.0],
        [0.0, 5.0],
    ],
])
empty_p = LG.readgeom("POLYGON EMPTY")
mp1 = LG.MultiPolygon([p2, p4])
empty_mp = LG.readgeom("MULTIPOLYGON EMPTY")
c = LG.GeometryCollection([p1, p2, r1, empty_l])
empty_c = LG.readgeom("GEOMETRYCOLLECTION EMPTY")

# Points, lines, and rings have zero area
@test GO.area(pt) == GO.signed_area(pt) == LG.area(pt) == 0
@test GO.area(empty_pt) == LG.area(empty_pt) == 0
@test GO.area(pt) isa Float64
@test GO.signed_area(pt, Float32) isa Float32
@test GO.signed_area(pt) isa Float64
@test GO.area(pt, Float32) isa Float32
@test GO.area(mpt) == GO.signed_area(mpt) == LG.area(mpt) == 0
@test GO.area(empty_mpt) == LG.area(empty_mpt) == 0
@test GO.area(l1) == GO.signed_area(l1) == LG.area(l1) == 0
@test GO.area(empty_l) == LG.area(empty_l) == 0
@test GO.area(ml1) == GO.signed_area(ml1) == LG.area(ml1) == 0
@test GO.area(empty_ml) == LG.area(empty_ml) == 0
@test GO.area(r1) == GO.signed_area(r1) == LG.area(r1) == 0
@test GO.area(empty_r) == LG.area(empty_r) == 0

# Polygons have non-zero area
# CCW polygons have positive signed area
@test GO.area(p1) == GO.signed_area(p1) == LG.area(p1)
@test GO.signed_area(p1) > 0
# Float32 calculations
@test GO.area(p1) isa Float64
@test GO.area(p1, Float32) isa Float32
# CW polygons have negative signed area
a2 = LG.area(p2)
@test GO.area(p2) == a2
@test GO.signed_area(p2) == -a2
# Winding order of holes doesn't affect sign of signed area
@test GO.signed_area(p3) == -a2
# Concave polygon correctly calculates area
a4 = LG.area(p4)
@test GO.area(p4) == a4
@test GO.signed_area(p4) == -a4
# Empty polygon
@test GO.area(empty_p) == LG.area(empty_p) == 0
@test GO.signed_area(empty_p) == 0

# Multipolygon calculations work
@test GO.area(mp1) == a2 + a4
@test GO.area(mp1, Float32) isa Float32
# Empty multipolygon
@test GO.area(empty_mp) == LG.area(empty_mp) == 0


# Geometry collection summed area
@test GO.area(c) == LG.area(c)
@test GO.area(c, Float32) isa Float32
# Empty collection
@test GO.area(empty_c) == LG.area(empty_c) == 0