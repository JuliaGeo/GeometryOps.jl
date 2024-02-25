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


@testset "Coverage" begin
    cell_extremes = (0.0, 20.0, 0.0, 20.0)
    cell_poly = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)]])
    cell_area = 400.0

    # points, lines, curves

    # polygon is the same as the cell
    p1 = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)]])
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
    # non-convex polygon split into two pieces
    p8 =  GI.Polygon([[(-10.0, 15.0), (10.0, 15.0), (10.0, 12.0), (-5.0, 12.0), (10.0, 9.0), (10.0, 6.0), (-10.0, 6.0), (-10.0, 15.0)]])
    @test GO.coverage(p8, cell_extremes...) ≈ LG.area(LG.intersection(p8, cell_poly))
    # counter-clockwise polygon

    # polygon with a hole
    

    # function test_rand_polys(n)
    #     xmin = 0.0
    #     xmax = 20.0
    #     ymin = 0.0
    #     ymax = 20.0
    #     seed = Xoshiro(1999)

    #     cell_poly = LG.Polygon([[
    #         [xmin, ymin],
    #         [xmin, ymax],
    #         [xmax, ymax],
    #         [xmax, ymin],
    #         [xmin, ymin],
    #     ]])

    #     for i in 1:n
    #         println(i)
    #         x1 = rand() + rand(seed, 4:15)
    #         y1 = rand() + rand(seed, 4:15)
    #         nverts1 = rand(seed, 4:15)
    #         avg_radius1 = rand(seed, 4:15)
    #         irregularity1 = rand(seed) * 0.25
    #         spikiness1 = rand(seed) * 0.25

    #         coords1 = generate_random_poly(
    #             x1,
    #             y1,
    #             nverts1,
    #             avg_radius1,
    #             irregularity1,
    #             spikiness1,
    #             seed,
    #         )
    #         poly1 = LG.Polygon(coords1)

    #         if LG.isValid(poly1)
    #             # Coords 1 with Grid Cell
    #             a1_sub = GO.coverage(poly1, xmin, xmax, ymin, ymax)
    #             a1_lib = LG.area(LG.intersection(cell_poly, poly1))
    #             @assert isapprox(a1_sub, a1_lib, atol = 1e-12) "$poly1"
    #         end
    #     end
    #     return
    # end
    # test_rand_polys(25)

end
