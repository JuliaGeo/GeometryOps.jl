const in_geom = 1#(true, false)
const on_geom = -1#(false, true)
const not_in_on_geom = 0#(false, false)

warn_msg = "Linestring isn't closed. Point cannot be 'in' linestring."
open_string = LG.LineString([
    [0.0, 0.0], [1.0, 1.0], [3.0, 1.25], [2.0, 3.0], [-1.0, 2.75]
])
closed_string = LG.LineString([
    [0.0, 0.0], [1.0, 1.0], [3.0, 1.25], [2.0, 3.0], [-1.0, 2.75], [0.0, 0.0]
])
rect_ring = LG.LinearRing([
    [-20.0, 0.0], [-20.0, -10.0], [-5.0, -10.0], [-5.0, 0.0], [-20.0, 0.0]
])
tri_ring = LG.LinearRing([
    [5.0, 0.0], [7.0, 1.995], [9.0, -1.0], [5.0, 0.0]
])
concave_out_spikes = LG.LinearRing([
    [0.0, 0.0], [0.0, 10.0], [20.0, 10.0], [20.0, 0.0], [15.0, -5.0],
    [10.0, 0.0], [5.0, -5.0], [0.0, 0.0],
])
concave_in_spikes = LG.LinearRing([
    [0.0, 0.0], [0.0, 10.0], [20.0, 10.0], [20.0, 0.0], [15.0, 5.0],
    [10.0, 0.0], [5.0, 5.0], [0.0, 0.0],
])
rect_poly = LG.Polygon([[
    [0.0, 0.0], [0.0, 10.0], [10.0, 10.0], [10.0, 0.0], [0.0, 0.0]
]])
diamond_poly = LG.Polygon([[
    [0.0, 0.0], [-5.0, 5.0], [0.0, 10.0], [5.0, 5.0], [0.0, 0.0],
]])
trap_with_hole = LG.Polygon([
    [[-10.0, 0.0], [-8.0, 5.0], [8.0, 5.0], [10.0, 0.0], [-10.0, 0.0]],
    [[-5.0, 2.0], [-5.0, 4.0], [-2.0, 4.0], [-2.0, 2.0], [-5.0, 2.0]]
])
concave_a = GI.Polygon([[
    (1.2938349167338743, -3.175128530227131),
    (-2.073885870841754, -1.6247711001754137),
    (-5.787437985975053, 0.06570713422599561),
    (-2.1308128111898093, 5.426689675486368),
    (2.3058074184797244, 6.926652158268195),
    (1.2938349167338743, -3.175128530227131),
]])
concave_b = GI.Polygon([[
    (-2.1902469793743924, -1.9576242117579579),
    (-4.726006206053999, 1.3907098941556428),
    (-3.165301985923147, 2.847612825874245),
    (-2.5529280962099428, 4.395492123980911),
    (0.5677700216973937, 6.344638314896882),
    (3.982554842356183, 4.853519613487035),
    (5.251193948893394, 0.9343031382106848),
    (5.53045582244555, -3.0101433691361734),
    (-2.1902469793743924, -1.9576242117579579),
]])


@testset "Point in Geom" begin
    # Line Strings
    @test (@test_logs (:warn, warn_msg) GO.point_in_geom((-12.0, -0.5), open_string)) == not_in_on_geom
    
    @test GO.point_in_geom((0.5, 0.5), closed_string) == on_geom
    @test GO.point_in_geom((2.0, 1.25), closed_string) == in_geom
    @test GO.point_in_geom((4.0, 0.0), closed_string) == not_in_on_geom

    # Linear Rings
    @test GO.point_in_geom((-12.0, -0.5), rect_ring) == in_geom
    @test GO.point_in_geom((20.0, 0.0), rect_ring) == not_in_on_geom
    @test GO.point_in_geom((-5.0, -10.0), rect_ring) == on_geom

    @test GO.point_in_geom((6.0, 0.0), tri_ring) == in_geom
    @test GO.point_in_geom((6.0, -0.25), tri_ring) == on_geom
    @test GO.point_in_geom((7.0, 2.0), tri_ring) == not_in_on_geom

    # Convex polygons
    @test GO.point_in_polygon((0.0, 0.0), rect_poly) == on_geom
    @test GO.point_in_polygon((0.0, 5.0), rect_poly) == on_geom
    @test GO.point_in_polygon((5.0, 10.0), rect_poly) == on_geom
    @test GO.point_in_polygon((2.5, 2.5), rect_poly) == in_geom
    @test GO.point_in_polygon((9.99, 9.99), rect_poly) == in_geom
    @test GO.point_in_polygon((20.0, 20.0), rect_poly) == not_in_on_geom

    @test GO.point_in_polygon((0.0, 0.0), diamond_poly) == on_geom
    @test GO.point_in_polygon((-2.5, 2.5), diamond_poly) == on_geom
    @test GO.point_in_polygon((2.5, 2.5), diamond_poly) == on_geom
    @test GO.point_in_polygon((0.0, 5.0), diamond_poly) == in_geom
    @test GO.point_in_polygon((4.99, 5.0), diamond_poly) == in_geom
    @test GO.point_in_polygon((20.0, 20.0), diamond_poly) == not_in_on_geom

    @test GO.point_in_polygon((-10.0, 0.0), trap_with_hole) == on_geom
    @test GO.point_in_polygon((-5.0, 2.0), trap_with_hole) == on_geom
    @test GO.point_in_polygon((-5.0, 3.0), trap_with_hole) == on_geom
    @test GO.point_in_polygon((-9.0, 0.01), trap_with_hole) == in_geom
    @test GO.point_in_polygon((-4.0, 3.0), trap_with_hole) == not_in_on_geom
    @test GO.point_in_polygon((20.0, 20.0), trap_with_hole) == not_in_on_geom
    
    # Concave polygons
    pt = (-2.1902469793743924, -1.9576242117579579)
    @test GO.point_in_polygon(pt, concave_a) == not_in_on_geom
    @test GO.point_in_polygon(pt, concave_b) == on_geom
    @test GO.point_in_polygon((0.0, 0.0), concave_a) == in_geom
    @test GO.point_in_polygon((0.0, 0.0), concave_b) == in_geom
end

@testset "Line in Geom" begin
    # Line Strings
    @test (@test_logs (:warn, warn_msg) GO.line_in_geom(
        LG.LineString([[0.0, 0.0], [0.0, 1.0]]),
        open_string,
    )) == not_in_on_geom

    # On the edge
    @test GO.line_in_geom(
        LG.LineString([[0.25, 0.25], [0.5, 0.5]]),
        closed_string
    ) == on_geom
    # Inside
    @test GO.line_in_geom(
        LG.LineString([[1.0, 1.25], [2, 1.3], [2.9, 1.25]]),
        closed_string
    ) == in_geom
    # Inside to outside
    @test GO.line_in_geom(
        LG.LineString([[1.0, 0.99], [2, 1.3], [2.9, 1.25]]),
        closed_string
    ) == not_in_on_geom
    # Outside of geom
    @test GO.line_in_geom(
        LG.LineString([[0.0, 0.0], [3.0, 1.25]]),
        closed_string
    ) == not_in_on_geom

    # Rings
    # Same geometry, sharing all edges
    @test GO.line_in_geom(
        LG.LineString([[-20.0, 0.0], [-20.0, -10.0], [-5.0, -10.0], [-5.0, 0.0], [-20.0, 0.0]]),
        rect_ring
    ) == on_geom

    # Same geometry, sharing all but last edge, which is inside
    @test GO.line_in_geom(
        LG.LineString([[-20.0, 0.0], [-20.0, -10.0], [-5.0, -10.0], [-5.0, 0.0], [-15.0, -5.0]]),
        rect_ring
    ) == on_geom
    # Within
    @test GO.line_in_geom(
        LG.LineString([[-10.0, -1.0], [-10.0, -9.0], [-7.0, -9.0], [-7.0, -1.0], [-10.0, -1.0]]),
        rect_ring
    ) == in_geom
    # Passing in to out
    @test GO.line_in_geom(
        LG.LineString([[-20.0, 0.0], [-20.0, -10.0], [-5.0, -10.0], [-5.0, 0.0], [-30.0, 0.0]]),
        rect_ring
    ) == not_in_on_geom

    horizontal_line = LG.LineString([[0.0, 0.0], [20.0, 0.0]])
    @test GO.line_in_geom(horizontal_line, concave_out_spikes) == on_geom
    @test GO.line_in_geom(horizontal_line, concave_in_spikes) == not_in_on_geom
    @test GO.line_in_geom(
        LG.LineString([[10.0, 0.0], [10.0, -0.0001]]),
        concave_out_spikes
    ) == not_in_on_geom

    # Polygons


end

@testset "Ring in Geom" begin

end

@testset "Polygon in Geom" begin

end