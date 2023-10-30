const in_geom = (true, false)
const on_geom = (false, true)
const not_in_on_geom = (false, false)


@testset "Point in Geom" begin
    # Convex polygons
    rect = LG.Polygon([[
        [0.0, 0.0], [0.0, 10.0], [10.0, 10.0], [10.0, 0.0], [0.0, 0.0]
    ]])
    @test GO.point_in_polygon((0.0, 0.0), rect) == on_geom
    @test GO.point_in_polygon((0.0, 5.0), rect) == on_geom
    @test GO.point_in_polygon((5.0, 10.0), rect) == on_geom
    @test GO.point_in_polygon((2.5, 2.5), rect) == in_geom
    @test GO.point_in_polygon((9.99, 9.99), rect) == in_geom
    @test GO.point_in_polygon((20.0, 20.0), rect) == not_in_on_geom

    diamond = LG.Polygon([[
        [0.0, 0.0], [-5.0, 5.0], [0.0, 10.0], [5.0, 5.0], [0.0, 0.0],
    ]])
    @test GO.point_in_polygon((0.0, 0.0), diamond) == on_geom
    @test GO.point_in_polygon((-2.5, 2.5), diamond) == on_geom
    @test GO.point_in_polygon((2.5, 2.5), diamond) == on_geom
    @test GO.point_in_polygon((0.0, 5.0), diamond) == in_geom
    @test GO.point_in_polygon((4.99, 5.0), diamond) == in_geom
    @test GO.point_in_polygon((20.0, 20.0), diamond) == not_in_on_geom

    trap_with_hole = LG.Polygon([
        [[-10.0, 0.0], [-8.0, 5.0], [8.0, 5.0], [10.0, 0.0], [-10.0, 0.0]],
        [[-5.0, 2.0], [-5.0, 4.0], [-2.0, 4.0], [-2.0, 2.0], [-5.0, 2.0]]
    ])
    @test GO.point_in_polygon((-10.0, 0.0), trap_with_hole) == on_geom
    @test GO.point_in_polygon((-5.0, 2.0), trap_with_hole) == on_geom
    @test GO.point_in_polygon((-5.0, 3.0), trap_with_hole) == on_geom
    @test GO.point_in_polygon((-9.0, 0.01), trap_with_hole) == in_geom
    @test GO.point_in_polygon((-4.0, 3.0), trap_with_hole) == not_in_on_geom
    @test GO.point_in_polygon((20.0, 20.0), trap_with_hole) == not_in_on_geom
    # Concave polygons
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
    pt = (-2.1902469793743924, -1.9576242117579579)
    @test GO.point_in_polygon(pt, concave_a) == (false, false)
    @test GO.point_in_polygon(pt, concave_b) == (false, true)
    @test GO.point_in_polygon((0.0, 0.0), concave_a) == (true, false)
    @test GO.point_in_polygon((0.0, 0.0), concave_b) == (true, false) 
end