@testset "Lines/Rings" begin
    # Line test intersects -----------------------------------------------------

    # Test for parallel lines
    l1 = GI.Line([(0.0, 0.0), (2.5, 0.0)])
    l2 = GI.Line([(0.0, 1.0), (2.5, 1.0)])
    @test !GO.intersects(l1, l2; meets = 0)
    @test !GO.intersects(l1, l2; meets = 1)
    @test isnothing(GO.intersection(l1, l2))

    # Test for non-parallel lines that don't intersect
    l1 = GI.Line([(0.0, 0.0), (2.5, 0.0)])
    l2 = GI.Line([(2.0, -3.0), (3.0, 0.0)])
    @test !GO.intersects(l1, l2; meets = 0)
    @test !GO.intersects(l1, l2; meets = 1)
    @test isnothing(GO.intersection(l1, l2))

    # Test for lines only touching at endpoint
    l1 = GI.Line([(0.0, 0.0), (2.5, 0.0)])
    l2 = GI.Line([(2.0, -3.0), (2.5, 0.0)])
    @test GO.intersects(l1, l2; meets = 0)
    @test !GO.intersects(l1, l2; meets = 1)
    @test all(GO.intersection(l1, l2) .≈ (2.5, 0.0))

    # Test for lines that intersect in the middle
    l1 = GI.Line([(0.0, 0.0), (5.0, 5.0)])
    l2 = GI.Line([(0.0, 5.0), (5.0, 0.0)])
    @test GO.intersects(l1, l2; meets = 0)
    @test GO.intersects(l1, l2; meets = 1)
    @test all(GO.intersection(l1, l2) .≈ (2.5, 2.5))

    # Line string test intersects ----------------------------------------------

    # Single element line strings crossing over each other
    l1 = LG.LineString([[5.5, 7.2], [11.2, 12.7]])
    l2 = LG.LineString([[4.3, 13.3], [9.6, 8.1]])
    @test GO.intersects(l1, l2; meets = 0)
    @test GO.intersects(l1, l2; meets = 1)
    go_inter = GO.intersection(l1, l2)
    lg_inter = LG.intersection(l1, l2)
    @test go_inter[1][1] .≈ GI.x(lg_inter)
    @test go_inter[1][2] .≈ GI.y(lg_inter)

    # Multi-element line strings crossing over on vertex
    l1 = LG.LineString([[0.0, 0.0], [2.5, 0.0], [5.0, 0.0]])
    l2 = LG.LineString([[2.0, -3.0], [3.0, 0.0], [4.0, 3.0]])
    @test GO.intersects(l1, l2; meets = 0)
    # TODO: Do we want this to be false? It is vertex of segment, not of whole line string
    @test !GO.intersects(l1, l2; meets = 1)
    go_inter = GO.intersection(l1, l2)
    @test length(go_inter) == 1
    lg_inter = LG.intersection(l1, l2)
    @test go_inter[1][1] .≈ GI.x(lg_inter)
    @test go_inter[1][2] .≈ GI.y(lg_inter)

    # Multi-element line strings crossing over with multiple intersections
    l1 = LG.LineString([[0.0, -1.0], [1.0, 1.0], [2.0, -1.0], [3.0, 1.0]])
    l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [3.0, 0.0]])
    @test GO.intersects(l1, l2; meets = 0)
    @test GO.intersects(l1, l2; meets = 1)
    go_inter = GO.intersection(l1, l2)
    @test length(go_inter) == 3
    lg_inter = LG.intersection(l1, l2)
    @test issetequal(
        Set(go_inter),
        Set(GO._tuple_point.(GI.getpoint(lg_inter)))
    )

    # Line strings far apart so extents don't overlap

    # Line strings close together that don't overlap

    # Line string with empty line string

    # Closed linear ring with open line string

    # Closed linear ring with closed linear ring

    # @test issetequal(
    #     Subzero.intersect_lines(l1, l2),
    #     Set([(0.5, -0.0), (1.5, 0), (2.5, -0.0)]),
    # )
    # l2 = [[[10., 10]]]
    # @test issetequal(
    #     Subzero.intersect_lines(l1, l2),
    #     Set{Tuple{Float64, Float64}}(),
    # )


end

@testset "Polygons" begin
    # Two polygons that intersect

    # Two polygons that don't intersect

    # Polygon that intersects with linestring

end

@testset "MultiPolygons" begin
    # Multi-polygon and polygon that intersect

    # Multi-polygon and polygon that don't intersect

    # Multi-polygon that intersects with linestring
    
end