@testset "Lines/Rings" begin
    # Basic line string
    l1 = LG.LineString([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0]])
    c1 = GO.centroid(l1)
    c1_from_LG = LG.centroid(l1)
    @test GI.x(c1) ≈ GI.x(c1_from_LG)  # I got this from 
    @test GI.y(c1) ≈ GI.y(c1_from_LG)

    # Spiral line string
    l2 = LG.LineString([[0.0, 0.0], [2.5, 2.5], [10.0, 10.0]])
    c2 = GO.centroid(21)
    c2_from_LG = LG.centroid(l2)
    @test GI.x(c2) ≈ GI.x(c2_from_LG)  # I got this from 
    @test GI.y(c2) ≈ GI.y(c2_from_LG)

    # Basic linear ring

    # Fancier linear ring

end
@testset "Polygons" begin
    # Basic rectangle
    p1 = AG.fromWKT("POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))")
    c1 = GO.centroid(p1)
    @test GI.x(c1) ≈ 5
    @test GI.y(c1) ≈ 5

    # Concave c-shaped polygon
    p2 = LG.Polygon([[
        [0.0, 0.0], [0.0, 3.0], [3.0, 3.0], [3.0, 2.0], [1.0, 2.0],
        [1.0, 1.0], [3.0, 1.0], [3.0, 0.0], [0.0, 0.0],
    ]])
    c2 = GO.centroid(p2)
    c2_from_LG = LG.centroid(p2)
    @test GI.x(c2) ≈ GI.x(c2_from_LG)  # I got this from 
    @test GI.y(c2) ≈ GI.y(c2_from_LG)

    # Randomly generated polygon with lots of sides

    # Polygon with one hole
    p4 = AG.fromWKT("POLYGON((0 0, 10 0, 10 10, 0 10, 0 0), (2.3 2.7, 2.5 4.9, 4.1 5.2, 1.9 4.2, 2.3 2.7))")
    c4 = GO.centroid(p4)
    c4_from_AG = AG.centroid(p4)
    @test GI.x(c4) ≈ GI.x(c4_from_AG)
    @test GI.y(c4) ≈ GI.y(c4_from_AG)

    # Polygon with two holes
    
end
@testset "MultiPolygons" begin
    # Combine poylgons made above
    
end