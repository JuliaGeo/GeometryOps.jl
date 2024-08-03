using Test
import GeoInterface as GI, 
    GeometryOps as GO, 
    LibGEOS as LG,
    ArchGDAL as AG

@testset "Lines/Rings" begin
    
    l1 = LG.LineString([[0.0, 0.0], [10.0, 0.0], [10.0, 10.0]])
    c1_from_LG = LG.centroid(l1)
    @test_implementations "Basic line string" l1 begin
        c1, len1 = GO.centroid_and_length(l1)
        @test c1[1] ≈ GI.x(c1_from_LG)
        @test c1[2] ≈ GI.y(c1_from_LG)
        @test len1 ≈ 20.0
    end

    # 
    l2 = LG.LineString([[0.0, 0.0], [2.5, -2.5], [-5.0, -3.0], [-4.0, 6.0], [10.0, 10.0], [12.0, -14.56]])
    c2_from_LG = LG.centroid(l2)
    @test_implementations "Spiral line string" l2 begin
        c2, len2 = GO.centroid_and_length(l2)
        @test c2[1] ≈ GI.x(c2_from_LG) 
        @test c2[2] ≈ GI.y(c2_from_LG)
        @test len2 ≈ 59.3090856788928

        # Test that non-closed line strings throw an error for centroid_and_area
        @test_throws AssertionError GO.centroid_and_area(l2)
    end

    # Basic linear ring - note that this still uses weighting by length
    r1 = LG.LinearRing([[0.0, 0.0], [3456.0, 7894.0], [6291.0, 1954.0], [0.0, 0.0]])
    c3_from_LG = LG.centroid(r1)

    @test_implementations "Basic linear ring" r1 begin
        c3 = GO.centroid(r1)
        @test c3[1] ≈ GI.x(c3_from_LG) 
        @test c3[2] ≈ GI.y(c3_from_LG)
    end
end

@testset "Polygons" begin
    # Basic rectangle
    p1 = AG.fromWKT("POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))")
    c1 = GO.centroid(p1)
    c1 .≈ (5, 5)
    @test GI.x(c1) ≈ 5
    @test GI.y(c1) ≈ 5

    # Concave c-shaped polygon
    p2 = LG.Polygon([[
        [11.0, 0.0], [11.0, 3.0], [14.0, 3.0], [14.0, 2.0], [12.0, 2.0],
        [12.0, 1.0], [14.0, 1.0], [14.0, 0.0], [11.0, 0.0],
    ]])
    c2, area2 = GO.centroid_and_area(p2)
    c2_from_LG = LG.centroid(p2)
    @test c2[1] ≈ GI.x(c2_from_LG)
    @test c2[2] ≈ GI.y(c2_from_LG)
    @test area2 ≈ LG.area(p2)

    # Randomly generated polygon with lots of sides
    p3 = LG.Polygon([[
        [14.567, 8.974], [13.579, 8.849], [12.076, 8.769], [11.725, 8.567],
        [11.424, 6.451], [10.187, 7.712], [8.187, 6.795], [8.065, 8.096],
        [6.827, 8.287], [7.1628, 8.9221], [5.8428, 10.468], [7.987, 11.734],
        [7.989, 12.081], [8.787, 11.930], [7.568, 13.926], [9.330, 13.340],
        [9.6817, 13.913], [10.391, 12.222], [12.150, 12.032], [14.567, 8.974],
    ]])
    c3, area3 = GO.centroid_and_area(p3)
    c3_from_LG = LG.centroid(p3)
    @test c3[1] ≈ GI.x(c3_from_LG)
    @test c3[2] ≈ GI.y(c3_from_LG)
    @test area3 ≈ LG.area(p3)

    # Polygon with one hole
    p4 = LG.Polygon([
        [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0], [0.0, 0.0]],
        [[2.3, 2.7], [2.5, 4.9], [4.1, 5.2], [4.2, 1.9], [2.3, 2.7]],
    ])
    c4, area4 = GO.centroid_and_area(p4)
    c4_from_LG = LG.centroid(p4)
    @test c4[1] ≈ GI.x(c4_from_LG)
    @test c4[2] ≈ GI.y(c4_from_LG)
    @test area4 ≈ LG.area(p4)

    # Polygon with two holes
    p5 = LG.Polygon([
        [[-10.0, -10.0], [-2.0, 0.0], [6.0, -10.0], [-10.0, -10.0]],
        [[-8.0, -8.0], [-8.0, -7.0], [-4.0, -7.0], [-4.0, -8.0], [-8.0, -8.0]],
        [[-3.0,-9.0], [3.0, -9.0], [3.0, -8.5], [-3.0, -8.5], [-3.0, -9.0]],
    ])
    c5 = GO.centroid(p5)
    c5_from_LG = LG.centroid(p5)
    @test c5[1] ≈ GI.x(c5_from_LG)
    @test c5[2] ≈ GI.y(c5_from_LG)

    # Same polygon as P5 but using a GeoInterface polygon
    p6 = GI.Polygon([
        [(-10.0, -10.0), (-2.0, 0.0), (6.0, -10.0), (-10.0, -10.0)],
        [(-8.0, -8.0), (-8.0, -7.0), (-4.0, -7.0), (-4.0, -8.0), (-8.0, -8.0)],
        [(-3.0, -9.0), (3.0, -9.0), (3.0, -8.5), (-3.0, -8.5), (-3.0, -9.0)],
    ])
    c6 = GO.centroid(p6)
    @test all(c5 .≈ c6)
end
@testset "MultiPolygons" begin
    # Combine polygons made above
    m1 = LG.MultiPolygon([
        [
            [[11.0, 0.0], [11.0, 3.0], [14.0, 3.0], [14.0, 2.0], [12.0, 2.0],
            [12.0, 1.0], [14.0, 1.0], [14.0, 0.0], [11.0, 0.0]],
        ],
        [
            [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0], [0.0, 0.0]],
            [[2.3, 2.7], [2.5, 4.9], [4.1, 5.2], [4.2, 1.9], [2.3, 2.7]],
        ]
    ])
    c1, area1 = GO.centroid_and_area(m1)
    c1_from_LG = LG.centroid(m1)
    @test c1[1] ≈ GI.x(c1_from_LG)
    @test c1[2] ≈ GI.y(c1_from_LG)
    @test area1 ≈ LG.area(m1)
end
