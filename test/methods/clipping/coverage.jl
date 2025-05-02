using Test
import GeoInterface as GI
import GeometryBasics as GB
import GeometryOps as GO
import LibGEOS as LG
using ..TestHelpers

cell_extremes = (0.0, 20.0, 0.0, 20.0)
cell_extent = GI.Extents.Extent(X = (0.0, 20.0), Y = (0.0, 20.0))
cell_poly = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)]])
cell_area = 400.0

# Basic test cases
pt1 = (1.0, 0.0)
mpt1 = GI.MultiPoint([pt1, (25.0, 0.0), (100.0, 50.0)])
l1 = GI.LineString([(-1.0, 5.0), (5.0, 10.0), (10.0, -2.0), (22.0, 5.0)])

@testset_implementations "Points, lines, curves" begin
    @test GO.coverage($pt1, cell_extremes...) == 0.0
    @test GO.coverage($mpt1, cell_extremes...) == 0.0
    @test GO.coverage($l1, cell_extremes...) == 0.0
end

# Basic polygon test cases
p1 = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(-10.0, -10.0), (-10.0, 30.0), (30.0, 30.0), (300.0, -10.0), (-10.0, -10.0)]])
p2b = GI.Polygon([[(-10, -10.0), (-10.0, 30.0), (30.0, 30.0), (300.0, -10.0), (-10.0, -10.0)]])
p3 = GI.Polygon([[(5.0, 5.0), (5.0, 15.0), (15.0, 15.0), (15.0, 5.0), (5.0, 5.0)]])
p4 = GI.Polygon([[(5.0, 5.0), (5.0, 25.0), (15.0, 25.0), (15.0, 5.0), (5.0, 5.0)]])
p5 = GI.Polygon([[(5.0, 5.0), (5.0, 25.0), (25.0, 25.0), (25.0, 5.0), (5.0, 5.0)]])
p6 = GI.Polygon([[(20.8826, 6.4239), (15.9663, 2.3014), (8.6078, 2.0995), (2.6849, 6.4088),
    (0.8449, 12.7452), (3.0813, 19.1654), (9.1906, 23.2520), (15.5835, 22.9101),
    (20.9143, 18.5933), (20.8826, 6.4239)]])
p7 = GI.Polygon([[(-5.0, 10.0), (-5.0, 25.0), (25.0, 25.0), (25.0, 10.0), (-5.0, 10.0)]])
p8 = GI.Polygon([[(-10.0, 15.0), (10.0, 15.0), (10.0, 12.0), (-5.0, 12.0), (-5.0, 9.0), (10.0, 9.0), (10.0, 6.0), (-10.0, 6.0), (-10.0, 15.0)]])
p9 = GI.Polygon([[(-10.0, 15.0), (-10.0, 6.0), (10.0, 6.0), (10.0, 9.0), (-5.0, 9.0), (-5.0, 12.0), (10.0, 12.0), (10.0, 15.0), (-10.0, 15.0),]])
p10 = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)],
    [(10.0, 10.0), (10.0, 15.0), (15.0, 15.0), (15.0, 10.0), (10.0, 10.0)],
    [(7.0, 7.0), (5.0, 7.0), (5.0, 5.0), (7.0, 7.0)],
])
p11 = GI.Polygon([[(-10.0, 15.0), (-10.0, -10.0), (15.0, -10.0), (15.0, 5.0), (10.0, 5.0),
    (10.0, -5.0), (-5.0, -5.0), (-5.0, 10.0), (5.0, 10.0), (5.0, 15.0), (-10.0, 15.0)]])

# New test cases for edge cases
# Point intersection
p12 = GI.Polygon([[(10.0, 10.0), (10.0, 10.0), (10.0, 10.0), (10.0, 10.0)]])
# Line along boundary
p13 = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (0.0, 0.0)]])
# Line through vertices
p14 = GI.Polygon([[(0.0, 0.0), (20.0, 20.0), (0.0, 0.0)]])
# Line through vertices with extra point
p15 = GI.Polygon([[(0.0, 0.0), (10.0, 10.0), (20.0, 20.0), (0.0, 0.0)]])
# Polygon touching at single vertex
p16 = GI.Polygon([[(20.0, 20.0), (30.0, 20.0), (30.0, 30.0), (20.0, 20.0)]])
# Polygon touching at edge
p17 = GI.Polygon([[(20.0, 0.0), (20.0, 20.0), (30.0, 20.0), (30.0, 0.0), (20.0, 0.0)]])
# Polygon with hole touching cell boundary
p18 = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)],
    [(0.0, 5.0), (5.0, 5.0), (5.0, 15.0), (0.0, 15.0), (0.0, 5.0)]])
# Polygon with hole completely inside cell
p19 = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)],
    [(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0), (5.0, 5.0)]])
# Polygon with hole partially outside cell
p20 = GI.Polygon([[(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0), (0.0, 0.0)],
    [(-5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (-5.0, 15.0), (-5.0, 5.0)]])

@testset_implementations "Polygons" [GI, GB, LG] begin
    # Basic polygon tests
    @test GO.coverage($p1, cell_extent) == cell_area
    @test GO.coverage($p1, cell_extremes...) == cell_area
    @test GO.coverage($p10, cell_extremes...) ≈ LG.area(LG.intersection($p10, cell_poly))
    @test GO.coverage($p2, cell_extremes...) == cell_area
    @test_implementations GO.coverage($p2b, cell_extremes...) == cell_area
    @test GO.coverage($p3, cell_extremes...) ≈ LG.area(LG.intersection(cell_poly, $p3))
    @test GO.coverage($p4, cell_extremes...) ≈ LG.area(LG.intersection(cell_poly, $p4))
    @test GO.coverage($p5, cell_extremes...) ≈ LG.area(LG.intersection(cell_poly, $p5))
    @test GO.coverage($p6, cell_extremes...) ≈ LG.area(LG.intersection($p6, cell_poly))
    @test GO.coverage($p7, cell_extremes...) ≈ LG.area(LG.intersection($p7, cell_poly))
    @test GO.coverage($p8, cell_extremes...) ≈ LG.area(LG.intersection($p8, cell_poly))
    @test GO.coverage($p9, cell_extremes...) ≈ LG.area(LG.intersection($p9, cell_poly))
    @test GO.coverage($p10, cell_extremes...) ≈ LG.area(LG.intersection($p10, cell_poly))
    @test GO.coverage($p11, cell_extremes...) ≈ LG.area(LG.intersection($p11, cell_poly))

    # Edge case tests
    @test GO.coverage($p12, cell_extremes...) == 0.0  # Point intersection
    @test GO.coverage($p13, cell_extremes...) == 0.0  # Line along boundary
    @test GO.coverage($p14, cell_extremes...) == 0.0  # Line through vertices
    @test GO.coverage($p15, cell_extremes...) == 0.0  # Line through vertices with extra point
    @test GO.coverage($p16, cell_extremes...) == 0.0  # Polygon touching at single vertex
    @test GO.coverage($p17, cell_extremes...) == 0.0  # Polygon touching at edge
    @test GO.coverage($p18, cell_extremes...) ≈ LG.area(LG.intersection($p18, cell_poly))  # Hole touching boundary
    @test GO.coverage($p19, cell_extremes...) ≈ LG.area(LG.intersection($p19, cell_poly))  # Hole inside cell
    @test GO.coverage($p20, cell_extremes...) ≈ LG.area(LG.intersection($p20, cell_poly))  # Hole partially outside
end

# Test with different cell sizes
small_cell = (0.0, 10.0, 0.0, 10.0)
large_cell = (0.0, 30.0, 0.0, 30.0)

@testset_implementations "Different cell sizes" [GI, GB, LG] begin
    @test GO.coverage($p1, small_cell...) == 100.0
    @test GO.coverage($p1, large_cell...) == cell_area
    @test GO.coverage($p3, small_cell...) ≈ LG.area(LG.intersection($p3, GI.Polygon([[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (0.0, 0.0)]])))
    @test GO.coverage($p3, large_cell...) ≈ LG.area(LG.intersection($p3, GI.Polygon([[(0.0, 0.0), (0.0, 30.0), (30.0, 30.0), (30.0, 0.0), (0.0, 0.0)]])))
end

# Test with non-square cells
rect_cell = (0.0, 30.0, 0.0, 10.0)
rect_cell_poly = GI.Polygon([[(0.0, 0.0), (0.0, 10.0), (30.0, 10.0), (30.0, 0.0), (0.0, 0.0)]])

@testset_implementations "Non-square cells" [GI, GB, LG] begin
    @test GO.coverage($p1, rect_cell...) ≈ LG.area(LG.intersection($p1, rect_cell_poly))
    @test GO.coverage($p3, rect_cell...) ≈ LG.area(LG.intersection($p3, rect_cell_poly))
    @test GO.coverage($p7, rect_cell...) ≈ LG.area(LG.intersection($p7, rect_cell_poly))
end

# Test with empty and degenerate polygons
empty_poly = GI.Polygon([])
degenerate_poly = GI.Polygon([[(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)]])

@testset_implementations "Empty and degenerate polygons" [GI, GB, LG] begin
    @test GO.coverage($empty_poly, cell_extremes...) == 0.0
    @test GO.coverage($degenerate_poly, cell_extremes...) == 0.0
end

# Test with multipolygons
mp1 = GI.MultiPolygon([p1, p2])
mp2 = GI.MultiPolygon([p3, p4])
mp3 = GI.MultiPolygon([p12, p13, p14])  # All zero area cases
mp4 = GI.MultiPolygon([p18, p19, p20])  # Mixed cases with holes

@testset_implementations "Multipolygons" [GI, GB, LG] begin
    @test GO.coverage($mp1, cell_extremes...) ≈ LG.area(LG.intersection($p1, cell_poly)) + LG.area(LG.intersection($p2, cell_poly))
    @test GO.coverage($mp2, cell_extremes...) ≈ LG.area(LG.intersection($p3, cell_poly)) + LG.area(LG.intersection($p4, cell_poly))
    @test GO.coverage($mp3, cell_extremes...) == 0.0
    @test GO.coverage($mp4, cell_extremes...) ≈ LG.area(LG.intersection($p18, cell_poly)) + LG.area(LG.intersection($p19, cell_poly)) + LG.area(LG.intersection($p20, cell_poly))
end
