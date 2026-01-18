using Test
import GeometryOps as GO
import GeoInterface as GI

@testset "ConvexConvexSutherlandHodgman" begin
    @testset "Basic intersection" begin
        # Two overlapping squares - intersection is 1x1 square
        square1 = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(result) ≈ 1.0 atol=1e-10
    end

    @testset "No intersection" begin
        # Disjoint squares
        square1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0), (5.0, 5.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(result) ≈ 0.0 atol=1e-10
    end

    @testset "One contains other" begin
        # Large square contains small square
        large = GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]])
        small = GI.Polygon([[(2.0, 2.0), (4.0, 2.0), (4.0, 4.0), (2.0, 4.0), (2.0, 2.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), large, small)
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(result) ≈ 4.0 atol=1e-10

        # Reverse order should give same result
        result2 = GO.intersection(GO.ConvexConvexSutherlandHodgman(), small, large)
        @test GO.area(result2) ≈ 4.0 atol=1e-10
    end

    @testset "Triangles" begin
        # Two overlapping triangles (both CCW winding)
        tri1 = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (2.0, 4.0), (0.0, 0.0)]])
        tri2 = GI.Polygon([[(0.0, 2.0), (2.0, -2.0), (4.0, 2.0), (0.0, 2.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), tri1, tri2)
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(result) > 0
    end

    @testset "Identical polygons" begin
        # Same polygon should return itself
        square = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square, square)
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(result) ≈ 4.0 atol=1e-10
    end

    @testset "Shared edge" begin
        # Two squares sharing an edge
        square1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        square2 = GI.Polygon([[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 0.0)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), square1, square2)
        @test GI.trait(result) isa GI.PolygonTrait
        # Shared edge only - area should be 0 or near 0
        @test GO.area(result) ≈ 0.0 atol=1e-10
    end

    @testset "Unsupported geometry types" begin
        square = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        point = GI.Point(0.5, 0.5)

        @test_throws ArgumentError GO.intersection(GO.ConvexConvexSutherlandHodgman(), square, point)
        @test_throws ArgumentError GO.intersection(GO.ConvexConvexSutherlandHodgman(), point, square)
    end

    @testset "Pentagons" begin
        # Helper to create regular polygon (CCW winding)
        function regular_polygon(n, radius, center_x, center_y)
            coords = Tuple{Float64,Float64}[]
            for i in 0:n-1
                θ = 2π * i / n
                push!(coords, (center_x + radius * cos(θ), center_y + radius * sin(θ)))
            end
            push!(coords, coords[1])  # close the ring
            return GI.Polygon([coords])
        end

        # Two overlapping pentagons
        pent1 = regular_polygon(5, 2.0, 0.0, 0.0)
        pent2 = regular_polygon(5, 2.0, 1.5, 0.0)

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), pent1, pent2)
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(result) > 0
        # Intersection should be smaller than either pentagon
        @test GO.area(result) < GO.area(pent1)
    end

    @testset "Hexagons" begin
        function regular_polygon(n, radius, center_x, center_y)
            coords = Tuple{Float64,Float64}[]
            for i in 0:n-1
                θ = 2π * i / n
                push!(coords, (center_x + radius * cos(θ), center_y + radius * sin(θ)))
            end
            push!(coords, coords[1])
            return GI.Polygon([coords])
        end

        # Two overlapping hexagons
        hex1 = regular_polygon(6, 2.0, 0.0, 0.0)
        hex2 = regular_polygon(6, 2.0, 2.0, 0.0)

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), hex1, hex2)
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(result) > 0
        @test GO.area(result) < GO.area(hex1)

        # Hexagon containing smaller hexagon
        hex_large = regular_polygon(6, 3.0, 0.0, 0.0)
        hex_small = regular_polygon(6, 1.0, 0.0, 0.0)

        result2 = GO.intersection(GO.ConvexConvexSutherlandHodgman(), hex_large, hex_small)
        @test GO.area(result2) ≈ GO.area(hex_small) atol=1e-10
    end

    @testset "Octagons" begin
        function regular_polygon(n, radius, center_x, center_y)
            coords = Tuple{Float64,Float64}[]
            for i in 0:n-1
                θ = 2π * i / n
                push!(coords, (center_x + radius * cos(θ), center_y + radius * sin(θ)))
            end
            push!(coords, coords[1])
            return GI.Polygon([coords])
        end

        # Two overlapping octagons
        oct1 = regular_polygon(8, 2.0, 0.0, 0.0)
        oct2 = regular_polygon(8, 2.0, 1.0, 1.0)

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), oct1, oct2)
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(result) > 0
        @test GO.area(result) < GO.area(oct1)

        # Identical octagons should return same area
        result2 = GO.intersection(GO.ConvexConvexSutherlandHodgman(), oct1, oct1)
        @test GO.area(result2) ≈ GO.area(oct1) atol=1e-10
    end

    @testset "Bordering rectangles with offset" begin
        # Two rectangles sharing an edge but offset vertically
        # rect1: [0,2] x [0,2], rect2: [2,4] x [0.5,2.5]
        # They share the edge at x=2 but are offset by 0.5 in y
        # Intersection should be zero area (just a line segment)
        rect1 = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
        rect2 = GI.Polygon([[(2.0, 0.5), (4.0, 0.5), (4.0, 2.5), (2.0, 2.5), (2.0, 0.5)]])

        result = GO.intersection(GO.ConvexConvexSutherlandHodgman(), rect1, rect2)
        @test GI.trait(result) isa GI.PolygonTrait
        @test GO.area(result) ≈ 0.0 atol=1e-10

        # Two rectangles sharing an edge but offset horizontally
        # rect3: [0,2] x [0,2], rect4: [0.5,2.5] x [2,4]
        rect3 = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
        rect4 = GI.Polygon([[(0.5, 2.0), (2.5, 2.0), (2.5, 4.0), (0.5, 4.0), (0.5, 2.0)]])

        result2 = GO.intersection(GO.ConvexConvexSutherlandHodgman(), rect3, rect4)
        @test GI.trait(result2) isa GI.PolygonTrait
        @test GO.area(result2) ≈ 0.0 atol=1e-10

        # Rectangles that just touch at a corner
        rect5 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        rect6 = GI.Polygon([[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)]])

        result3 = GO.intersection(GO.ConvexConvexSutherlandHodgman(), rect5, rect6)
        @test GI.trait(result3) isa GI.PolygonTrait
        @test GO.area(result3) ≈ 0.0 atol=1e-10
    end

    @testset "Spherical helpers" begin
        using GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic

        @testset "_point_in_convex_spherical_polygon" begin
            transform = UnitSphereFromGeographic()

            # CCW square near equator
            square_pts = UnitSphericalPoint{Float64}[
                transform((0.0, 0.0)),
                transform((2.0, 0.0)),
                transform((2.0, 2.0)),
                transform((0.0, 2.0))
            ]

            inside_pt = transform((1.0, 1.0))
            outside_pt = transform((5.0, 5.0))
            edge_pt = transform((1.0, 0.0))  # On edge

            @test GO._point_in_convex_spherical_polygon(inside_pt, square_pts) == true
            @test GO._point_in_convex_spherical_polygon(outside_pt, square_pts) == false
            # Edge point should be considered inside (>= 0 check)
            @test GO._point_in_convex_spherical_polygon(edge_pt, square_pts) == true
        end
    end

    @testset "Spherical ConvexConvexSutherlandHodgman" begin
        using GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic

        function spherical_polygon(coords)
            transform = UnitSphereFromGeographic()
            points = [transform((lon, lat)) for (lon, lat) in coords]
            push!(points, points[1])
            return GI.Polygon([points])
        end

        spherical_area(poly) = GO.area(GO.Spherical(), poly)

        @testset "Basic intersection - small region" begin
            square1 = spherical_polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)])
            square2 = spherical_polygon([(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)])

            result = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                square1, square2
            )
            @test GI.trait(result) isa GI.PolygonTrait
            @test spherical_area(result) > 0
        end

        @testset "No intersection" begin
            square1 = spherical_polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
            square2 = spherical_polygon([(10.0, 10.0), (11.0, 10.0), (11.0, 11.0), (10.0, 11.0)])

            result = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                square1, square2
            )
            @test spherical_area(result) ≈ 0.0 atol=1e-10
        end

        @testset "Partial overlap" begin
            # Two overlapping squares - not containment, just partial overlap
            square1 = spherical_polygon([(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0)])
            square2 = spherical_polygon([(1.5, 1.5), (4.5, 1.5), (4.5, 4.5), (1.5, 4.5)])

            result = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                square1, square2
            )
            @test GI.trait(result) isa GI.PolygonTrait
            @test spherical_area(result) > 0
            # Intersection should be smaller than both inputs
            @test spherical_area(result) < spherical_area(square1)
            @test spherical_area(result) < spherical_area(square2)
        end

        @testset "Triangles" begin
            # Two overlapping triangles - like the planar test
            tri1 = spherical_polygon([(0.0, 0.0), (4.0, 0.0), (2.0, 4.0)])
            tri2 = spherical_polygon([(0.0, 2.0), (2.0, -2.0), (4.0, 2.0)])

            result = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                tri1, tri2
            )
            @test GI.trait(result) isa GI.PolygonTrait
            @test spherical_area(result) > 0
        end

        @testset "Near pole" begin
            tri1 = spherical_polygon([(0.0, 85.0), (120.0, 85.0), (240.0, 85.0)])
            tri2 = spherical_polygon([(60.0, 85.0), (180.0, 85.0), (300.0, 85.0)])

            result = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                tri1, tri2
            )
            @test GI.trait(result) isa GI.PolygonTrait
            @test spherical_area(result) > 0
        end

        @testset "Input validation" begin
            planar_poly = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])

            @test_throws ArgumentError GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                planar_poly, planar_poly
            )
        end

        @testset "One contains other" begin
            large = spherical_polygon([(-5.0, -5.0), (5.0, -5.0), (5.0, 5.0), (-5.0, 5.0)])
            small = spherical_polygon([(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)])

            # Both orderings should return the small polygon's area
            result1 = GO.intersection(GO.ConvexConvexSutherlandHodgman(GO.Spherical()), large, small)
            result2 = GO.intersection(GO.ConvexConvexSutherlandHodgman(GO.Spherical()), small, large)

            @test spherical_area(result1) ≈ spherical_area(small) rtol=1e-3
            @test spherical_area(result2) ≈ spherical_area(small) rtol=1e-3
        end

        @testset "Adjacent cells with shared edge" begin
            # Adjacent cells sharing an edge should have zero-area intersection
            # This tests that points exactly on the clip edge (orient=0) are handled
            # correctly without introducing numerical errors from intersection computation

            # Small adjacent cells (1° x 1°)
            cell1 = spherical_polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
            cell2 = spherical_polygon([(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)])

            result = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                cell1, cell2
            )
            @test spherical_area(result) == 0.0

            # Reverse order should also be zero
            result_rev = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                cell2, cell1
            )
            @test spherical_area(result_rev) == 0.0

            # Larger adjacent cells (10° x 10°)
            rect1 = spherical_polygon([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
            rect2 = spherical_polygon([(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)])

            result_large = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                rect1, rect2
            )
            @test spherical_area(result_large) == 0.0

            # Vertically adjacent cells
            top = spherical_polygon([(0.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)])
            bottom = spherical_polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])

            result_vert = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                top, bottom
            )
            @test spherical_area(result_vert) == 0.0
        end

        @testset "Cells touching at corner only" begin
            # Cells that share only a corner point should have zero intersection
            cell1 = spherical_polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
            cell2 = spherical_polygon([(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)])

            result = GO.intersection(
                GO.ConvexConvexSutherlandHodgman(GO.Spherical()),
                cell1, cell2
            )
            @test spherical_area(result) == 0.0
        end
    end
end

# Spherical Sutherland-Hodgman Tests
@testset "ConvexConvexSutherlandHodgman - Spherical" begin
    using GeometryOps.UnitSpherical: UnitSphereFromGeographic

    # Transform lon/lat to UnitSphericalPoint
    _transform = UnitSphereFromGeographic()
    lonlat_to_point(lon, lat) = _transform((lon, lat))

    alg = GO.ConvexConvexSutherlandHodgman(GO.Spherical())

    @testset "Basic spherical intersection" begin
        # Two overlapping 2°×2° spherical cells
        poly_a = GI.Polygon([[
            lonlat_to_point(0.0, 0.0),
            lonlat_to_point(2.0, 0.0),
            lonlat_to_point(2.0, 2.0),
            lonlat_to_point(0.0, 2.0),
            lonlat_to_point(0.0, 0.0),
        ]])

        poly_b = GI.Polygon([[
            lonlat_to_point(1.0, 1.0),
            lonlat_to_point(3.0, 1.0),
            lonlat_to_point(3.0, 3.0),
            lonlat_to_point(1.0, 3.0),
            lonlat_to_point(1.0, 1.0),
        ]])

        result = GO.intersection(alg, poly_a, poly_b; target=GI.PolygonTrait())
        area = GO.area(GO.Spherical(), result)

        # Should be approximately 1°×1° intersection area
        @test area > 0
        @test area < GO.area(GO.Spherical(), poly_a)
    end

    @testset "Adjacent polygons (shared edge) - THE FIX" begin
        # Two 1°×1° cells sharing lon=126 edge
        # This is the main bug that was reported
        poly_a = GI.Polygon([[
            lonlat_to_point(125.0, 53.0),
            lonlat_to_point(126.0, 53.0),
            lonlat_to_point(126.0, 54.0),
            lonlat_to_point(125.0, 54.0),
            lonlat_to_point(125.0, 53.0),
        ]])

        poly_b = GI.Polygon([[
            lonlat_to_point(126.0, 53.0),
            lonlat_to_point(127.0, 53.0),
            lonlat_to_point(127.0, 54.0),
            lonlat_to_point(126.0, 54.0),
            lonlat_to_point(126.0, 53.0),
        ]])

        result_ab = GO.intersection(alg, poly_a, poly_b; target=GI.PolygonTrait())
        result_ba = GO.intersection(alg, poly_b, poly_a; target=GI.PolygonTrait())

        area_ab = GO.area(GO.Spherical(), result_ab)
        area_ba = GO.area(GO.Spherical(), result_ba)

        # Adjacent polygons should have zero/negligible intersection
        @test area_ab < 1e-10
        @test area_ba < 1e-10

        # Operation should be symmetric
        @test area_ab ≈ area_ba atol=1e-10
    end

    @testset "Vertex on edge (no overlap)" begin
        # Subject polygon (lon 125-126, lat 53-54)
        poly_a = GI.Polygon([[
            lonlat_to_point(125.0, 53.0),
            lonlat_to_point(126.0, 53.0),
            lonlat_to_point(126.0, 54.0),
            lonlat_to_point(125.0, 54.0),
            lonlat_to_point(125.0, 53.0),
        ]])

        # Polygon with vertex at (126.0, 53.5) - ON poly_a's lon=126 edge but outside
        poly_b = GI.Polygon([[
            lonlat_to_point(126.0, 53.5),  # ON poly_a's lon=126 edge
            lonlat_to_point(127.0, 53.0),
            lonlat_to_point(127.0, 54.0),
            lonlat_to_point(126.0, 53.5),
        ]])

        result = GO.intersection(alg, poly_a, poly_b; target=GI.PolygonTrait())
        area = GO.area(GO.Spherical(), result)

        # Should be zero, not the area of poly_b!
        @test area < 1e-10
    end

    @testset "Overlapping spherical polygons" begin
        poly_a = GI.Polygon([[
            lonlat_to_point(125.0, 53.0),
            lonlat_to_point(127.0, 53.0),
            lonlat_to_point(127.0, 55.0),
            lonlat_to_point(125.0, 55.0),
            lonlat_to_point(125.0, 53.0),
        ]])

        poly_b = GI.Polygon([[
            lonlat_to_point(126.0, 54.0),
            lonlat_to_point(128.0, 54.0),
            lonlat_to_point(128.0, 56.0),
            lonlat_to_point(126.0, 56.0),
            lonlat_to_point(126.0, 54.0),
        ]])

        result = GO.intersection(alg, poly_a, poly_b; target=GI.PolygonTrait())
        area = GO.area(GO.Spherical(), result)

        # Should be approximately 1°×1° = area of overlap region
        expected_poly = GI.Polygon([[
            lonlat_to_point(126.0, 54.0),
            lonlat_to_point(127.0, 54.0),
            lonlat_to_point(127.0, 55.0),
            lonlat_to_point(126.0, 55.0),
            lonlat_to_point(126.0, 54.0),
        ]])
        expected_area = GO.area(GO.Spherical(), expected_poly)

        @test area ≈ expected_area rtol=0.05
    end

    @testset "Containment" begin
        outer = GI.Polygon([[
            lonlat_to_point(120.0, 50.0),
            lonlat_to_point(130.0, 50.0),
            lonlat_to_point(130.0, 60.0),
            lonlat_to_point(120.0, 60.0),
            lonlat_to_point(120.0, 50.0),
        ]])

        inner = GI.Polygon([[
            lonlat_to_point(124.0, 54.0),
            lonlat_to_point(126.0, 54.0),
            lonlat_to_point(126.0, 56.0),
            lonlat_to_point(124.0, 56.0),
            lonlat_to_point(124.0, 54.0),
        ]])

        result = GO.intersection(alg, outer, inner; target=GI.PolygonTrait())
        inner_area = GO.area(GO.Spherical(), inner)
        result_area = GO.area(GO.Spherical(), result)

        @test result_area ≈ inner_area rtol=0.01
    end

    @testset "Original bug report case" begin
        # Subject polygon (lon 125-126, lat 53-54)
        subject = GI.Polygon([[
            lonlat_to_point(125.0, 53.0),
            lonlat_to_point(126.0, 53.0),
            lonlat_to_point(126.0, 54.0),
            lonlat_to_point(125.0, 54.0),
            lonlat_to_point(125.0, 53.0),
        ]])

        # Clip polygon - adjacent, with vertex at (126.0, 53.23) ON subject's edge
        clip = GI.Polygon([[
            lonlat_to_point(126.0, 53.23),   # ON subject's lon=126 edge!
            lonlat_to_point(126.86, 52.32),
            lonlat_to_point(127.86, 53.25),
            lonlat_to_point(126.95, 54.15),
            lonlat_to_point(126.0, 53.23),
        ]])

        result = GO.intersection(alg, subject, clip; target=GI.PolygonTrait())
        result_area = GO.area(GO.Spherical(), result)
        clip_area = GO.area(GO.Spherical(), clip)

        # BUG FIX: result_area should be ~0, NOT clip_area
        # The ratio should be ~0, not ~1
        @test result_area < clip_area * 0.01  # Less than 1% of clip area
    end
end
