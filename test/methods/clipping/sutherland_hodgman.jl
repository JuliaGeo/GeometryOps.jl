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
end
