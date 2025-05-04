using Test
using GeometryOps
using GeoInterface
using JSON

# Path to test data
const TEST_DATA_DIR = joinpath(@__DIR__, "data")

# Helper to read JSON test files
function read_test_geojson(filename)
    path = joinpath(TEST_DATA_DIR, filename)
    data = JSON.parsefile(path)
    return GI.read(data)
end

@testset "cut_at_antimeridian" begin
    @testset "Simple cases" begin
        # Test with a point (shouldn't be affected)
        point = GI.Point(0.0, 0.0)
        @test cut_at_antimeridian(point) === point
        
        # Test with a polygon that doesn't cross the antimeridian
        poly = GI.Polygon([[(10.0, 10.0), (20.0, 10.0), (20.0, 20.0), (10.0, 20.0), (10.0, 10.0)]])
        @test cut_at_antimeridian(poly) === poly
    end
    
    @testset "Crossing the antimeridian" begin
        # Test with a polygon that crosses the antimeridian
        poly = GI.Polygon([[(170.0, 40.0), (170.0, 50.0), (-170.0, 50.0), (-170.0, 40.0), (170.0, 40.0)]])
        result = cut_at_antimeridian(poly)
        
        # Should result in a MultiPolygon with two parts
        @test GI.geomtrait(result) isa GI.MultiPolygonTrait
        @test length(GI.getgeom(result)) == 2
        
        # Check that each part is on one side of the antimeridian
        part1, part2 = GI.getgeom(result)
        part1_coords = collect(GI.getpoint(GI.getexterior(part1)))
        part2_coords = collect(GI.getpoint(GI.getexterior(part2)))
        
        # All coordinates in part1 should be on one side
        @test all(p -> p[1] >= 0, part1_coords) || all(p -> p[1] <= 0, part1_coords)
        # All coordinates in part2 should be on the other side
        @test all(p -> p[1] >= 0, part2_coords) || all(p -> p[1] <= 0, part2_coords)
        # The two parts should be on opposite sides
        @test (all(p -> p[1] >= 0, part1_coords) && all(p -> p[1] <= 0, part2_coords)) ||
              (all(p -> p[1] <= 0, part1_coords) && all(p -> p[1] >= 0, part2_coords))
    end
    
    @testset "Great circle vs flat interpolation" begin
        # Test the difference between great circle and flat interpolation
        poly = GI.Polygon([[(170.0, 10.0), (170.0, 50.0), (-170.0, 50.0), (-170.0, 10.0), (170.0, 10.0)]])
        
        result_great_circle = cut_at_antimeridian(poly, great_circle=true)
        result_flat = cut_at_antimeridian(poly, great_circle=false)
        
        # They should produce different results for polygons with large latitude spans
        # Extract coordinates for comparison
        gc_part1 = collect(GI.getpoint(GI.getexterior(GI.getgeom(result_great_circle)[1])))
        flat_part1 = collect(GI.getpoint(GI.getexterior(GI.getgeom(result_flat)[1])))
        
        # Check that the interpolation points are different
        # We'll compare the points that should be on the antimeridian
        # These will be at different latitudes for great circle vs flat
        gc_antimeridian_points = filter(p -> abs(abs(p[1]) - 180.0) < 1e-6, gc_part1)
        flat_antimeridian_points = filter(p -> abs(abs(p[1]) - 180.0) < 1e-6, flat_part1)
        
        # There should be at least one point on the antimeridian
        @test !isempty(gc_antimeridian_points)
        @test !isempty(flat_antimeridian_points)
        
        # The latitudes should be different
        gc_lats = sort([p[2] for p in gc_antimeridian_points])
        flat_lats = sort([p[2] for p in flat_antimeridian_points])
        
        # At least one latitude should differ by more than a small epsilon
        @test any(abs(gc_lats[i] - flat_lats[i]) > 1e-6 for i in 1:min(length(gc_lats), length(flat_lats)))
    end
    
    @testset "Different central meridians" begin
        # Test with a different central meridian
        poly = GI.Polygon([[(-10.0, 40.0), (-10.0, 50.0), (10.0, 50.0), (10.0, 40.0), (-10.0, 40.0)]])
        
        # Cut at 0° meridian
        result = cut_at_antimeridian(poly, left_edge=-90.0, center_edge=0.0, right_edge=90.0)
        
        # Should result in a MultiPolygon with two parts
        @test GI.geomtrait(result) isa GI.MultiPolygonTrait
        @test length(GI.getgeom(result)) == 2
        
        # Check that each part is on one side of the 0° meridian
        part1, part2 = GI.getgeom(result)
        part1_coords = collect(GI.getpoint(GI.getexterior(part1)))
        part2_coords = collect(GI.getpoint(GI.getexterior(part2)))
        
        # All coordinates in part1 should be on one side of 0°
        @test all(p -> p[1] >= 0, part1_coords) || all(p -> p[1] <= 0, part1_coords)
        # All coordinates in part2 should be on the other side of 0°
        @test all(p -> p[1] >= 0, part2_coords) || all(p -> p[1] <= 0, part2_coords)
        # The two parts should be on opposite sides
        @test (all(p -> p[1] >= 0, part1_coords) && all(p -> p[1] <= 0, part2_coords)) ||
              (all(p -> p[1] <= 0, part1_coords) && all(p -> p[1] >= 0, part2_coords))
    end
    
    @testset "Ported Python tests" begin
        test_cases = [
            "simple",
            "split",
            "complex-split",
            "north-pole",
            "south-pole"
        ]
        
        for case in test_cases
            @testset "$case" begin
                input = read_test_geojson(joinpath("input", "$case.json"))
                
                # Test flat interpolation
                result_flat = cut_at_antimeridian(input, great_circle=false)
                expected_flat = read_test_geojson(joinpath("flat", "$case.json"))
                
                # Test spherical interpolation
                result_spherical = cut_at_antimeridian(input, great_circle=true)
                expected_spherical = read_test_geojson(joinpath("spherical", "$case.json"))
                
                # Check that the result has the same geometry type as expected
                @test GI.geomtrait(result_flat) == GI.geomtrait(expected_flat)
                @test GI.geomtrait(result_spherical) == GI.geomtrait(expected_spherical)
                
                # For MultiPolygons, check that they have the same number of parts
                if GI.geomtrait(expected_flat) isa GI.MultiPolygonTrait
                    @test length(GI.getgeom(result_flat)) == length(GI.getgeom(expected_flat))
                end
                
                if GI.geomtrait(expected_spherical) isa GI.MultiPolygonTrait
                    @test length(GI.getgeom(result_spherical)) == length(GI.getgeom(expected_spherical))
                end
            end
        end
    end
    
    @testset "Natural Earth data" begin
        # This would test with NaturalEarth data for different central meridians
        # For now, we'll skip this test as it requires the NaturalEarth package
        @test_skip "Test with NaturalEarth data for different central meridians"
    end
end
