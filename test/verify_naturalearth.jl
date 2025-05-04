using GeometryOps
using GeoInterface
using NaturalEarth

# Load NaturalEarth admin_0_countries at 10m resolution
countries = naturalearth("admin_0_countries", 10)

# Test function with different central meridians
function test_all_central_meridians(countries)
    for central_meridian in -180:30:180
        left_edge = central_meridian - 180
        right_edge = central_meridian + 180
        
        println("Testing central meridian: $central_meridian")
        
        # Process each country
        for i in 1:length(countries)
            country = countries[i]
            
            # Apply the cut_at_antimeridian function
            result = cut_at_antimeridian(
                country, 
                left_edge=Float64(left_edge), 
                center_edge=Float64(central_meridian), 
                right_edge=Float64(right_edge)
            )
            
            # Verify that the result is valid
            # Here we just check that we don't get errors
            if GI.geomtrait(result) isa GI.MultiPolygonTrait
                @assert length(GI.getgeom(result)) > 0
            end
        end
        
        println("Successfully processed all countries with central meridian: $central_meridian")
    end
end

test_all_central_meridians(countries)
