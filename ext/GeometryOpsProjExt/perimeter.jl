import GeometryOps: perimeter, applyreduce, TraitTarget, WithTrait
import GeoInterface as GI

function perimeter(m::Geodesic, geom, ::Type{T} = Float64; init = zero(T), kwargs...) where T
    # Create a Proj geodesic object using the ellipsoid parameters from the Geodesic manifold
    proj_geodesic = Ref(Proj.geod_geodesic(m.semimajor_axis, 1/m.inv_flattening))
    proj_polygon = Ref(Proj._null(Proj.geod_polygon))
    
    function _perimeter_geodesic_inner(trait, geom)
        @assert GI.npoint(geom) >= 2 "Geodesic perimeter requires at least 2 points"
        
        # Initialize the polygon
        proj_polygon[] = Proj._null(Proj.geod_polygon)
        Proj.geod_polygon_init(proj_polygon, 1)
        
        # Add all points to the polygon
        for point in GI.getpoint(trait, geom)
            lat, lon = GI.y(point), GI.x(point)  # Proj expects lat, lon order
            Proj.geod_polygon_addpoint(proj_geodesic, proj_polygon, lat, lon)
        end
        
        # Compute the polygon properties
        # geod_polygon_compute returns (num_vertices, perimeter, area)
        area_result, perimeter_result = Proj.geod_polygon_compute(proj_geodesic[], proj_polygon[], false, true)
        
        return T(perimeter_result)
    end
    
    return applyreduce(
        WithTrait(_perimeter_geodesic_inner), 
        +, 
        TraitTarget(GI.AbstractCurveTrait), 
        geom; init, kwargs...
    )
end