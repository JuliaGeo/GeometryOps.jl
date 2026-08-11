import GeometryOps: perimeter, area, applyreduce, TraitTarget, WithTrait
import GeoInterface as GI
import GeoFormatTypes

# PROJ expects flattening, but GeometryOps stores inverse flattening; a sphere encodes both as zero.
_flattening(inv_flattening) = iszero(inv_flattening) ? zero(inv_flattening) : inv(inv_flattening)

function GeometryOps._area_auto(
    ::GI.AbstractGeographicTrait,
    crs::Union{GeoFormatTypes.GeoFormat,Proj.CRS},
    geom,
    ::Type{T};
    threaded=false,
    kwargs...,
) where T
    # A geographic CRS provides its ellipsoid and angular units, unlike the degree-only public Geodesic API.
    proj_crs = convert(Proj.CRS, crs)
    geodesic = _geodesic_manifold(proj_crs)
    longitude_scale, latitude_scale = _geographic_degree_scales(proj_crs)
    _geodesic_area(geodesic, geom, T, longitude_scale, latitude_scale; threaded, kwargs...)
end

function GeometryOps._area_auto(
    ::GI.UnknownTrait,
    crs::Union{GeoFormatTypes.GeoFormat,Proj.CRS},
    geom,
    ::Type{T};
    threaded=false,
    kwargs...,
) where T
    # UnknownTrait is a projected-trait subtype, so classify its attached CRS before using that fallback.
    proj_crs = convert(Proj.CRS, crs)
    if Proj.is_geographic(proj_crs)
        return GeometryOps._area_auto(GI.GeographicTrait(), proj_crs, geom, T; threaded, kwargs...)
    elseif Proj.is_projected(proj_crs)
        return GeometryOps._area_auto(GI.ProjectedTrait(), proj_crs, geom, T; threaded, kwargs...)
    end
    throw(ArgumentError("CRS $(crs) is neither geographic nor projected"))
end

function _geodesic_manifold(crs::Proj.CRS)
    ellipsoid = Proj.proj_get_ellipsoid(crs)
    ellipsoid == C_NULL && throw(ArgumentError("CRS has no ellipsoid"))
    try
        semimajor_axis = Ref{Cdouble}()
        semiminor_axis = Ref{Cdouble}()
        is_semiminor_computed = Ref{Cint}()
        inv_flattening = Ref{Cdouble}()
        Proj.proj_ellipsoid_get_parameters(
            ellipsoid,
            semimajor_axis,
            semiminor_axis,
            is_semiminor_computed,
            inv_flattening,
        ) == 0 && throw(ArgumentError("unable to get ellipsoid parameters from CRS"))
        return GeometryOps.Geodesic(; semimajor_axis=semimajor_axis[], inv_flattening=inv_flattening[])
    finally
        Proj.proj_destroy(ellipsoid)
    end
end

function _geographic_degree_scales(crs::Proj.CRS)
    coordinate_system = Proj.proj_crs_get_coordinate_system(crs)
    coordinate_system == C_NULL && throw(ArgumentError("CRS has no coordinate system"))
    try
        longitude_scale = nothing
        latitude_scale = nothing
        for axis in 0:(Proj.proj_cs_get_axis_count(coordinate_system) - 1)
            direction = Ref{Cstring}()
            unit_conversion = Ref{Cdouble}()
            Proj.proj_cs_get_axis_info(
                coordinate_system,
                axis,
                C_NULL,
                C_NULL,
                direction,
                unit_conversion,
                C_NULL,
                C_NULL,
                C_NULL,
            ) == 0 && throw(ArgumentError("unable to get coordinate-system axis $(axis + 1) from CRS"))
            direction[] == C_NULL && throw(ArgumentError("coordinate-system axis $(axis + 1) has no direction"))
            degrees_per_unit = unit_conversion[] * (180 / π)
            isfinite(degrees_per_unit) && degrees_per_unit > 0 ||
                throw(ArgumentError("coordinate-system axis $(axis + 1) has invalid angular unit conversion factor $(unit_conversion[])"))
            axis_direction = unsafe_string(direction[])
            # Axis direction identifies whether this angular scale belongs to longitude or latitude.
            if axis_direction == "east"
                isnothing(longitude_scale) || throw(ArgumentError("CRS has multiple longitude axes"))
                longitude_scale = degrees_per_unit
            elseif axis_direction == "west"
                isnothing(longitude_scale) || throw(ArgumentError("CRS has multiple longitude axes"))
                longitude_scale = -degrees_per_unit
            elseif axis_direction == "north"
                isnothing(latitude_scale) || throw(ArgumentError("CRS has multiple latitude axes"))
                latitude_scale = degrees_per_unit
            elseif axis_direction == "south"
                isnothing(latitude_scale) || throw(ArgumentError("CRS has multiple latitude axes"))
                latitude_scale = -degrees_per_unit
            end
        end
        isnothing(longitude_scale) && throw(ArgumentError("CRS has no longitude axis"))
        isnothing(latitude_scale) && throw(ArgumentError("CRS has no latitude axis"))
        return longitude_scale, latitude_scale
    finally
        Proj.proj_destroy(coordinate_system)
    end
end

function perimeter(m::Geodesic, geom, ::Type{T} = Float64; init = zero(T), kwargs...) where T
    # Create a Proj geodesic object using the ellipsoid parameters from the Geodesic manifold
    proj_geodesic = Ref(Proj.geod_geodesic(m.semimajor_axis, _flattening(m.inv_flattening)))
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

function _geodesic_area(
    m::Geodesic,
    geom,
    ::Type{T},
    longitude_scale,
    latitude_scale;
    threaded=false,
    init=zero(T),
    kwargs...,
) where T
    function _ring_area(ring)
        proj_geodesic = Ref(Proj.geod_geodesic(m.semimajor_axis, _flattening(m.inv_flattening)))
        proj_polygon = Ref(Proj._null(Proj.geod_polygon))
        Proj.geod_polygon_init(proj_polygon, 0)

        for point in GI.getpoint(ring)
            lat, lon = latitude_scale * GI.y(point), longitude_scale * GI.x(point)
            Proj.geod_polygon_addpoint(proj_geodesic, proj_polygon, lat, lon)
        end

        signed_area, _ = Proj.geod_polygon_compute(proj_geodesic[], proj_polygon[], false, true)
        return T(signed_area)
    end

    function _area_geodesic_inner(::GI.PolygonTrait, poly)
        GI.isempty(poly) && return zero(T)
        exterior_area = abs(_ring_area(GI.getexterior(poly)))
        hole_area = sum(hole -> abs(_ring_area(hole)), GI.gethole(poly); init = zero(T))
        return exterior_area - hole_area
    end
    _area_geodesic_inner(::GI.AbstractGeometryTrait, geom) = zero(T)

    return applyreduce(
        WithTrait(_area_geodesic_inner),
        +,
        GeometryOps._AREA_TARGETS,
        geom; threaded, init, kwargs...
    )
end

GeometryOps.area(m::Geodesic, geom, ::Type{T} = Float64; kwargs...) where T =
    _geodesic_area(m, geom, T, 1, 1; kwargs...)