module GeometryOpsTestHelpersArchGDALExt

using GeometryOpsTestHelpers
using GeoInterface
using ArchGDAL

function __init__()
    # Register ArchGDAL in the test modules list
    push!(GeometryOpsTestHelpers.TEST_MODULES, ArchGDAL)
end

# Monkey-patch ArchGDAL to handle polygon conversion correctly
@eval ArchGDAL begin
    function GeoInterface.convert(
        ::Type{T},
        type::GeoInterface.PolygonTrait,
        geom,
    ) where {T<:IGeometry}
        f = get(lookup_method, typeof(type), nothing)
        isnothing(f) && error(
            "Cannot convert an object of $(typeof(geom)) with the $(typeof(type)) trait (yet). Please report an issue.",
        )
        poly = createpolygon()
        foreach(GeoInterface.getring(geom)) do ring
            xs = GeoInterface.x.(GeoInterface.getpoint(ring)) |> collect
            ys = GeoInterface.y.(GeoInterface.getpoint(ring)) |> collect
            subgeom = unsafe_createlinearring(xs, ys)
            result = GDAL.ogr_g_addgeometrydirectly(poly, subgeom)
            @ogrerr result "Failed to add linearring."
        end
        return poly
    end
end

end # module GeometryOpsTestHelpersArchGDALExt
