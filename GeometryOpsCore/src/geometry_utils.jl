
_linearring(geom::GI.LineString) = GI.LinearRing(parent(geom); extent=geom.extent, crs=geom.crs)
_linearring(geom::GI.LinearRing) = geom


function get_geometries(x; geometrycolumn=nothing)
    # Check if already an AbstractArray to avoid unnecessary collection
    if x isa AbstractArray
        # Handle offset axes if needed
        if Base.has_offset_axes(x)
            return vec(x)
        else
            return x
        end
    elseif GI.isgeometry(x)
        return x
    elseif GI.isgeometrycollection(x)
        # Handle GeometryCollection properly with collect
        return collect(GI.getgeom(x))
    elseif GI.isfeature(x)
        return GI.geometry(x)
    elseif GI.isfeaturecollection(x)
        return [GI.geometry(f) for f in GI.getfeature(x)]
    elseif Tables.istable(x)
        # Handle multiple geometry columns with kwarg
        geom_col = if geometrycolumn !== nothing
            geometrycolumn
        else
            geom_cols = GI.geometrycolumns(x)
            isempty(geom_cols) ? nothing : first(geom_cols)
        end
        
        if geom_col !== nothing && Tables.hascolumn(x, geom_col)
            return Tables.getcolumn(x, geom_col)
        else
            error("No valid geometry column found in table")
        end
    else
        c = collect(x)
        if c isa AbstractArray && GI.trait(first(c)) isa GI.AbstractGeometryTrait
            return c
        else
            throw(ArgumentError("""
                Expected a geometry, feature, feature collection, table with geometry column, or iterable of geometries.
                Got $(typeof(x)).
                
                The input must be one of:
                - A GeoInterface geometry (has trait <: AbstractGeometryTrait)
                - A GeoInterface feature (has trait FeatureTrait) 
                - A GeoInterface feature collection (has trait FeatureCollectionTrait)
                - A Tables.jl table with a geometry column
                - An iterable containing geometries
                """))
        end
    end
end
