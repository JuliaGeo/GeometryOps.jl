export get_geometries

_linearring(geom::GI.LineString) = GI.LinearRing(parent(geom); extent=geom.extent, crs=geom.crs)
_linearring(geom::GI.LinearRing) = geom


function get_geometries(x; geometrycolumn=nothing)
    # Check if already an AbstractArray to avoid unnecessary collection
    if x isa AbstractArray
        # Handle offset axes if needed - collect to ensure 1-based indexing
        if Base.has_offset_axes(x)
            return collect(x)
        else
            return x
        end
    elseif GI.trait(x) isa GI.GeometryCollectionTrait
        # Handle GeometryCollection properly with collect
        # Check this BEFORE isgeometry since GeometryCollections are also geometries
        # Recursively process in case elements need unwrapping
        geoms = collect(GI.getgeom(x))
        return get_geometries(geoms; geometrycolumn=geometrycolumn)
    elseif GI.isgeometry(x)
        return x
    elseif GI.isfeature(x)
        return GI.geometry(x)
    elseif GI.isfeaturecollection(x)
        return [GI.geometry(f) for f in GI.getfeature(x)]
    elseif Tables.istable(x)
        # Handle multiple geometry columns with kwarg
        cols = Tables.columns(x)
        colnames = Tables.columnnames(cols)

        geom_col = if geometrycolumn !== nothing
            geometrycolumn
        else
            geom_cols = GI.geometrycolumns(x)
            isempty(geom_cols) ? nothing : first(geom_cols)
        end

        if geom_col !== nothing && geom_col in colnames
            return collect(Tables.getcolumn(cols, geom_col))
        else
            error("No valid geometry column found in table")
        end
    else
        # the abstract array case has already been handled in the first
        # branch of the main if statement.  this branch is only for non-array iterables.
        collected = collect(x)
        if collected isa AbstractArray && GI.trait(first(collected)) isa GI.AbstractGeometryTrait
            return collected
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
