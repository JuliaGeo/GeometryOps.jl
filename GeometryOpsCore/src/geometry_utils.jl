
_linearring(geom::GI.LineString) = GI.LinearRing(parent(geom); extent=geom.extent, crs=geom.crs)
_linearring(geom::GI.LinearRing) = geom


function get_geometries(x)
    if GI.isgeometry(x)
        return x
    # elseif GI.isgeometrycollection(x)
    #     return GI.getgeom(x)
    elseif GI.isfeature(x)
        return GI.geometry(x)
    elseif GI.isfeaturecollection(x)
        return [GI.geometry(f) for f in GI.getfeature(x)]
    elseif Tables.istable(x) && Tables.hascolumn(x, first(GI.geometrycolumns(x)))
        return Tables.getcolumn(x, first(GI.geometrycolumns(x)))
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
