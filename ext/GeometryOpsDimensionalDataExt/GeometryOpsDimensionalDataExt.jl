module GeometryOpsDimensionalDataExt

import DimensionalData as DD
import GeometryOps as GO
import GeoInterface as GI

# Polygonize a `DimArray` (or any `AbstractDimArray`, e.g. a `Raster`) using its
# `X`/`Y` lookup values rather than the raw integer axes.  We build interval
# bounds from the lookups so that the resulting polygons live in coordinate space.
function GO.polygonize(A::DD.AbstractDimArray; dims=(DD.X(), DD.Y()), crs=GI.crs(A), kw...)
    lookups = DD.lookup(A, dims)
    bounds_vecs = if all(DD.isintervals, lookups)
        map(DD.intervalbounds, lookups)
    else
        @warn "`polygonize` is not strictly correct for `Points` sampling, as polygons cover space by definition. Treating as `Intervals`, but this may not be appropriate."
        map(lookups) do l
            DD.intervalbounds(DD.set(l, DD.Intervals()))
        end
    end
    return GO.polygonize(bounds_vecs..., A; crs, kw...)
end

end
