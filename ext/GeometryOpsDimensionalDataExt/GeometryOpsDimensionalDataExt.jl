module GeometryOpsDimensionalDataExt

import DimensionalData as DD
import GeometryOps as GO
import GeoInterface as GI

# Polygonize an `AbstractDimArray` in the coordinate space of its `X`/`Y` lookups
# (via their interval bounds) rather than the raw integer axes.
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
