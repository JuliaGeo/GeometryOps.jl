module GeometryOpsDimensionalDataExt

import DimensionalData as DD
import GeometryOps as GO
import GeoInterface as GI

function GO.polygonize(A::DD.AbstractDimArray{T, N}; dims=(DD.XDim, DD.YDim), crs=GI.crs(A), kw...) where {T, N}
    # Extract the lookups of the dimensions we care about
    lookups = DD.lookup(A, dims)
    # Convert them to interval-bound vectors
    bounds_vecs = if DD.isintervals(lookups)
        map(DD.intervalbounds, lookups)
    else
        @warn "`polygonsize` is not possible for `Points` sampling, as polygons cover space by definition. Treating as `Intervals`, but this may not be appropriate"
        map(lookups) do l
            DD.intervalbounds(DD.set(l, DD.Intervals()))
        end
    end

    # This tree switches on the array type.
    # If ndims < 2, then we can't polygonize anyway, so we error.
    # If ndims > 2, then we return a DimArray across the remaining dimensions
    # not provided in `dims`.  Each slice in `dims` (usually X and Y) 
    if N < 2
        error("Array had $N dimensions, but it must have at least two!")
    elseif N == 2
        # If ndims == 2, then we polygonize directly, and return the output as requested.
        Ap = PermutedDimsArray(A, dims)
        GO.polygonize(bounds_vecs..., Ap; crs, kw...)
    else
        map(eachslice(A; dims = DD.otherdims(A, dims))) do a
            ap = PermutedDimsArray(a, dims)
            GO.polygonize(bounds_vecs..., a; crs, kw...)
        end
    end
            
end


end
