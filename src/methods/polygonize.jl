# # Polygonizing raster data

export polygonize

# The methods in this file are able to convert a raster image into a set of polygons,
# by contour detection using a clockwise Moore neighborhood method.

## Example

#=
Here's a basic example, using the `Makie.peaks()` function.  First, let's investigate the nature of the function:

```@example polygonize
using Makie, GeometryOps
n = 49
xs, ys = LinRange(-3, 3, n), LinRange(-3, 3, n)
zs = Makie.peaks(n)
z_max_value = maximum(abs.(extrema(zs)))
f, a, p = heatmap(
    xs, ys, zs; 
    axis = (; aspect = DataAspect(), title = "Exact function")
)
cb = Colorbar(f[1, 2], p; label = "Z-value")
f 
```

Now, we can use the `polygonize` function to convert the raster data into polygons.

For this particular example, we chose a range of z-values between 0.8 and 3.2, 
which would provide two distinct polyogns with holes.

```@example polygonize
polygons = polygonize(xs, ys, 0.8 .< zs .< 3.2)
```
This returns a list of `GeometryBasics.Polygon`, which can be plotted immediately, 
or wrapped directly in a `GeometryBasics.MultiPolygon`.  Let's see how these look:

```@example polygonize
f, a, p = poly(polygons; label = "Polygonized polygons", axis = (; aspect = DataAspect()))
```

Finally, let's plot the Makie contour lines on top, to see how well the polygonization worked:
```@example polygonize
contour!(a, xs, ys, zs; labels = true, levels = [0.8, 3.2], label = "Contour lines")
f
```

## Implementation

The implementation follows:
=# 

"""
    polygonize(A::AbstractMatrix{Bool}; kw...)
    polygonize(f, A::AbstractMatrix; kw...)
    polygonize(xs, ys, A::AbstractMatrix{Bool}; kw...)
    polygonize(f, xs, ys, A::AbstractMatrix; kw...)

Polygonize an `AbstractMatrix` of values, currently to a single class of polygons.

For `AbstractArray{Bool}` function `f` is not needed. 

For other matrix eltypes, function `f` should return `true` or `false` 
based on the matrix values, translating to inside or outside the polygons.
These will return a single `MultiPolygon` of the `true` values. 

For `AbtractArray{<:Integer}` multiple `multipolygon`s are calculated
for each value in the array (or passed in `values` keyword), and returned
as a `FeatureCollection`.

If `xs` and `ys` are ranges, they are used as the pixel center points.
If they are `Vector` of `Tuple` they are used as the lower and upper bounds of each pixel.

# Keywords

- `minpoints`: ignore polygons with less than `minpoints` points.
- `values`: the values to turn into polygons for `Integer` arrays. 
    By default these are `union(A)`

# Example

```julia
using GeometryOps
multipolygon = polygonize(>(0.6), rand(100), minpoints=3)

using GeometryOps
featurecollection = polygonize(rand(1:4, 100) * (fill(1, 100))')

```
"""
polygonize(A::AbstractMatrix{Bool}; kw...) = polygonize(identity, A; kw...)
polygonize(f::Base.Callable, A::AbstractMatrix; kw...) = polygonize(f, axes(A)..., A; kw...)
polygonize(A::AbstractMatrix; kw...) = polygonize(axes(A)..., A; kw...)
polygonize(xs::AbstractVector, ys::AbstractVector, A::AbstractMatrix{Bool}; kw...) =
    polygonize(identity, xs, ys, A)
function polygonize(xs::AbstractVector, ys::AbstractVector, A::AbstractMatrix{<:Integer}; 
    values=Base.union(A),
    kw...
)
    # Create one feature per value
    features = map(values) do value
        multipolygon = polygonize(==(value), xs, ys, A)
        GI.Feature(multipolygon; properties=(; value))
    end 

    return GI.FeatureCollection(features)
end
function polygonize(f::Base.Callable, xs::AbstractRange, ys::AbstractRange, A::AbstractMatrix; 
    kw...
)
    # Make vectors of pixel bounds
    xhalf = step(xs) / 2
    yhalf = step(ys) / 2
    # Make bounds ranges first to avoid floating point error making gaps or overlaps
    xbounds = first(xs) - xhalf : step(xs) : last(xs) + xhalf
    ybounds = first(ys) - yhalf : step(ys) : last(ys) + yhalf
    Tx = eltype(xbounds)
    Ty = eltype(ybounds)
    xvec = Vector{Tuple{Tx,Tx}}(undef, length(xs))
    yvec = Vector{Tuple{Ty,Ty}}(undef, length(ys))
    for i in eachindex(xvec)
        xvec[i] = xbounds[i], xbounds[i+1]
    end
    for i in eachindex(yvec)
        yvec[i] = ybounds[i], ybounds[i+1]
    end
    return polygonize(f, xvec, yvec, A; kw...)
end

function polygonize(f, xs::AbstractVector{T}, ys::AbstractVector{T}, A::AbstractMatrix; 
    minpoints=0,
) where T
    # Define buffers for edges and rings
    edges = Dict{T,Vector{T}}()
    rings = Vector{T}[]

    @assert (length(xs), length(ys)) == size(A)

    # First get all the valid edges between filled and empty pixels
    si, sj = map(last, axes(A))
    @inbounds for i in axes(A, 1), j in axes(A, 2)
        if f(A[i, j]) # This is a pixel inside a polygon
            # xs and ys hold pixel bounds
            x1, x2 = xs[i]
            y1, y2 = ys[j]
            # We check the Von Neumann neighborhood to
            # decide what edges are needed, if any.

            (j == 1 || !f(A[i, j-1])) && push!(get!(() -> T[], edges, (x1, y1)), (x2, y1))
            (i == 1 || !f(A[i-1, j])) && push!(get!(() -> T[], edges, (x1, y2)), (x1, y1))
            (j == sj || !f(A[i, j+1])) && push!(get!(() -> T[], edges, (x2, y2)), (x1, y2))
            (i == si || !f(A[i+1, j])) && push!(get!(() -> T[], edges, (x2, y1)), (x2, y2))
        end
    end

    # Then we join the edges into rings
    # Sorting now lets us use `searchsortedlast` for speed later
    edgekeys = collect(keys(edges))
    while length(edges) > 0
        found = false
        local firstpoint
        while found == false && length(edgekeys) > 0
            # Take the first edge from the array
            firstpoint::T = last(edgekeys)
            if haskey(edges, firstpoint) && !isempty(edges[firstpoint])
                found = true
            else
                pop!(edgekeys)
            end
        end
        found == false && break
        nextpoints = edges[firstpoint]
        nextpoint = pop!(nextpoints)
        if length(nextpoints) == 1
            pop!(edgekeys)
        end
        edge = (firstpoint, nextpoint)
        ring = T[firstpoint, nextpoint]
        push!(rings, ring)
        while length(edges) > 0
            # Find an edge that matches the next point
            candidates = edges[nextpoint]
            # When there are two possible edges, 
            # choose the edge that has turned a corner
            if length(candidates) > 1 
                if (nextpoint[1] == candidates[2][1] || nextpoint[2] == candidates[2][2]) 
                    newedge = (nextpoint, pop!(candidates))
                else
                    newedge = (nexpoint, popfirst!(candidates))
                end
            else
                newedge = (nextpoint, candidates[1])
                delete!(edges, nextpoint)
            end
            edge = newedge
            # TODO: Here we actually need to check which edge maintains
            # the winding direction
            nextpoint = last(edge)
            # Close the ring if we get to the start
            if nextpoint == firstpoint
                push!(ring, nextpoint)
                break
            # Split off another ring if we touch the current ring anywhere else
            else
                i = findfirst(==(nextpoint), ring)
                if isnothing(i)
                    # Keep adding points to the ring
                    push!(ring, nextpoint)
                else
                    splitring = ring[i:lastindex(ring)]
                    deleteat!(ring, i:lastindex(ring))
                    push!(ring, nextpoint)
                    push!(splitring, nextpoint)
                    push!(rings, splitring)
                    continue
                end
            end
        end
    end

    # Define wrapped LinearRings, with embedded extents
    # so we only calculate them once
    linearrings = map(rings) do ring
        extent = GI.extent(GI.LinearRing(ring))
        GI.LinearRing(ring; extent)
    end

    # Separate exteriors from holes by winding direction
    direction = last(last(xs)) - first(first(xs)) * last(last(ys)) - first(first(ys))
    exterior_inds = if direction > 0 
        .!isclockwise.(linearrings)
    else
        isclockwise.(linearrings)
    end
    holes = linearrings[.!exterior_inds]
    polygons = map(view(linearrings, exterior_inds)) do lr
        GI.Polygon([lr]; extent=GI.extent(lr))
    end

    # Then we add the holes to the polygons they are inside of
    unused = fill(true, length(holes))
    for poly in polygons
        exterior = GI.Polygon(StaticArrays.SVector(GI.getexterior(poly)))
        for i in eachindex(holes)
            unused[i] || continue
            hole = holes[i]
            if covers(poly, hole)
                # Hole is in the exterior, so add it to the ring list
                push!(poly.geom, hole)
                # remove i
                unused[i] = false
                break
            end
        end
    end

    if isempty(polygons)
        # TODO: this really should return an emtpty MultiPolygon but
        # GeoInterface wrappers cant do that yet, which is not ideal...
        @warn "No polgons found, check your data or try another function for `f`"
        return nothing
    else
        # Otherwise return a wrapped MultiPolygon
        return GI.MultiPolygon(polygons)
    end
end

