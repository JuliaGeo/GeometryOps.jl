# # Polygonizing raster data

export polygonize

#=
The methods in this file are able to convert a raster image into a set of polygons,
by contour detection using a clockwise Moore neighborhood method.

The resulting polygons are snapped to the boundaries of the cells of the input raster,
so they will look different from traditional contours from a plotting package.

## Example

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
This returns a `GI.MultiPolygon`, which is directly plottable.  Let's see how these look:

```@example polygonize
f, a, p = poly(polygons; label = "Polygonized polygons", axis = (; aspect = DataAspect()))
```

Finally, let's plot the Makie contour lines on top, to see how the polygonization compares:
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

For `AbtractArray{<:Integer}`, a `multipolygon` is calculated
for each value in the array (or passed in `values` keyword), 
and these multipolygons and their associated values are returned
as a `FeatureCollection`.  You can convert this into a DataFrame
by calling `DataFrame(polygonize(...))`.

If `xs` and `ys` are ranges, they are used as the pixel/cell center points.
If they are `Vector` of `Tuple` they are used as the lower and upper bounds of each pixel/cell.

# Keywords

- `minpoints`: ignore polygons with less than `minpoints` points.
- `values`: the values to turn into polygons for `Integer` arrays. 
    By default these are `unique(A)`.

# Example

```julia
using GeometryOps
multipolygon = polygonize(>(0.6), rand(100, 100); minpoints=3)

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

function updateval(dict, key, val)
    if haskey(dict, key)
        existingval = dict[key][1][1]
        newval = ((existingval, val), (true, true))
        dict[key] = newval 
    else
        newval = ((val, val), (true, false))
        dict[key] = newval 
    end
end

function polygonize(f, xs::AbstractVector{T}, ys::AbstractVector{T}, A::AbstractMatrix; 
    minpoints=0,
) where T
    # Define buffers for edges and rings
    edges = Dict{T,Tuple{Tuple{T,T},Tuple{Bool,Bool}}}()
    rings = Vector{T}[]

    @assert (length(xs), length(ys)) == size(A)


    # First we collect all the edges around target pixels
    fi, fj = map(first, axes(A))
    li, lj = map(last, axes(A))
    @inbounds for i in axes(A, 1), j in axes(A, 2)
        if f(A[i, j]) # This is a pixel inside a polygon
            # xs and ys hold pixel bounds
            x1, x2 = xs[i]
            y1, y2 = ys[j]

            # We check the Von Neumann neighborhood to
            # decide what edges are needed, if any.
            (j == fi || !f(A[i, j-1])) && updateval(edges, (x1, y1), (x2, y1))
            (i == fj || !f(A[i-1, j])) && updateval(edges, (x1, y2), (x1, y1))
            (j == lj || !f(A[i, j+1])) && updateval(edges, (x2, y2), (x1, y2))
            (i == li || !f(A[i+1, j])) && updateval(edges, (x2, y1), (x2, y2))
        end
    end

    # Keep dict keys separately in a vector for performance
    edgekeys = collect(keys(edges))
    # We don't delete keys we just reduce length with nkeys
    nkeys = length(edgekeys)

    # Now create rings from the edges, 
    # looping until there are no edge keys left
    while nkeys > 0
        found = false
        local firstpoint, nextpoints, pointstatus

        # Loop until we find a key that hasn't been removed,
        # decrementing nkeys as we go.
        while nkeys > 0
            # Take the first edge from the array
            firstpoint::T = edgekeys[nkeys]
            nextpoints, pointstatus = edges[firstpoint]
            if any(pointstatus)
                found = true
                break
            else
                nkeys -= 1
            end
        end

        # If we found nothing this time, we are done
        found == false && break

        # Check if there are one or two lines going through this node
        # and take one of them, then update the status
        if pointstatus[2]
            nextpoint = nextpoints[2]
            edges[firstpoint] = ((nextpoints[1], map(zero, nextpoint)), (true, false))
        else
            nkeys -= 1
            nextpoint = nextpoints[1]
            edges[firstpoint] = ((map(zero, nextpoint), map(zero, nextpoint)), (false, false))
        end
        currentpoint = firstpoint
        ring = T[currentpoint, nextpoint]
        push!(rings, ring)
        # println()
        # @show currentpoint, nextpoint, pointstatus
        
        # Loop until we close a the ring and break
        while true
            # Find an edge that matches the next point
            (c1, c2), pointstatus = edges[nextpoint]
            # @show c1, c2, pointstatus
            # When there are two possible edges, 
            # choose the edge that has turned the furthest right
            if pointstatus[2]
                selectedpoint, remainingpoint = if currentpoint[1] == nextpoint[1] # vertical
                    wasincreasing = nextpoint[2] > currentpoint[2]
                    firstisstraight = nextpoint[1] == c1[1]
                    firstisleft = nextpoint[1] < c1[1]
                    if firstisstraight
                        secondisleft = nextpoint[1] > c2[1]
                        xor(wasincreasing, secondisleft) ? (c1, c2) : (c2, c1)
                    elseif firstisleft
                        wasincreasing ? (c2, c1) : (c1, c2)
                    else # firstisright
                        wasincreasing ? (c1, c2) : (c2, c1)
                    end
                else # horizontal
                    wasincreasing = nextpoint[1] > currentpoint[1]
                    firstisstraight = nextpoint[2] == c1[2]
                    firstisleft = nextpoint[2] < c1[2]
                    if firstisstraight
                        secondisleft = nextpoint[2] > c2[2]
                        xor(wasincreasing, secondisleft) ? (c1, c2) : (c2, c1)
                    elseif firstisleft
                        wasincreasing ? (c2, c1) : (c1, c2)
                    else # firstisright
                        wasincreasing ? (c1, c2) : (c2, c1)
                    end
                end
                edges[nextpoint] = ((remainingpoint, map(zero, remainingpoint)), (true, false))
                currentpoint, nextpoint = nextpoint, selectedpoint
            else
                edges[nextpoint] = ((map(zero, c1), map(zero, c1)), (false, false))
                currentpoint, nextpoint = nextpoint, c1
                # Write empty points, they are cleaned up later
            end
            # @show currentpoint, nextpoint, pointstatus
            # Close the ring if we get to the start
            if nextpoint == firstpoint
                push!(ring, nextpoint)
                break
            else
                i = findfirst(==(nextpoint), ring)
                if !isnothing(i)
                    # We found a touching point in the middle, 
                    # so we need to split the ring into two rings
                    splitring = ring[i:lastindex(ring)]
                    deleteat!(ring, i:lastindex(ring))
                    push!(splitring, nextpoint)
                    push!(rings, splitring)
                end
                push!(ring, nextpoint)
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
    foundholes = 0
    for poly in polygons
        exterior = GI.Polygon(StaticArrays.SVector(GI.getexterior(poly)))
        for i in eachindex(holes)
            unused[i] || continue
            hole = holes[i]
            if covers(poly, hole)
                foundholes += 1
                # Hole is in the exterior, so add it to the ring list
                push!(poly.geom, hole)
                # remove i
                unused[i] = false
                break
            end
        end
    end
    @show foundholes length(holes) length(polygons)

    # holepolygons = map(view(holes, unused)) do lr
    #     GI.Polygon([lr]; extent=GI.extent(lr))
    # end
    # append!(polygons, holepolygons) 
    # @assert foundholes == length(holes)

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

