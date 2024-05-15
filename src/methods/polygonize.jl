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

Returns a `MultiPolygon` for `Bool` values and `f` return values, and
a `FeatureCollection` of `Feature`s holding `MultiPolygon` for all other values.

Function `f` should return either `true` or `false` or a transformation
of values into simpler groups, especially useful for floating point arrays.

If `xs` and `ys` are ranges, they are used as the pixel center points.
If they are `Vector` of `Tuple` they are used as the lower and upper bounds of each pixel.

# Keywords

- `minpoints`: ignore polygons with less than `minpoints` points.
- `values`: the values to turn into polygons. By default these are `union(A)`,
    If function `f` is passed these refer to the return values of `f`, by
    default `union(map(f, A)`. If values `Bool`, false is ignored and a single
    `MultiPolygon` is returned rather than a `FeatureCollection`.

# Example

```julia
using GeometryOps, GLMakie, Shapefile
A = rand(100, 100)
multipolygon = polygonize(>(0.5), A);

p = Makie.plot(multipolygon.geom[1])
Makie.plot!.(multipolygon.geom[2:end])
```
"""
polygonize(A::AbstractMatrix{Bool}; kw...) = polygonize(identity, A; kw...)
polygonize(f::Base.Callable, A::AbstractMatrix; kw...) = polygonize(f, axes(A)..., A; kw...)
polygonize(A::AbstractMatrix; kw...) = polygonize(axes(A)..., A; kw...)
polygonize(xs::AbstractVector, ys::AbstractVector, A::AbstractMatrix{Bool}; kw...) =
    _polygonize(identity, xs, ys, A)
function polygonize(xs::AbstractVector, ys::AbstractVector, A::AbstractMatrix; 
    values=sort!(Base.union(A)), kw...
)
    _polygonize_featurecollection(identity, xs, ys, A; values, kw...) 
end
function polygonize(f::Base.Callable, xs::AbstractRange, ys::AbstractRange, A::AbstractMatrix; 
    values=_default_values(f, A), kw...
)
    if isnothing(values)
        _polygonize(f, xs, ys, A; kw...) 
    else
        _polygonize_featurecollection(f, xs, ys, A; kw...) 
    end
end

function _polygonize_featurecollection(f::Base.Callable, xs::AbstractRange, ys::AbstractRange, A::AbstractMatrix; 
    values=_default_values(f, A), kw...
)
    # Create one feature per value
    features = map(values) do value
        multipolygon = _polygonize(x -> isequal(f(x), value), xs, ys, A; kw...)
        GI.Feature(multipolygon; properties=(; value))
    end 

    return GI.FeatureCollection(features)
end

function _polygonize(f::Base.Callable, xs::AbstractRange, ys::AbstractRange, A::AbstractMatrix; 
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
    return _polygonize(f, xvec, yvec, A; kw...)
end

function _default_values(f, A)
    # Get union of f return values with resolved eltype
    values = map(identity, sort!(Base.union(Iterators.map(f, A))))
    # We ignore pure Bool
    return eltype(values) == Bool ? nothing : collect(skipmissing(values))
end

@inline function updateval(dict, key, val)
    if haskey(dict, key)
        existingvals = dict[key]
        existingval = existingvals[1]
        newval = (existingval, val)
        dict[key] = newval 
    else
        newval = (val, map(typemax, val))
        dict[key] = newval 
    end
end

function _polygonize(f, xs::AbstractVector{T}, ys::AbstractVector{T}, A::AbstractMatrix; 
    minpoints=0,
) where T
    # Define buffers for edges and rings
    edges = Dict{T,Tuple{T,T}}()
    rings = Vector{T}[]

    (length(xs), length(ys)) == size(A) || throw(ArgumentError("length of xs and ys must match the array size"))

    # First we collect all the edges around target pixels
    fi, fj = map(first, axes(A))
    li, lj = map(last, axes(A))
    for i in axes(A, 1), j in axes(A, 2)
        if f(A[i, j]) # This is a pixel inside a polygon
            # xs and ys hold pixel bounds
            x1, x2 = xs[i]
            y1, y2 = ys[j]

            # We check the Von Neumann neighborhood to
            # decide what edges are needed, if any.
            (j == fi || !f(A[i, j-1])) && updateval(edges, (x1, y1), (x2, y1)) # S
            (i == fj || !f(A[i-1, j])) && updateval(edges, (x1, y2), (x1, y1)) # W
            (j == lj || !f(A[i, j+1])) && updateval(edges, (x2, y2), (x1, y2)) # N
            (i == li || !f(A[i+1, j])) && updateval(edges, (x2, y1), (x2, y2)) # E
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
            nextpoints = edges[firstpoint]
            pointstatus = map(!=(typemax(first(firstpoint))) ∘ first, nextpoints)
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
            edges[firstpoint] = (nextpoints[1], map(typemax, nextpoint))
        else
            nkeys -= 1
            nextpoint = nextpoints[1]
            edges[firstpoint] = (map(typemax, nextpoint), map(typemax, nextpoint))
        end

        # Start a new ring
        currentpoint = firstpoint
        ring = [currentpoint, nextpoint]
        push!(rings, ring)
        
        # Loop until we close a the ring and break
        while true
            # Find an edge that matches the next point
            (c1, c2) = possiblepoints = edges[nextpoint]
            pointstatus = map(!=(typemax(first(firstpoint))) ∘ first, possiblepoints)
            # When there are two possible edges, 
            # choose the edge that has turned the furthest right
            if pointstatus[2]
                selectedpoint, remainingpoint = if currentpoint[1] == nextpoint[1] # vertical
                    wasincreasing = nextpoint[2] > currentpoint[2]
                    firstisstraight = nextpoint[1] == c1[1]
                    firstisleft = nextpoint[1] > c1[1]
                    if firstisstraight
                        secondisleft = nextpoint[1] > c2[1]
                        if secondisleft 
                            wasincreasing ? (c2, c1) : (c1, c2)
                        else
                            wasincreasing ? (c1, c2) : (c2, c1) 
                        end
                    elseif firstisleft
                        wasincreasing ? (c1, c2) : (c2, c1)
                    else # firstisright
                        wasincreasing ? (c2, c1) : (c1, c2)
                    end
                else # horizontal
                    wasincreasing = nextpoint[1] > currentpoint[1]
                    firstisstraight = nextpoint[2] == c1[2]
                    firstisleft = nextpoint[2] > c1[2]
                    if firstisstraight
                        secondisleft = nextpoint[2] > c2[2]
                        if secondisleft 
                            wasincreasing ? (c1, c2) : (c2, c1) 
                        else
                            wasincreasing ? (c2, c1) : (c1, c2)
                        end
                    elseif firstisleft
                        wasincreasing ? (c2, c1) : (c1, c2)
                    else # firstisright
                        wasincreasing ? (c1, c2) : (c2, c1)
                    end
                end

                # Update edges
                edges[nextpoint] = (remainingpoint, map(typemax, remainingpoint))
                currentpoint, nextpoint = nextpoint, selectedpoint
            else
                # Update edges
                edges[nextpoint] = (map(typemax, c1), map(typemax, c1))
                currentpoint, nextpoint = nextpoint, c1
                # Write empty points, they are cleaned up later
            end
            push!(ring, nextpoint)
            # If the ring is closed, start a new one
            nextpoint == firstpoint && break
        end
    end

    # Define wrapped LinearRings, with embedded extents
    # so we only calculate them once
    linearrings = map(rings) do ring
        extent = GI.extent(GI.LinearRing(ring))
        GI.LinearRing(ring; extent)
    end

    # Separate exteriors from holes by winding direction
    direction = (last(last(xs)) - first(first(xs))) * (last(last(ys)) - first(first(ys)))
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
    assigned = fill(false, length(holes))
    for i in eachindex(holes)
        hole = holes[i]
        prepared_hole = GI.LinearRing(holes[i]; extent=GI.extent(holes[i]))
        for poly in polygons
            exterior = GI.Polygon(StaticArrays.SVector(GI.getexterior(poly)); extent=GI.extent(poly))
            if covers(exterior, prepared_hole)
                # Hole is in the exterior, so add it to the polygon
                push!(poly.geom, hole)
                assigned[i] = true
                break
            end
        end
    end

    assigned_holes = count(assigned)
    assigned_holes == length(holes) || @warn "Not all holes were assigned to polygons, $(length(holes) - assigned_holes) where missed from $(length(holes)) holes and $(length(polygons)) polygons"

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

