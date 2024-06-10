# # Polygonizing raster data

export polygonize

#=
The methods in this file convert a raster image into a set of polygons, 
by contour detection using a clockwise Moore neighborhood method.

The resulting polygons are snapped to the boundaries of the cells of the input raster,
so they will look different from traditional contours from a plotting package.

The main entry point is the [`polygonize`](@ref) function.

```@docs
polygonize
```

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

Returns a `MultiPolygon` for `Bool` values and `f` return values, and
a `FeatureCollection` of `Feature`s holding `MultiPolygon` for all other values.


Function `f` should return either `true` or `false` or a transformation
of values into simpler groups, especially useful for floating point arrays.

If `xs` and `ys` are ranges, they are used as the pixel/cell center points.
If they are `Vector` of `Tuple` they are used as the lower and upper bounds of each pixel/cell.

# Keywords

- `minpoints`: ignore polygons with less than `minpoints` points.
- `values`: the values to turn into polygons. By default these are `union(A)`,
    If function `f` is passed these refer to the return values of `f`, by
    default `union(map(f, A)`. If values `Bool`, false is ignored and a single
    `MultiPolygon` is returned rather than a `FeatureCollection`.

# Example

```julia
using GeometryOps
A = rand(100, 100)
multipolygon = polygonize(>(0.5), A);
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
    xvec = similar(Vector{Tuple{Tx,Tx}}, xs)
    yvec = similar(Vector{Tuple{Ty,Ty}}, ys)
    for (xind, i) in enumerate(eachindex(xvec))
        xvec[i] = xbounds[xind], xbounds[xind+1]
    end
    for (yind, i) in enumerate(eachindex(yvec))
        yvec[i] = ybounds[yind], ybounds[yind+1]
    end
    return _polygonize(f, xvec, yvec, A; kw...)
end
function _polygonize(f, xs::AbstractVector{T}, ys::AbstractVector{T}, A::AbstractMatrix; 
    minpoints=0,
) where T<:Tuple
    (length(xs), length(ys)) == size(A) || throw(ArgumentError("length of xs and ys must match the array size"))

    # Extract the CRS of the array (if it is some kind of geo array / raster)
    crs = GI.crs(A)
    # Define buffers for edges and rings
    rings = Vector{T}[]

    strait = true
    turning = false

    # Get edges from the array A
    edges = _pixel_edges(f, xs, ys, A)
    # Keep dict keys separately in a vector for performance
    edgekeys = collect(keys(edges))
    # We don't delete keys we just reduce length with nkeys
    nkeys = length(edgekeys)

    # Now create rings from the edges, 
    # looping until there are no edge keys left
    while nkeys > 0
        found = false
        local firstnode, nextnodes, nodestatus

        map_partial(x,y) = map(!=(typemax(first(x))) ∘ first, y)
        # Loop until we find a key that hasn't been removed,
        # decrementing nkeys as we go.
        while nkeys > 0
            # Take the first node from the array
            firstnode::T = edgekeys[nkeys]
            nextnodes = edges[firstnode]
            nodestatus = map_partial(firstnode, nextnodes)
            if any(nodestatus)
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
        if nodestatus[2]
            nextnode = nextnodes[2]
            edges[firstnode] = (nextnodes[1], map(typemax, nextnode))
        else
            nkeys -= 1
            nextnode = nextnodes[1]
            edges[firstnode] = (map(typemax, nextnode), map(typemax, nextnode))
        end

        # Start a new ring
        currentnode = firstnode
        ring = [currentnode, nextnode]
        push!(rings, ring)
        
        # Loop until we close a the ring and break
        while true
            # Find a node that matches the next node
            (c1, c2) = possiblenodes = edges[nextnode]
            nodestatus = map_partial(firstnode, possiblenodes)
            if nodestatus[2]
                # When there are two possible node, 
                # choose the node that is the furthest to the left
                # We also need to check if we are on a straight line
                # to avoid adding unnecessary points.
                selectednode, remainingnode, straightline = if currentnode[1] == nextnode[1] # vertical
                    wasincreasing = nextnode[2] > currentnode[2]
                    firstisstraight = nextnode[1] == c1[1]
                    firstisleft = nextnode[1] > c1[1]
                    secondisstraight = nextnode[1] == c2[1]
                    secondisleft = nextnode[1] > c2[1]
                    if firstisstraight
                        if secondisleft 
                            if wasincreasing 
                                (c2, c1, turning)
                            else
                                (c1, c2, straight)
                            end
                        else
                            if wasincreasing 
                                (c1, c2, straight)
                            else
                                (c2, c1, secondisstraight)
                            end
                        end
                    elseif firstisleft
                        if wasincreasing 
                            (c1, c2, turning)
                        else
                            (c2, c1, secondisstraight)
                        end
                    else # firstisright
                        if wasincreasing 
                            (c2, c1, secondisstraight)
                        else
                            (c1, c2, turning)
                        end
                    end
                else # horizontal
                    wasincreasing = nextnode[1] > currentnode[1]
                    firstisstraight = nextnode[2] == c1[2]
                    firstisleft = nextnode[2] > c1[2]
                    secondisleft = nextnode[2] > c2[2]
                    secondisstraight = nextnode[2] == c2[2]
                    if firstisstraight
                        if secondisleft 
                            if wasincreasing 
                                (c1, c2, straight)
                            else
                                (c2, c1, turning)
                            end
                        else
                            if wasincreasing 
                                (c2, c1, turning)
                            else
                                (c1, c2, straight)
                            end
                        end
                    elseif firstisleft
                        if wasincreasing 
                            (c2, c1, secondisstraight)
                        else
                            (c1, c2, turning)
                        end
                    else # firstisright
                        if wasincreasing 
                            (c1, c2, turning)
                        else
                            (c2, c1, secondisstraight)
                        end
                    end
                end
                # Update edges
                edges[nextnode] = (remainingnode, map(typemax, remainingnode))
            else
                # Here we simply choose the first (and only valid) node
                selectednode = c1
                # Replace the edge nodes with empty nodes, they will be skipped later
                edges[nextnode] = (map(typemax, c1), map(typemax, c1))
                # Check if we are on a straight line
                straightline = currentnode[1] == nextnode[1] == c1[1] || 
                               currentnode[2] == nextnode[2] == c1[2]
            end

            # Update the current and next nodes with the next and selected nodes
            currentnode, nextnode = nextnode, selectednode
            # Update the current node or add a new node to the ring 
            if straightline
                # replace the last node we don't need it
                ring[end] = nextnode
            else
                # add a new node, we have turned a corner
                push!(ring, nextnode)
            end
            # If the ring is closed, break the loop and start a new one
            nextnode == firstnode && break
        end
    end

    # Define wrapped LinearRings, with embedded extents
    # so we only calculate them once
    linearrings = map(rings) do ring
        extent = GI.extent(GI.LinearRing(ring))
        GI.LinearRing(ring; extent, crs)
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
        GI.Polygon([lr]; extent=GI.extent(lr), crs)
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
        return GI.MultiPolygon(polygons; crs, extent = mapreduce(GI.extent, Extents.union, polygons))
    end
end

function _polygonize_featurecollection(f::Base.Callable, xs::AbstractRange, ys::AbstractRange, A::AbstractMatrix; 
    values=_default_values(f, A), kw...
)
    crs = GI.crs(A)
    # Create one feature per value
    features = map(values) do value
        multipolygon = _polygonize(x -> isequal(f(x), value), xs, ys, A; kw...)
        GI.Feature(multipolygon; properties=(; value), extent = GI.extent(multipolygon), crs)
    end 

    return GI.FeatureCollection(features; extent = mapreduce(GI.extent, Extents.union, features), crs)
end

function _default_values(f, A)
    # Get union of f return values with resolved eltype
    values = map(identity, sort!(Base.union(Iterators.map(f, A))))
    # We ignore pure Bool
    return eltype(values) == Bool ? nothing : collect(skipmissing(values))
end

function update_edge!(dict, key, node)
    newnodes = (node, map(typemax, node))
    # Get or write in one go, to skip a hash lookup
    existingnodes = get!(() -> newnodes, dict, key)
    # If we actually fetched an existing node, update it
    if existingnodes[1] != node
        dict[key] = (existingnodes[1], node)
    end
end

function _pixel_edges(f, xs::AbstractVector{T}, ys::AbstractVector{T}, A) where T<:Tuple
    edges = Dict{T,Tuple{T,T}}()
    # First we collect all the edges around target pixels
    fi, fj = map(first, axes(A))
    li, lj = map(last, axes(A))
    for j in axes(A, 2)
        y1, y2 = ys[j]
        for i in axes(A, 1) 
            if f(A[i, j]) # This is a pixel inside a polygon
                # xs and ys hold pixel bounds
                x1, x2 = xs[i]
                # We check the Von Neumann neighborhood to
                # decide what edges are needed, if any.
                (j == fi || !f(A[i, j-1])) && update_edge!(edges, (x1, y1), (x2, y1)) # S
                (i == fj || !f(A[i-1, j])) && update_edge!(edges, (x1, y2), (x1, y1)) # W
                (j == lj || !f(A[i, j+1])) && update_edge!(edges, (x2, y2), (x1, y2)) # N
                (i == li || !f(A[i+1, j])) && update_edge!(edges, (x2, y1), (x2, y2)) # E
            end
        end
    end
    return edges
end


