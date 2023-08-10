# # Polygonizing raster data
export polygonize
# The methods in this file are able to convert a raster image into a set of polygons,
# by contour detection using a clockwise Moore neighborhood method.

# The main entry point is the [`polygonize`](@ref) function.

# ```@doc
# polygonize
# ```

# ## Example

# Here's a basic implementation, using the `Makie.peaks()` function.  First, let's investigate the nature of the function:
# ```@example polygonize
# using Makie, GeometryOps
# n = 49
# xs, ys = LinRange(-3, 3, n), LinRange(-3, 3, n)
# zs = Makie.peaks(n)
# z_max_value = maximum(abs.(extrema(zs)))
# f, a, p = heatmap(
#     xs, ys, zs;
#     axis = (; aspect = DataAspect(), title = "Exact function")
# )
# cb = Colorbar(f[1, 2], p; label = "Z-value")
# f
# ```

# Now, we can use the `polygonize` function to convert the raster data into polygons.

# For this particular example, we chose a range of z-values between 0.8 and 3.2,
# which would provide two distinct polyogns with holes.

# ```@example polygonize
# polygons = polygonize(xs, ys, 0.8 .< zs .< 3.2)
# ```
# This returns a list of `GeometryBasics.Polygon`, which can be plotted immediately,
# or wrapped directly in a `GeometryBasics.MultiPolygon`.  Let's see how these look:

# ```@example polygonize
# f, a, p = poly(polygons; label = "Polygonized polygons", axis = (; aspect = DataAspect()))
# ```

# Finally, let's plot the Makie contour lines on top, to see how well the polygonization worked:
# ```@example polygonize
# contour!(a, zs; labels = true, levels = [0.8, 3.2], label = "Contour lines")
# f
# ```

# ## Implementation

# The implementation follows:

"""
    polygonize(A; minpoints=10)
    polygonize(xs, ys, A; minpoints=10)

Convert matrix `A` to polygons.

If `xs` and `ys` are passed in they are used as the pixel center points.

# Keywords
- `minpoints`: ignore polygons with less than `minpoints` points.
"""
polygonize(A::AbstractMatrix; kw...) = polygonize(axes(A)..., A; kw...)
function polygonize(xs::AbstractRange, ys::AbstractRange, A::AbstractMatrix{Bool})
    # Make ranges of pixel bounds
    xbounds = first(xs) - step(xs) / 2 : step(xs) : last(xs) + step(xs) / 2 
    ybounds = first(ys) - step(ys) / 2 : step(ys) : last(ys) + step(ys) / 2 
    return _polygonize(xbounds, ybounds, A)
end
function polygonize(xs::AbstractVector, ys::AbstractVector, A::AbstractMatrix{Bool})
    _polygonize(xs, ys, A)
end

function _polygonize(xs::AbstractVector{T}, ys::AbstractVector{T}, A::AbstractMatrix{Bool}) where T
    # Define buffers for edges and rings
    edges = Tuple{Tuple{T,T},Tuple{T,T}}[]
    rings = Vector{Tuple{T,T}}[]

    # First get all the valid edges between filled and empty pixels
    si, sj = map(last, axes(A))
    for i in axes(A, 1), j in axes(A, 2)
        if A[i, j] # This is a pixel inside a polygon
            # xs and ys hold pixel bounds
            x1, x2 = xs[i], xs[i + 1] 
            y1, y2 = ys[j], ys[j + 1] 
            # We check the Von Neumann neighborhood to
            # decide what edges are needed, if any.
            (j == 1 || !A[i, j-1]) && push!(edges, ((x1, y1), (x2, y1)))
            (i == 1 || !A[i-1, j]) && push!(edges, ((x1, y2), (x1, y1)))
            (j == sj || !A[i, j+1]) && push!(edges, ((x2, y2), (x1, y2)))
            (i == si || !A[i+1, j]) && push!(edges, ((x2, y1), (x2, y2)))
        end
    end

    # Then we join the edges into rings
    # Sorting now lets us use `searchsortedlast` for speed later
    sort!(edges)
    while length(edges) > 0
        # Take the last edge from the array
        edge = pop!(edges)
        firstpoint::Tuple{T,T} = first(edge)
        nextpoint::Tuple{T,T} = last(edge)
        ring = [firstpoint, nextpoint]
        push!(rings, ring)
        while length(edges) > 0
            # Find an edge that matches the next point
            i = searchsortedlast(edges, (nextpoint, nextpoint); by=first)
            newedge = edges[i]
            # When there are two possible edges, 
            # choose the edge that has turned a corner
            if (i > 1) && (otheredge = edges[i - 1]; otheredge[1] == newedge[1]) &&
                (edge[2][1] == newedge[2][1] || edge[2][2] == newedge[2][2]) 
                newedge = otheredge
                deleteat!(edges, i - 1)
            else
                deleteat!(edges, i)
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
            elseif nextpoint in ring
                i = findfirst(==(nextpoint), ring)
                splitring = ring[i:lastindex(ring)]
                deleteat!(ring, i:lastindex(ring))
                push!(ring, nextpoint)
                push!(splitring, nextpoint)
                push!(rings, splitring)
                continue
            # Otherwise keep adding points to the ring
            else
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
    direction = last(xs) - first(xs) * last(ys) - first(ys)
    exterior_inds = if direction > 0 
        .!isclockwise.(linearrings)
    else
        isclockwise.(linearrings)
    end
    holes = map(x -> x, linearrings[.!exterior_inds])
    exteriors = map(x -> x, linearrings[exterior_inds])
    polygons = map(x -> [x], linearrings[exterior_inds])

    # Then we add the holes to the polygons they are inside of
    blacklist = Set{Int}()
    for rings in polygons
        for i in eachindex(holes)
            i in blacklist && continue
            exterior = rings[1]
            hole = holes[i]
            if polygon_in_polygon(hole, exterior)
                # Hole is in the exterior, so add it to the ring list
                push!(rings, hole)
                # And blacklist it so we don't check it again
                push!(blacklist, i)
            end
        end
    end

    # Finally, return wrapped Polygons
    return GI.Polygon.(polygons)
    # return (x -> GI.Polygon([x])).(exteriors), (x -> GI.Polygon([x])).(holes)
end

