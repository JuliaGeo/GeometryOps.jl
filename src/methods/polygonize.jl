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
function polygonize(xs::AbstractRange, ys::AbstractRange, A::AbstractMatrix{Bool}; minpoints=3)
    # Define buffers for edges and rings
    edges = Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}[] 
    rings = Vector{Tuple{Float64,Float64}}[]

    # First get all the valid edges between 
    # filled an empty pixels
    halfxstep, halfystep = step(xs) / 2, step(ys) /2
    si, sj = map(last, axes(A))
    for i in axes(A, 1), j in axes(A, 2)
        if A[i, j] # This is a pixel inside a polygon
            # xs and ys hold pixel centers
            x, y = xs[i], ys[j]
            # So we need to offset the edges by half
            x1, y1 = x - halfxstep, y - halfystep
            x2, y2 = x + halfxstep, y + halfystep
            # Then we check the Von Neumann neighborhood to
            # decide what edges are needed, if any.
            j >= 1 && !A[i, j-1] && push!(edges, ((x1, y1), (x2, y1)))
            i >= 1 && !A[i-1, j] && push!(edges, ((x1, y2), (x1, y1)))
            j <= sj && !A[i, j+1] && push!(edges, ((x2, y2), (x1, y2)))
            i <= si && !A[i+1, j] && push!(edges, ((x2, y1), (x2, y2)))
        end
    end

    # Then we join the edges into rings
    sort!(edges)
    while length(edges) > 0
        edge = pop!(edges)
        firstpoint::Tuple{Float64,Float64} = first(edge)
        nextpoint::Tuple{Float64,Float64} = last(edge)
        ring = [firstpoint, nextpoint]
        push!(rings, ring)
        while true
            # Find an edge that matches the next point
            i = searchsortedlast(edges, (nextpoint, nextpoint); by=first)
            edge = popat!(edges, i)
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
    clockwise = isclockwise.(linearrings)
    polygons = map(x -> [x], linearrings[.!clockwise])
    holes = map(x -> x, linearrings[clockwise])

    # Then we add holes to the polygons
    blacklist = Set{Int}()
    polylist = Set{Int}()
    foreach(polygons) do rings
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
end
