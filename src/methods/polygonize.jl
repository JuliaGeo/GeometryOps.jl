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
    edges = Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}[] 
    xstep, ystep = step(xs) / 2, step(ys) /2
    rings = Vector{Tuple{Float64,Float64}}[]
    # Get all the valid edges
    for i in 2:size(A, 1)-1, j in 2:size(A, 2)-1
        if A[i, j]
            x, y = xs[i], ys[j]
            x1, y1 = x - xstep, y - ystep
            x2, y2 = x + xstep, y + ystep
            A[i, j-1] || push!(edges, ((x1, y1), (x2, y1)))
            A[i-1, j] || push!(edges, ((x1, y2), (x1, y1)))
            A[i, j+1] || push!(edges, ((x2, y2), (x1, y2)))
            A[i+1, j] || push!(edges, ((x2, y1), (x2, y2)))
        end
    end
    # Then join them into polygons
    ring = Tuple{Float64,Float64}[]
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
            # Close if we get to the start
            if nextpoint == firstpoint
                push!(ring, nextpoint)
                break
            # Split if we touch the ring
            elseif nextpoint in ring
                i = findfirst(==(nextpoint), ring)
                splitring = ring[i:lastindex(ring)]
                deleteat!(ring, i:lastindex(ring))
                push!(ring, nextpoint)
                push!(splitring, nextpoint)
                push!(rings, splitring)
                continue
            # Otherwise keep adding points
            else
                push!(ring, nextpoint)
            end
        end
    end
    i = 1
    blacklist = Set{Int}()
    polylist = Set{Int}()
    lrings = map(rings) do ring
        extent = GI.extent(GI.LinearRing(ring))
        GI.LinearRing(ring; extent)
    end
    clk = isclockwise.(lrings)
    exteriors = map(x -> [x], lrings[.!clk])
    holes = map(x -> x, lrings[clk])
    blacklist = Set{Int}()
    foreach(exteriors) do e
        for i in eachindex(holes)
            i in blacklist && continue
            h = holes[i]
            if polygon_in_polygon(h, e[1])
                push!(e, h)
                push!(blacklist, i)
            end
        end
    end
    return GI.Polygon.(exteriors)
end
