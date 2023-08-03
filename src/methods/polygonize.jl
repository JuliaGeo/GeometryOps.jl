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
    xstep, ystep = step(xs), step(ys)
    rings = Vector{Tuple{Float64,Float64}}[]
    for i in 2:size(A, 1)-1, j in 2:size(A, 2)-1
        if A[i, j]
            x1, y1 = xs[i], ys[j]
            x2, y2 = x1 + xstep, y1 + ystep
            A[i, j-1] || push!(edges, ((x1, y1), (x2, y1)))
            A[i, j+1] || push!(edges, ((x1, y2), (x2, y2)))
            A[i-1, j] || push!(edges, ((x1, y1), (x1, y2)))
            A[i+1, j] || push!(edges, ((x2, y1), (x2, y2)))
        end
        @show edges
    end
    ring = Tuple{Float64,Float64}[]
    dir = 1
    while length(edges) > 0
        edge = pop!(edges)
        firstpoint = first(edge)
        nextpoint = last(edge)
        ring = [firstpoint, nextpoint]
        push!(rings, ring)
        while length(edges) > 0
            i = findfirst(e -> nextpoint in e, edges)
            isnothing(i) && break
            edge = popat!(edges, i)
            newpoint = otherpoint(edge, nextpoint)
            if newpoint == firstpoint
                push!(ring, newpoint)
                break
            # elseif newpoint in ring
            #     i = findfirst(==(newpoint), ring)
            #     splitring = ring[i:lastindex(ring)]
            #     deleteat!(ring, i:lastindex(ring))
            #     push!(splitring, newpoint)
            #     push!(rings, splitring)
            #     nextpoint = last(ring)
            #     continue
            else
                push!(ring, newpoint)
                nextpoint = newpoint
            end
        end
    end
    map(rings) do ring
        GB.Polygon(map(GB.Point, ring))
    end
end

otherpoint(edge, point) = first(edge) == point ? last(edge) : first(edge)

_position(mode, dir) = _deltas(mode)[dir]


## rotate direction clockwise
rot_clockwise(dir) = dir % 4 + 1
## rotate direction counterclockwise
rot_counterclockwise(dir) = (dir + 2) % 4 + 1
