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
            A[i, j+1] || push!(edges, ((x1, y2), (x2, y2)))
            A[i-1, j] || push!(edges, ((x1, y1), (x1, y2)))
            A[i+1, j] || push!(edges, ((x2, y1), (x2, y2)))
        end
    end
    # Then join them into polygons
    ring = Tuple{Float64,Float64}[]
    while length(edges) > 0
        edge = pop!(edges)
        firstpoint = first(edge)
        nextpoint = last(edge)
        ring = [firstpoint, nextpoint]
        push!(rings, ring)
        while true
            # Find an edge that matches the next point
            i = findfirst(e -> nextpoint in e, edges)
            edge = popat!(edges, i)
            nextpoint = otherpoint(edge, nextpoint)
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
    polygons = map(x -> [GI.LinearRing(x)], rings)
    while i <= length(rings)
        push!(blacklist, i)
        ring_i = GI.LinearRing(rings[i]) 
        for j in 1:length(rings)
            j in blacklist && continue
            ring_j = GI.LinearRing(rings[j]) 
            if polygon_in_polygon(ring_j, ring_i)
                push!(blacklist, j)
                push!(polygons[i], ring_j)
                push!(polylist, i)
            elseif polygon_in_polygon(ring_i, rings_j)
                push!(blacklist, j)
                push!(polygons[j], ring_i)
                push!(polylist, j)
            end
        end
    end
    return GI.Polygon.(polygons[sort!(collect(polylist))])
end

otherpoint(edge, point) = first(edge) == point ? last(edge) : first(edge)

_position(mode, dir) = _deltas(mode)[dir]


## rotate direction clockwise
rot_clockwise(dir) = dir % 4 + 1
## rotate direction counterclockwise
rot_counterclockwise(dir) = (dir + 2) % 4 + 1
