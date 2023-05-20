# # Barycentric coordinates
# Generalized barycentric coordinates are a generalization of barycentric coordinates, 
# which are typically used in triangles, to arbitrary polygons. 

# They provide a way to express a point within a polygon as a weighted average 
# of the polygon's vertices.

# In the case of a triangle, barycentric coordinates are a set of three numbers 
# $(λ_1, λ_2, λ_3)$, each associated with a vertex of the triangle. Any point within 
# the triangle can be expressed as a weighted average of the vertices, where the 
# weights are the barycentric coordinates. The weights sum to 1, and each is non-negative.

# For a polygon with n vertices, generalized barycentric coordinates are a set of 
# $n$ numbers $(λ_1, λ_2, ..., λ_n)$, each associated with a vertex of the polygon. 
# Any point within the polygon can be expressed as a weighted average of the vertices, 
# where the weights are the generalized barycentric coordinates. 

# As with the triangle case, the weights sum to 1, and each is non-negative.

# ## Methods to find barycentric coordinates

function mean_value_barycentric_coordinates(polypoints::Vector{<: Point{N1, T1}}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real}
    n = length(polypoints)
    λ = zeros(promote_type(T1, T2), n)
    if (N1 == N2 && T1 == T2)
        poly = polypoints
    else
        poly = Makie.to_ndim.((Point{N2, promote_type(T1, T2)},), polypoints, 0)
    end
    # if !GeometryOps.contains(poly, point)
    #     return λ
    # end
    # Describe this loop
    # The loop computes barycentric coordinates by the mean-value method.
    # The mean-value method is a method for computing barycentric coordinates
    # for a point in a polygon. It is based on the observation that the
    # barycentric coordinates of a point in a polygon are proportional to the
    # areas of the triangles formed by the point and each pair of edges of the
    # polygon. The mean-value method computes the areas of these triangles by
    # computing the areas of the triangles formed by the point and each pair of
    # edges of the polygon, and then averaging these areas.
    for i in 1:n

        prev = poly[mod1(i-2, n)]
        curr = poly[i]
        next = poly[mod1(i, n)]

        α1 = angle(prev, point, curr)
        α2 = angle(curr, point, next)
        d1 = distance(point, curr)
        λ[i] = (tan(α1 / 2) + tan(α2 / 2)) / abs(d1)

    end
    λ /= sum(λ)
    return λ
end

function angle(a, b, c)
    ab = a - b
    cb = c - b
    acos_param = dot(ab, cb) / (norm(ab) * norm(cb))
    return if abs(acos_param) > 1
        0.0
    else
        acos(acos_param)
    end
end

function distance(a, b)
    return norm(a - b)
end

# ## Example
# This example was taken from [this page of CGAL's documentation](https://doc.cgal.org/latest/Barycentric_coordinates_2/index.html).

# ```@example barycentric
# polygon_points = [
# Point(0.03, 0.05, 0.00), Point(0.07, 0.04, 0.02), Point(0.10, 0.04, 0.04),
# Point(0.14, 0.04, 0.06), Point(0.17, 0.07, 0.08), Point(0.20, 0.09, 0.10),
# Point(0.22, 0.11, 0.12), Point(0.25, 0.11, 0.14), Point(0.27, 0.10, 0.16),
# Point(0.30, 0.07, 0.18), Point(0.31, 0.04, 0.20), Point(0.34, 0.03, 0.22),
# Point(0.37, 0.02, 0.24), Point(0.40, 0.03, 0.26), Point(0.42, 0.04, 0.28),
# Point(0.44, 0.07, 0.30), Point(0.45, 0.10, 0.32), Point(0.46, 0.13, 0.34),
# Point(0.46, 0.19, 0.36), Point(0.47, 0.26, 0.38), Point(0.47, 0.31, 0.40),
# Point(0.47, 0.35, 0.42), Point(0.45, 0.37, 0.44), Point(0.41, 0.38, 0.46),
# Point(0.38, 0.37, 0.48), Point(0.35, 0.36, 0.50), Point(0.32, 0.35, 0.52),
# Point(0.30, 0.37, 0.54), Point(0.28, 0.39, 0.56), Point(0.25, 0.40, 0.58),
# Point(0.23, 0.39, 0.60), Point(0.21, 0.37, 0.62), Point(0.21, 0.34, 0.64),
# Point(0.23, 0.32, 0.66), Point(0.24, 0.29, 0.68), Point(0.27, 0.24, 0.70),
# Point(0.29, 0.21, 0.72), Point(0.29, 0.18, 0.74), Point(0.26, 0.16, 0.76),
# Point(0.24, 0.17, 0.78), Point(0.23, 0.19, 0.80), Point(0.24, 0.22, 0.82),
# Point(0.24, 0.25, 0.84), Point(0.21, 0.26, 0.86), Point(0.17, 0.26, 0.88),
# Point(0.12, 0.24, 0.90), Point(0.07, 0.20, 0.92), Point(0.03, 0.15, 0.94),
# Point(0.01, 0.10, 0.97), Point(0.02, 0.07, 1.00)]
# f, a, p = poly(polygon_points; color = last.(polygon_points), colormap = cgrad(:jet, 18; categorical = true), shading = false)
# ```
# You can see that the polygon triangulation doesn't do any justice to the actual structure.  Let's visualize it using mean value coordinates:
# ```@example barycentric
# polygon_zs = last.(polygon_points)
# xrange = LinRange(0, 0.5, 1000)
# yrange = LinRange(0, 0.5, 1000)
# @time mean_value_coordinate_field = mean_value_barycentric_coordinates.(
#     (Point2{Float64}.(polygon_points),), 
#     Point2{Float64}.(xrange, yrange')
# )
# ```
# This returned a list of vectors of weights, one for each point in the grid.  

# We can use these vectors to extrapolate and interpolate from the z values of the polygon to some z value at each point in the grid:
# ```@example barycentric
# mean_values = map(mean_value_coordinate_field) do λs
#     sum(λs .* polygon_zs)
# end
# ```
# Now, we can visualize this!
# ```@example barycentric
# f, a, p = heatmap(xrange, yrange, mean_values; colormap = cgrad(:jet, 18; categorical = true))
# poly!(Point2f.(polygon_points); color = :transparent, strokecolor = :black, strokewidth = 1.3)
# Colorbar(f[1, 2], p)
# f
# ```