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
    ## if !GeometryOps.contains(poly, point)
    ##     return λ
    ## end
    ## Describe this loop
    ## The loop computes barycentric coordinates by the mean-value method.
    ## The mean-value method is a method for computing barycentric coordinates
    ## for a point in a polygon. It is based on the observation that the
    ## barycentric coordinates of a point in a polygon are proportional to the
    ## areas of the triangles formed by the point and each pair of edges of the
    ## polygon. The mean-value method computes the areas of these triangles by
    ## computing the areas of the triangles formed by the point and each pair of
    ## edges of the polygon, and then averaging these areas.
    for i in 1:n

        prev = poly[mod1(i-1, n)]
        curr = poly[i]
        next = poly[mod1(i+1, n)]

        α1 = angle(prev, point, curr)
        α2 = angle(curr, point, next)
        d1 = distance(point, curr)
        λ[i] = (tan(α1 / 2) + tan(α2 / 2)) / abs(d1)

    end
    ## Normalize the vector to sum to 1
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

function mean_value_barycentric_coordinates(poly::Polygon{N1, T1}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real}
    ext = decompose(Point{N2, T2}, poly.exterior)
    ints = decompose.((Point{N2, T2},), poly.interiors)
    n = length(ext) + if isempty(ints)
        0
    else
        sum(length.(ints))
    end
    λ = zeros(promote_type(T1, T2), n)
    ## if !GeometryOps.contains(poly, point)
    ##     return λ
    ## end
    ## Describe this loop
    ## The loop computes barycentric coordinates by the mean-value method.
    ## The mean-value method is a method for computing barycentric coordinates
    ## for a point in a polygon. It is based on the observation that the
    ## barycentric coordinates of a point in a polygon are proportional to the
    ## areas of the triangles formed by the point and each pair of edges of the
    ## polygon. The mean-value method computes the areas of these triangles by
    ## computing the areas of the triangles formed by the point and each pair of
    ## edges of the polygon, and then averaging these areas.
    current_ind = 1
    l_ext = length(ext)
    for i in 1:l_ext

        prev = ext[mod1(i-1, l_ext)]
        curr = ext[i]
        next = ext[mod1(i+1, l_ext)]

        α1 = angle(prev, point, curr)
        α2 = angle(curr, point, next)
        d1 = distance(point, curr)
        λ[current_ind] = (tan(α1 / 2) + tan(α2 / 2)) / abs(d1)
        current_ind += 1

    end
    for hole in ints
        l_hole = length(hole)
        for i in 1:l_hole

            prev = hole[mod1(i-1, l_hole)]
            curr = hole[i]
            next = hole[mod1(i+1, l_hole)]

            α1 = angle(prev, point, curr)
            α2 = angle(curr, point, next)
            d1 = distance(point, curr)
            λ[current_ind] = (tan(α1 / 2) + tan(α2 / 2)) / abs(d1)
            current_ind += 1

        end
    end
    ## Normalize the vector to sum to 1
    λ /= sum(λ)
    return λ
end

# ## Example
# This example was taken from [this page of CGAL's documentation](https://doc.cgal.org/latest/Barycentric_coordinates_2/index.html).

# ```@example barycentric
# using GeometryBasics, GeometryOps, Makie
# polygon_points = Point3f[
# (0.03, 0.05, 0.00), (0.07, 0.04, 0.02), (0.10, 0.04, 0.04),
# (0.14, 0.04, 0.06), (0.17, 0.07, 0.08), (0.20, 0.09, 0.10),
# (0.22, 0.11, 0.12), (0.25, 0.11, 0.14), (0.27, 0.10, 0.16),
# (0.30, 0.07, 0.18), (0.31, 0.04, 0.20), (0.34, 0.03, 0.22),
# (0.37, 0.02, 0.24), (0.40, 0.03, 0.26), (0.42, 0.04, 0.28),
# (0.44, 0.07, 0.30), (0.45, 0.10, 0.32), (0.46, 0.13, 0.34),
# (0.46, 0.19, 0.36), (0.47, 0.26, 0.38), (0.47, 0.31, 0.40),
# (0.47, 0.35, 0.42), (0.45, 0.37, 0.44), (0.41, 0.38, 0.46),
# (0.38, 0.37, 0.48), (0.35, 0.36, 0.50), (0.32, 0.35, 0.52),
# (0.30, 0.37, 0.54), (0.28, 0.39, 0.56), (0.25, 0.40, 0.58),
# (0.23, 0.39, 0.60), (0.21, 0.37, 0.62), (0.21, 0.34, 0.64),
# (0.23, 0.32, 0.66), (0.24, 0.29, 0.68), (0.27, 0.24, 0.70),
# (0.29, 0.21, 0.72), (0.29, 0.18, 0.74), (0.26, 0.16, 0.76),
# (0.24, 0.17, 0.78), (0.23, 0.19, 0.80), (0.24, 0.22, 0.82),
# (0.24, 0.25, 0.84), (0.21, 0.26, 0.86), (0.17, 0.26, 0.88),
# (0.12, 0.24, 0.90), (0.07, 0.20, 0.92), (0.03, 0.15, 0.94),
# (0.01, 0.10, 0.97), (0.02, 0.07, 1.00)]
# f, a, p = poly(polygon_points; color = last.(polygon_points), colormap = cgrad(:jet, 18; categorical = true), shading = false)
# ```
# You can see that the polygon triangulation doesn't do any justice to the actual structure.  Let's visualize it using mean value coordinates:
# ```@example barycentric
# polygon_zs = last.(polygon_points)
# xrange = LinRange(0, 0.5, 1000)
# yrange = LinRange(0, 0.5, 1000)
# @benchmark mean_value_coordinate_field = mean_value_barycentric_coordinates.(
#     (GeometryBasics.Polygon(Point2f.(polygon_points)),), 
#     Point2f.(xrange, yrange')
# )
# # @time mean_value_coordinate_field_threadsx = ThreadsX.map(Point2{Float64}.(xrange, yrange')) do p
# #     mean_value_barycentric_coordinates(polygon_points, p)
# # end # 7s v/s 2s for serial processing...why?
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
# f, a, p = heatmap(xrange, yrange, mean_values; colormap = cgrad(:jet, 18; categorical = true), axis = (; aspect = DataAspect()))
# poly!(Point2f.(polygon_points); color = :transparent, strokecolor = :black, strokewidth = 1.3)
# Colorbar(f[1, 2], p)
# f
# ```


"""
    abstract type AbstractBarycentricCoordinateMethod

Abstract supertype for barycentric coordinate methods.  
The subtypes may serve as dispatch types, or may cache 
some information about the target polygon.  

## API
The following methods must be implemented for all subtypes:
- `barycentric_coordinates!(λs::Vector{<: Real}, method::AbstractBarycentricCoordinateMethod, polypoints::Vector{<: Point{N1, T1}}, point::Point{N2, T2})`
"""
abstract type AbstractBarycentricCoordinateMethod end

Base.@propagate_inbounds function barycentric_coordinates!(λs::Vector{<: Real}, method::AbstractBarycentricCoordinateMethod, polypoints::Vector{<: Point{N1, T1}}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real}
    @boundscheck @assert length(λs) == length(polypoints)
    @boundscheck @assert length(polypoints) >= 3

    @error("Not implemented yet for method $(method).")
end

Base.@propagate_inbounds function barycentric_coordinates(method::AbstractBarycentricCoordinateMethod, polypoints::Vector{<: Point{N1, T1}}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real}
    λs = zeros(promote_type(T1, T2), length(polypoints))
    barycentric_coordinates!(λs, method, polypoints, point)
    return λs
end

struct MeanValue <: AbstractBarycentricCoordinateMethod 
end

struct Wachspress <: AbstractBarycentricCoordinateMethod
end

# ```@example barycentric
# n = 200
# angles = range(0, 2pi, length = n)
# x_ext = 2 .* cos.(angles .+ pi/n)
# y_ext = 2 .* sin.(angles .+ pi/n)
# x_int = cos.(angles)
# y_int = sin.(angles)
# z_ext = (x_ext .- 0.5).^2 + (y_ext .- 0.5).^2 .+ 0.05.*randn.()
# z_int = (x_int .- 0.5).^2 + (y_int .- 0.5).^2 .+ 0.05.*randn.()
# z = vcat(z_ext, z_int)
# circ_poly = Polygon(Point2f.(x_ext, y_ext), [Point2f.(x_int, y_int)])
# ext = GeoInterface.extent(circ_poly)
# xrange = LinRange(ext.X..., 1000)
# yrange = LinRange(ext.Y..., 1000)
# @time mean_value_coordinate_field = mean_value_barycentric_coordinates.(
#     (circ_poly,),
#     Point2f.(xrange, yrange')
# )

# mean_values = ThreadsX.map(mean_value_coordinate_field) do λs
#     sum(λs .* z)
# end
# f, a, p = heatmap(xrange, yrange, mean_values; colormap = :viridis, axis = (; aspect = DataAspect()))
# poly!(a, circ_poly; color = :transparent, strokecolor = :black, strokewidth = 1.3)
# f
# ```
