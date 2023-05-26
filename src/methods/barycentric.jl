# # Barycentric coordinates
# Generalized barycentric coordinates are a generalization of barycentric coordinates, 
# which are typically used in triangles, to arbitrary polygons. 

# They provide a way to express a point within a polygon as a weighted average 
# of the polygon's vertices.

# In the case of a triangle, barycentric coordinates are a set of three numbers 
# $(λ_1, λ_2, λ_3)$, each associated with a vertex of the triangle. Any point within 
# the triangle can be expressed as a weighted average of the vertices, where the 
# weights are the barycentric coordinates. The weights sum to 1, and each is non-negative.

# For a polygon with $n$ vertices, generalized barycentric coordinates are a set of 
# $n$ numbers $(λ_1, λ_2, ..., λ_n)$, each associated with a vertex of the polygon. 
# Any point within the polygon can be expressed as a weighted average of the vertices, 
# where the weights are the generalized barycentric coordinates. 

# As with the triangle case, the weights sum to 1, and each is non-negative.

# ## Barycentric-coordinate API
# In most cases, we actually want barycentric interpolation and have no interest in the coordinates themselves.  However, the coordinates can be useful for debugging, and so we provide an API for computing them as well.
#


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


Base.@propagate_inbounds function barycentric_coordinates!(λs::Vector{<: Real}, method::AbstractBarycentricCoordinateMethod, polypoints::AbstractVector{<: Point{N1, T1}}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real}
    @boundscheck @assert length(λs) == length(polypoints)
    @boundscheck @assert length(polypoints) >= 3

    @error("Not implemented yet for method $(method).")
end

Base.@propagate_inbounds function barycentric_coordinates(method::AbstractBarycentricCoordinateMethod, polypoints::AbstractVector{<: Point{N1, T1}}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real}
    λs = zeros(promote_type(T1, T2), length(polypoints))
    barycentric_coordinates!(λs, method, polypoints, point)
    return λs
end

Base.@propagate_inbounds function barycentric_interpolate(method::AbstractBarycentricCoordinateMethod, polypoints::AbstractVector{<: Point{N, T1}}, values::AbstractVector{V}, point::Point{N, T2}) where {N, T1 <: Real, T2 <: Real, V}
    @boundscheck @assert length(values) == length(polypoints)
    @boundscheck @assert length(polypoints) >= 3
    λs = barycentric_coordinates(method, polypoints, point)
    return sum(λs .* values)
end



struct MeanValue <: AbstractBarycentricCoordinateMethod 
end

struct Wachspress <: AbstractBarycentricCoordinateMethod
end
# function mean_value_barycentric_coordinates(::PolygonTrait, ::PointTrait, )

# ## Example
# This example was taken from [this page of CGAL's documentation](https://doc.cgal.org/latest/Barycentric_coordinates_2/index.html).

# ```@example barycentric
using GeometryBasics, GeometryOps, Makie
polygon_points = Point3f[
(0.03, 0.05, 0.00), (0.07, 0.04, 0.02), (0.10, 0.04, 0.04),
(0.14, 0.04, 0.06), (0.17, 0.07, 0.08), (0.20, 0.09, 0.10),
(0.22, 0.11, 0.12), (0.25, 0.11, 0.14), (0.27, 0.10, 0.16),
(0.30, 0.07, 0.18), (0.31, 0.04, 0.20), (0.34, 0.03, 0.22),
(0.37, 0.02, 0.24), (0.40, 0.03, 0.26), (0.42, 0.04, 0.28),
(0.44, 0.07, 0.30), (0.45, 0.10, 0.32), (0.46, 0.13, 0.34),
(0.46, 0.19, 0.36), (0.47, 0.26, 0.38), (0.47, 0.31, 0.40),
(0.47, 0.35, 0.42), (0.45, 0.37, 0.44), (0.41, 0.38, 0.46),
(0.38, 0.37, 0.48), (0.35, 0.36, 0.50), (0.32, 0.35, 0.52),
(0.30, 0.37, 0.54), (0.28, 0.39, 0.56), (0.25, 0.40, 0.58),
(0.23, 0.39, 0.60), (0.21, 0.37, 0.62), (0.21, 0.34, 0.64),
(0.23, 0.32, 0.66), (0.24, 0.29, 0.68), (0.27, 0.24, 0.70),
(0.29, 0.21, 0.72), (0.29, 0.18, 0.74), (0.26, 0.16, 0.76),
(0.24, 0.17, 0.78), (0.23, 0.19, 0.80), (0.24, 0.22, 0.82),
(0.24, 0.25, 0.84), (0.21, 0.26, 0.86), (0.17, 0.26, 0.88),
(0.12, 0.24, 0.90), (0.07, 0.20, 0.92), (0.03, 0.15, 0.94),
(0.01, 0.10, 0.97), (0.02, 0.07, 1.00)]
using MakieThemes
Makie.set_theme!(MakieThemes.bbc())
f, a, p = poly(polygon_points; color = last.(polygon_points), colormap = cgrad(:jet, 18; categorical = true), shading = false, axis = (; aspect = DataAspect(), title = "Makie mesh based polygon rendering", subtitle = "Makie"))
cb = Colorbar(f[1, 2], p.plots[1])
hidedecorations!(a)
f
ax_bbox = a.finallimits[]
ext = GeoInterface.Extent(NamedTuple{(:X, :Y)}(zip(minimum(ax_bbox), maximum(ax_bbox))))
poly_rast = Rasters.rasterize(GeometryBasics.Polygon(Point2f.(polygon_points)); ext = ext, size = tuple(round.(Int, widths(a.scene.px_area[]))...), fill = RGBAf(0,0,0,0))
@time mean_value_coordinate_field = hormann_mean_value_coordinates.(
    (Point2f.(polygon_points),), 
    Point2f.(collect(poly_rast.dims[1]), collect(poly_rast.dims[2])')
)
zs = last.(polygon_points)
mean_values = map(mean_value_coordinate_field) do λs
    sum(λs .* zs)
end

@time mean_values = hormann_mean_value_interpolation.(
    (Point2f.(polygon_points),),
    (zs,),
    Point2f.(collect(poly_rast.dims[1]), collect(poly_rast.dims[2])')
)
fig, ax, mvplt = heatmap(collect(poly_rast.dims[1]), collect(poly_rast.dims[2]), mean_values; colormap = cgrad(:jet, 18; categorical = true), axis = (; aspect = DataAspect(), title = "Barycentric coordinate based rendering", subtitle = "Mean value method"), colorrange = Makie.distinct_extrema_nan(zs))
hidedecorations!(ax)
cb = Colorbar(fig[1, 2], mvplt)
poly!(ax, GeometryBasics.Polygon(Point2f[(ext.X[1], ext.Y[1]), (ext.X[2], ext.Y[1]), (ext.X[2], ext.Y[2]), (ext.X[1], ext.Y[2]), (ext.X[1], ext.Y[1])], [reverse(Point2f.(polygon_points))]); color = :white, xautolimits = false, yautolimits = false)
fig
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



# ```@example barycentric
n = 200
angles = range(0, 2pi, length = n)[1:end-1]
x_ext = 2 .* cos.(angles .+ pi/n)
y_ext = 2 .* sin.(angles .+ pi/n)
x_int = cos.(angles)
y_int = sin.(angles)
z_ext = (x_ext .- 0.5).^2 + (y_ext .- 0.5).^2 .+ 0.05.*randn.()
z_int = (x_int .- 0.5).^2 + (y_int .- 0.5).^2 .+ 0.05.*randn.()
z = vcat(z_ext, z_int)
circ_poly = Polygon(Point2f.(x_ext, y_ext), [Point2f.(x_int, y_int)])
circ_poly = Polygon(Point{3, Float64}.(polygon_points))
ext = GeoInterface.extent(circ_poly)
xrange = LinRange(ext.X..., 1000) .|> Float64
yrange = LinRange(ext.Y..., 1000) .|> Float64
@time mean_value_coordinate_field = hormann_mean_value_interpolation.(
    (circ_poly,),
    (z,),
    Point2f.(xrange, yrange')
)
itp_xs = [xrange[i] for i in 1:length(xrange), j in 1:length(yrange)]
itp_ys = [yrange[j] for i in 1:length(xrange), j in 1:length(yrange)]
vals = itp(vec(itp_xs), vec(itp_ys); method = NaturalNeighbours.Sibson(1))
heatmap(xrange, yrange, reshape(vals, (length(xrange), length(yrange))); colormap = cgrad(:jet, 18; categorical = true), axis = (; aspect = DataAspect(), title = "NaturalNeighbours.jl Laplace interpolation"))
poly!(Point2f.(polygon_points); color = :transparent, strokecolor = :black, strokewidth = 1.3)
Colorbar(Makie.current_figure()[1, 2], Makie.current_axis().scene.plots[1])
Makie.current_figure()
itp = NaturalNeighbours.interpolate(vcat(x_ext, x_int), vcat(y_ext, y_int), z; derivatives = true)
itp_xs = [xrange[i] for i in 1:length(xrange), j in 1:length(yrange)]
itp_ys = [yrange[j] for i in 1:length(xrange), j in 1:length(yrange)]
vals = itp(vec(itp_xs), vec(itp_ys); method = NaturalNeighbours.Laplace())
heatmap(xrange, yrange, reshape(vals, (length(xrange), length(yrange))); axis = (; aspect = DataAspect(), title = "NaturalNeighbours.jl Laplace interpolation"))
mean_values = ThreadsX.map(mean_value_coordinate_field) do λs
    sum(λs .* z)
end

function NaturalNeighbours.interpolate(points::AbstractVector{<: Point{2, T1}}, zs::AbstractVector{T2}; kwargs...) where {T1 <: Real, T2 <: Real}
    return NaturalNeighbours.interpolate(getindex.(points, 1), getindex.(points, 2), zs; kwargs...)
end

function NaturalNeighbours.interpolate(points::AbstractVector{<: Point{3, T}}; kwargs...) where {T <: Real}
    return NaturalNeighbours.interpolate(getindex.(points, 1), getindex.(points, 2), getindex.(points, 3); kwargs...)
end

function NaturalNeighbours.interpolate(poly::Polygon{3, T}; kwargs...) where {T <: Real}
    # Decompose the exterior and interiors of the polygon into a list of points
    poly_ext_points = decompose(Point{3, T}, poly.exterior)
    poly_int_points = decompose.((Point{3, T},), poly.interiors)
    final_point_list = Point{3, T}[]
    # If the exterior or interior points are closed, we need to remove the last point,
    # to avoid degeneracy and duplicate points.
    if poly_ext_points[end] == poly_ext_points[begin]
        append!(final_point_list, @view(poly_ext_points[begin:end-1]))
    else
        append!(final_point_list, poly_ext_points)
    end
    for hole in poly_int_points
        if hole[end] == hole[begin]
            append!(final_point_list, @view(hole[begin:end-1]))
        else
            append!(final_point_list, hole)
        end
    end
    return NaturalNeighbours.interpolate(final_point_list; kwargs...)
end

# f, a, p = heatmap(xrange, yrange, mean_values; colormap = :viridis, axis = (; aspect = DataAspect()))
# poly!(a, circ_poly; color = :transparent, strokecolor = :black, strokewidth = 1.3)
# f
# ```

# ```
# x = randn(50)
# y = randn(50)
# z = -sqrt.(x .^ 2 .+ y .^ 2) .+ 0.1 .* randn.()
#
# xrange = LinRange(ext.X..., 1000)
# yrange = LinRange(ext.Y..., 1000)
# @time mean_value_coordinate_field = mean_value_barycentric_coordinates.(
#     (Point2f.(x, y),),
#     Point2f.(xrange, yrange')
# )
#
# mean_values = ThreadsX.map(mean_value_coordinate_field) do λs
#     sum(λs .* z)
# end
#
# contour(getindex.(λs, 2))
# heatmap(mean_values)
# ```

xs = randn(100)
ys = randn(100)
zs = sin.(xs) .* cos.(ys)

itp = NaturalNeighbours.interpolate(xs, ys, zs; derivatives = true)
xrange = LinRange(-3, 3, 1000) .|> Float32
yrange = LinRange(-3, 3, 1000) .|> Float32
itp_xs = [xrange[i] for i in 1:length(xrange), j in 1:length(yrange)]
itp_ys = [yrange[j] for i in 1:length(xrange), j in 1:length(yrange)]
vals = itp(vec(itp_xs), vec(itp_ys); method = NaturalNeighbours.Sibson(1))
heatmap(xrange, yrange, reshape(vals, (length(xrange), length(yrange))))
scatter!(xs, ys)
Makie.current_figure()

"""
    _det(s1::Point2{T1}, s2::Point2{T2}) where {T1 <: Real, T2 <: Real}

Returns the determinant of the matrix formed by the two points `s1` and `s2`.
Specifically, this is: 
```julia
s1[1] * s2[2] - s1[2] * s2[1]
```

## Extended help

## Doctests

```jldoctest
julia> _det((1,0), (0,1))
1

julia> _det(Point2f(1, 2), Point2f(3, 4))
-2.0f0
```
"""
function _det(s1::Point2{T1}, s2::Point2{T2}) where {T1 <: Real, T2 <: Real}
    return s1[1] * s2[2] - s1[2] * s2[1]
end

"""
    t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)

Returns the "T-value" as described in Hormann's presentation on how to calculate
the mean-value coordinate.

See [`hormann_mean_value_interpolation`](@ref) for more details.
"""
function t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)
    return _det(sᵢ, sᵢ₊₁) / muladd(rᵢ, rᵢ₊₁, dot(sᵢ, sᵢ₊₁))
end

"""
    hormann_mean_value_coordinates(points::AbstractVector{<: Point2}, point::Point2)::Vector

Returns a vector of the length of `points`, containing the mean value coordinates of `point` with respect to `points`.

# Arguments

- `points::AbstractVector{<: Point{2, T1}}`: A vector of points defining a polygon.  The end point should not repeat.
- `point::Point{2, T2}`: The point to compute the mean value coordinates for.

"""
function hormann_mean_value_coordinates(points::AbstractVector{<: Point{2, T1}}, point::Point{2, T2}) where {T1 <: Real, T2 <: Real}
    # First, allocate the output array.
    # 
    n_points = length(points)
    NumType = promote_type(T1, T2)
    λs = zeros(NumType, n_points)
    sᵢ₋₁ = points[end] - point
    sᵢ   = points[begin] - point
    sᵢ₊₁ = points[begin+1] - point
    rᵢ₋₁ = norm(sᵢ₋₁) # radius / Euclidean distance between points.
    rᵢ   = norm(sᵢ  ) # radius / Euclidean distance between points.
    rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.

    λs[1] = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 

    for i in 2:n_points
        # Increment counters + set variables
        sᵢ₋₁ = sᵢ
        sᵢ   = sᵢ₊₁
        sᵢ₊₁ = points[mod1(i+1, n_points)] - point
        rᵢ₋₁ = rᵢ
        rᵢ   = rᵢ₊₁
        rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
        λs[i] = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ
    end
    λs ./= sum(λs) 
    return λs
end

function hormann_mean_value_interpolation(points::AbstractVector{<: Point{2, T1}}, values::AbstractVector{T2}, point::Point{2, T3}) where {T1 <: Real, T2, T3 <: Real}
    n_points = length(points)
    sᵢ₋₁ = points[end] - point
    sᵢ   = points[begin] - point
    sᵢ₊₁ = points[begin+1] - point
    rᵢ₋₁ = norm(sᵢ₋₁) # radius / Euclidean distance between points.
    rᵢ   = norm(sᵢ  ) # radius / Euclidean distance between points.
    rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
    # Now, we set the interpolated value to the first point's value, multiplied
    # by the weight computed relative to the first point in the polygon. 
    wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
    wₜₒₜ = wᵢ
    interpolated_value = values[begin] * wᵢ
    for i in 2:n_points
        # Increment counters + set variables
        sᵢ₋₁ = sᵢ
        sᵢ   = sᵢ₊₁
        sᵢ₊₁ = points[mod1(i+1, n_points)] - point
        rᵢ₋₁ = rᵢ
        rᵢ   = rᵢ₊₁
        rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
        wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
        wₜₒₜ += wᵢ
        interpolated_value += values[i] * wᵢ
    end
    interpolated_value /= wₜₒₜ
    return interpolated_value
end

function hormann_mean_value_interpolation(poly::GeometryBasics.Polygon{2, T1}, values::AbstractVector{T2}, point::Point{2, T3}) where {T1 <: Real, T2, T3 <: Real}
    ext = decompose(Point{2, T2}, poly.exterior)
    ints = decompose.((Point{2, T2},), poly.interiors)
    current_index = 1
    l_ext = length(ext)

    sᵢ₋₁ = ext[end] - point
    sᵢ   = ext[begin] - point
    sᵢ₊₁ = ext[begin+1] - point
    rᵢ₋₁ = norm(sᵢ₋₁) # radius / Euclidean distance between points.
    rᵢ   = norm(sᵢ  ) # radius / Euclidean distance between points.
    rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
    # Now, we set the interpolated value to the first point's value, multiplied
    # by the weight computed relative to the first point in the polygon. 
    wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
    wₜₒₜ = wᵢ
    interpolated_value = values[begin] * wᵢ

    for i in 2:l_ext
        # Increment counters + set variables
        sᵢ₋₁ = sᵢ
        sᵢ   = sᵢ₊₁
        sᵢ₊₁ = ext[mod1(i+1, l_ext)] - point
        rᵢ₋₁ = rᵢ
        rᵢ   = rᵢ₊₁
        rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
        wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
        # Updates - first the interpolated value,
        interpolated_value += values[current_index] * wᵢ
        # then the accumulators for total weight and current index.
        wₜₒₜ += wᵢ
        current_index += 1

    end
    for hole in ints
        l_hole = length(hole)
        sᵢ₋₁ = hole[end] - point
        sᵢ   = hole[begin] - point
        sᵢ₊₁ = hole[begin+1] - point
        rᵢ₋₁ = norm(sᵢ₋₁) # radius / Euclidean distance between points.
        rᵢ   = norm(sᵢ  ) # radius / Euclidean distance between points.
        rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
        # Now, we set the interpolated value to the first point's value, multiplied
        # by the weight computed relative to the first point in the polygon. 
        wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ

        interpolated_value += values[current_index] * wᵢ

        wₜₒₜ += wᵢ
        current_index += 1
    
        for i in 2:l_hole
            # Increment counters + set variables
            sᵢ₋₁ = sᵢ
            sᵢ   = sᵢ₊₁
            sᵢ₊₁ = hole[mod1(i+1, l_hole)] - point
            rᵢ₋₁ = rᵢ
            rᵢ   = rᵢ₊₁
            rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
            wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
            interpolated_value += values[current_index] * wᵢ
            wₜₒₜ += wᵢ
            current_index += 1
        end
    end
    return interpolated_value / wₜₒₜ

end