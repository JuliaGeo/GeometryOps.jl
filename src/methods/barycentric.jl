# # Barycentric coordinates

export barycentric_coordinates, barycentric_coordinates!, barycentric_interpolate
export MeanValue

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


# ## Example
# This example was taken from [this page of CGAL's documentation](https://doc.cgal.org/latest/Barycentric_coordinates_2/index.html).
#=
```@example barycentric
import GeometryOps as GO
using GeometryOps.GeometryBasics
using Makie
using CairoMakie
# Define a polygon
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
# Plot it!
# First, we'll plot the polygon using Makie's rendering:
f, a1, p1 = poly(
    polygon_points; 
    color = last.(polygon_points), colormap = cgrad(:jet, 18; categorical = true), 
    axis = (; 
        aspect = DataAspect(), title = "Makie mesh based polygon rendering", subtitle = "CairoMakie"
    ), 
    figure = (; resolution = (800, 400),)
)

Makie.update_state_before_display!(f) # We have to call this explicitly, to get the axis limits correct
# Now that we've plotted the first polygon,
# we can render it using barycentric coordinates.
a1_bbox = a1.finallimits[] # First we get the extent of the axis
ext = GeometryOps.GI.Extent(NamedTuple{(:X, :Y)}(zip(minimum(a1_bbox), maximum(a1_bbox))))

a2 = Axis(
        f[1, 2], 
        aspect = DataAspect(), 
        title = "Barycentric coordinate based polygon rendering", subtitle = "GeometryOps",
        limits = (ext.X, ext.Y)
    )
p2box = poly!( # Now, we plot a cropping rectangle around the axis so we only show the polygon
    a2, 
    GeometryOps.GeometryBasics.Polygon( # This is a rectangle with an internal hole shaped like the polygon.
        Point2f[(ext.X[1], ext.Y[1]), (ext.X[2], ext.Y[1]), (ext.X[2], ext.Y[2]), (ext.X[1], ext.Y[2]), (ext.X[1], ext.Y[1])], 
        [reverse(Point2f.(polygon_points))]
    ); 
    color = :white, xautolimits = false, yautolimits = false
)
hidedecorations!(a1)
hidedecorations!(a2)
cb = Colorbar(f[2, :], p1.plots[1]; vertical = false, flipaxis = true)
# Finally, we perform barycentric interpolation on a grid,
xrange = LinRange(ext.X..., widths(a2.scene.px_area[])[1] * 4) # 2 rendered pixels per "physical" pixel
yrange = LinRange(ext.Y..., widths(a2.scene.px_area[])[2] * 4) # 2 rendered pixels per "physical" pixel
@time mean_values = barycentric_interpolate.(
    (MeanValue(),), # The barycentric coordinate algorithm (MeanValue is the only one for now)
    (Point2f.(polygon_points),), # The polygon points as `Point2f`
    (last.(polygon_points,),),   # The values per polygon point - can be anything which supports addition and division
    Point2f.(xrange, yrange')    # The points at which to interpolate
)
# and render!
hm = heatmap!(
    a2, xrange, yrange, mean_values;
    colormap = p1.colormap, # Use the same colormap as the original polygon plot
    colorrange = p1.plots[1].colorrange[], # Access the rendered mesh plot's colorrange directly
    transformation = (; translation = Vec3f(0,0,-1)), # This gets the heatmap to render "behind" the previously plotted polygon
    xautolimits = false, yautolimits = false
)
f
```

## Barycentric-coordinate API
In some cases, we actually want barycentric interpolation, and have no interest 
in the coordinates themselves.  

However, the coordinates can be useful for debugging, and when performing 3D rendering,
multiple barycentric values (depth, uv) are needed for depth buffering.
=#

const _VecTypes = Union{Tuple{Vararg{T, N}}, GeometryBasics.StaticArraysCore.StaticArray{Tuple{N}, T, 1}} where {N, T}

"""
    abstract type AbstractBarycentricCoordinateMethod

Abstract supertype for barycentric coordinate methods.  
The subtypes may serve as dispatch types, or may cache 
some information about the target polygon.  

## API
The following methods must be implemented for all subtypes:
- `barycentric_coordinates!(λs::Vector{<: Real}, method::AbstractBarycentricCoordinateMethod, exterior::Vector{<: Point{2, T1}}, point::Point{2, T2})`
- `barycentric_interpolate(method::AbstractBarycentricCoordinateMethod, exterior::Vector{<: Point{2, T1}}, values::Vector{V}, point::Point{2, T2})::V`
- `barycentric_interpolate(method::AbstractBarycentricCoordinateMethod, exterior::Vector{<: Point{2, T1}}, interiors::Vector{<: Vector{<: Point{2, T1}}} values::Vector{V}, point::Point{2, T2})::V`
The rest of the methods will be implemented in terms of these, and have efficient dispatches for broadcasting.
"""
abstract type AbstractBarycentricCoordinateMethod end


Base.@propagate_inbounds function barycentric_coordinates!(λs::Vector{<: Real}, method::AbstractBarycentricCoordinateMethod, polypoints::AbstractVector{<: Point{N1, T1}}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real}
    @boundscheck @assert length(λs) == length(polypoints)
    @boundscheck @assert length(polypoints) >= 3

    @error("Not implemented yet for method $(method).")
end
Base.@propagate_inbounds barycentric_coordinates!(λs::Vector{<: Real}, polypoints::AbstractVector{<: Point{N1, T1}}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real} = barycentric_coordinates!(λs, MeanValue(), polypoints, point)

Base.@propagate_inbounds function barycentric_coordinates(method::AbstractBarycentricCoordinateMethod, polypoints::AbstractVector{<: Point{N1, T1}}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real}
    λs = zeros(promote_type(T1, T2), length(polypoints))
    barycentric_coordinates!(λs, method, polypoints, point)
    return λs
end
Base.@propagate_inbounds barycentric_coordinates(polypoints::AbstractVector{<: Point{N1, T1}}, point::Point{N2, T2}) where {N1, N2, T1 <: Real, T2 <: Real} = barycentric_coordinates(MeanValue(), polypoints, point)

Base.@propagate_inbounds function barycentric_interpolate(method::AbstractBarycentricCoordinateMethod, polypoints::AbstractVector{<: Point{N, T1}}, values::AbstractVector{V}, point::Point{N, T2}) where {N, T1 <: Real, T2 <: Real, V}
    @boundscheck @assert length(values) == length(polypoints)
    @boundscheck @assert length(polypoints) >= 3
    λs = barycentric_coordinates(method, polypoints, point)
    return sum(λs .* values)
end
Base.@propagate_inbounds barycentric_interpolate(polypoints::AbstractVector{<: Point{N, T1}}, values::AbstractVector{V}, point::Point{N, T2}) where {N, T1 <: Real, T2 <: Real, V} = barycentric_interpolate(MeanValue(), polypoints, values, point)

Base.@propagate_inbounds function barycentric_interpolate(method::AbstractBarycentricCoordinateMethod, exterior::AbstractVector{<: Point{N, T1}}, interiors::AbstractVector{<: Point{N, T1}}, values::AbstractVector{V}, point::Point{N, T2}) where {N, T1 <: Real, T2 <: Real, V}
    @boundscheck @assert length(values) == length(exterior) + isempty(interiors) ? 0 : sum(length.(interiors))
    @boundscheck @assert length(exterior) >= 3
    λs = barycentric_coordinates(method, exterior, interiors, point)
    return sum(λs .* values)
end
Base.@propagate_inbounds barycentric_interpolate(exterior::AbstractVector{<: Point{N, T1}}, interiors::AbstractVector{<: Point{N, T1}}, values::AbstractVector{V}, point::Point{N, T2}) where {N, T1 <: Real, T2 <: Real, V} = barycentric_interpolate(MeanValue(), exterior, interiors, values, point)

Base.@propagate_inbounds function barycentric_interpolate(method::AbstractBarycentricCoordinateMethod, polygon::Polygon{2, T1}, values::AbstractVector{V}, point::Point{2, T2}) where {T1 <: Real, T2 <: Real, V}
    exterior = decompose(Point{2, promote_type(T1, T2)}, polygon.exterior)
    if isempty(polygon.interiors)
        @boundscheck @assert length(values) == length(exterior)
        return barycentric_interpolate(method, exterior, values, point)
    else # the poly has interiors
        interiors = reverse.(decompose.((Point{2, promote_type(T1, T2)},), polygon.interiors))
        @boundscheck @assert length(values) == length(exterior) + sum(length.(interiors))
        return barycentric_interpolate(method, exterior, interiors, values, point)
    end
end
Base.@propagate_inbounds barycentric_interpolate(polygon::Polygon{2, T1}, values::AbstractVector{V}, point::Point{2, T2}) where {T1 <: Real, T2 <: Real, V} = barycentric_interpolate(MeanValue(), polygon, values, point)

# 3D polygons are considered to have their vertices in the XY plane, 
# and the Z coordinate must represent some value.  This is to say that
# the Z coordinate is interpreted as an M coordinate.
Base.@propagate_inbounds function barycentric_interpolate(method::AbstractBarycentricCoordinateMethod, polygon::Polygon{3, T1}, point::Point{2, T2}) where {T1 <: Real, T2 <: Real}
    exterior_point3s = decompose(Point{3, promote_type(T1, T2)}, polygon.exterior)
    exterior_values = getindex.(exterior_point3s, 3)
    exterior_points = Point2f.(exterior_point3s)
    if isempty(polygon.interiors)
        return barycentric_interpolate(method, exterior_points, exterior_values, point)
    else # the poly has interiors
        interior_point3s = decompose.((Point{3, promote_type(T1, T2)},), polygon.interiors)
        interior_values = collect(Iterators.flatten((getindex.(point3s, 3) for point3s in interior_point3s)))
        interior_points = map(point3s -> Point2f.(point3s), interior_point3s)
        return barycentric_interpolate(method, exterior_points, interior_points, vcat(exterior_values, interior_values), point)
    end
end
Base.@propagate_inbounds barycentric_interpolate(polygon::Polygon{3, T1}, point::Point{2, T2}) where {T1 <: Real, T2 <: Real} = barycentric_interpolate(MeanValue(), polygon, point)

# This method is the one which supports GeoInterface.
Base.@propagate_inbounds function barycentric_interpolate(method::AbstractBarycentricCoordinateMethod, polygon, values::AbstractVector{V}, point) where V
    @assert GeoInterface.trait(polygon) isa GeoInterface.PolygonTrait
    @assert GeoInterface.trait(point) isa GeoInterface.PointTrait
    passable_polygon = GeoInterface.convert(GeometryBasics, polygon)
    @assert passable_polygon isa GeometryBasics.Polygon "The polygon was converted to a $(typeof(passable_polygon)), which is not a `GeometryBasics.Polygon`."
    ## first_poly_point = GeoInterface.getpoint(GeoInterface.getexterior(polygon))
    passable_point = GeoInterface.convert(GeometryBasics, point)
    return barycentric_interpolate(method, passable_polygon, Point2(passable_point))
end
Base.@propagate_inbounds barycentric_interpolate(polygon, values::AbstractVector{V}, point) where V = barycentric_interpolate(MeanValue(), polygon, values, point)

"""
    weighted_mean(weight::Real, x1, x2)

Returns the weighted mean of `x1` and `x2`, where `weight` is the weight of `x1`.

Specifically, calculates `x1 * weight + x2 * (1 - weight)`.

!!! note
    The idea for this method is that you can override this for custom types, like Color types, in extension modules.
"""
function weighted_mean(weight::WT, x1, x2) where {WT <: Real}
    return muladd(x1, weight, x2 * (oneunit(WT) - weight))
end


"""
    MeanValue() <: AbstractBarycentricCoordinateMethod

This method calculates barycentric coordinates using the mean value method.

## References

"""
struct MeanValue <: AbstractBarycentricCoordinateMethod 
end

# Before we go to the actual implementation, there are some quick and simple utility functions
# that we need to implement.  These are mainly for convenience and code brevity.

"""
    _det(s1::Point2{T1}, s2::Point2{T2}) where {T1 <: Real, T2 <: Real}

Returns the determinant of the matrix formed by `hcat`'ing two points `s1` and `s2`.

Specifically, this is: 
```julia
s1[1] * s2[2] - s1[2] * s2[1]
```
"""
function _det(s1::_VecTypes{2, T1}, s2::_VecTypes{2, T2}) where {T1 <: Real, T2 <: Real}
    return s1[1] * s2[2] - s1[2] * s2[1]
end

"""
    t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)

Returns the "T-value" as described in Hormann's presentation [^HormannPresentation] on how to calculate
the mean-value coordinate.  

Here, `sᵢ` is the vector from vertex `vᵢ` to the point, and `rᵢ` is the norm (length) of `sᵢ`.
`s` must be `Point` and `r` must be real numbers.

```math
tᵢ = \\frac{\\mathrm{det}\\left(sᵢ, sᵢ₊₁\\right)}{rᵢ * rᵢ₊₁ + sᵢ ⋅ sᵢ₊₁}
```

[^HormannPresentation]: K. Hormann and N. Sukumar. Generalized Barycentric Coordinates in Computer Graphics and Computational Mechanics. Taylor & Fancis, CRC Press, 2017.
```

"""
function t_value(sᵢ::_VecTypes{N, T1}, sᵢ₊₁::_VecTypes{N, T1}, rᵢ::T2, rᵢ₊₁::T2) where {N, T1 <: Real, T2 <: Real}
    return _det(sᵢ, sᵢ₊₁) / muladd(rᵢ, rᵢ₊₁, dot(sᵢ, sᵢ₊₁))
end


function barycentric_coordinates!(λs::Vector{<: Real}, ::MeanValue, polypoints::AbstractVector{<: Point{2, T1}}, point::Point{2, T2}) where {T1 <: Real, T2 <: Real}
    @boundscheck @assert length(λs) == length(polypoints)
    @boundscheck @assert length(polypoints) >= 3
    n_points = length(polypoints)
    ## Initialize counters and register variables
    ## Points - these are actually vectors from point to vertices
    ##  polypoints[i-1], polypoints[i], polypoints[i+1]
    sᵢ₋₁ = polypoints[end] - point
    sᵢ   = polypoints[begin] - point
    sᵢ₊₁ = polypoints[begin+1] - point
    ## radius / Euclidean distance between points.
    rᵢ₋₁ = norm(sᵢ₋₁) 
    rᵢ   = norm(sᵢ  )
    rᵢ₊₁ = norm(sᵢ₊₁)
    ## Perform the first computation explicitly, so we can cut down on 
    ## a mod in the loop.
    λs[1] = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
    ## Loop through the rest of the vertices, compute, store in λs
    for i in 2:n_points
        ## Increment counters + set variables
        sᵢ₋₁ = sᵢ
        sᵢ   = sᵢ₊₁
        sᵢ₊₁ = polypoints[mod1(i+1, n_points)] - point
        rᵢ₋₁ = rᵢ
        rᵢ   = rᵢ₊₁
        rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
        λs[i] = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ
    end
    ## Normalize λs to the 1-norm (sum=1)
    λs ./= sum(λs) 
    return λs
end

# ```julia
# function barycentric_coordinates(::MeanValue, polypoints::NTuple{N, Point{2, T2}}, point::Point{2, T1},) where {N, T1, T2}
#     ## Initialize counters and register variables
#     ## Points - these are actually vectors from point to vertices
#     ##  polypoints[i-1], polypoints[i], polypoints[i+1]
#     sᵢ₋₁ = polypoints[end] - point
#     sᵢ   = polypoints[begin] - point
#     sᵢ₊₁ = polypoints[begin+1] - point
#     ## radius / Euclidean distance between points.
#     rᵢ₋₁ = norm(sᵢ₋₁) 
#     rᵢ   = norm(sᵢ  )
#     rᵢ₊₁ = norm(sᵢ₊₁)
#     λ₁ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
#     λs = ntuple(N) do i
#         if i == 1
#             return λ₁
#         end
#         ## Increment counters + set variables
#         sᵢ₋₁ = sᵢ
#         sᵢ   = sᵢ₊₁
#         sᵢ₊₁ = polypoints[mod1(i+1, N)] - point
#         rᵢ₋₁ = rᵢ
#         rᵢ   = rᵢ₊₁
#         rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
#         return (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ
#     end
#
#     ∑λ = sum(λs)
#
#     return ntuple(N) do i
#         λs[i] / ∑λ
#     end
# end
# ```

# This performs an inplace accumulation, using less memory and is faster.
# That's particularly good if you are using a polygon with a large number of points...
function barycentric_interpolate(::MeanValue, polypoints::AbstractVector{<: Point{2, T1}}, values::AbstractVector{V}, point::Point{2, T2}) where {T1 <: Real, T2 <: Real, V}
    @boundscheck @assert length(values) == length(polypoints)
    @boundscheck @assert length(polypoints) >= 3
    
    n_points = length(polypoints)
    ## Initialize counters and register variables
    ## Points - these are actually vectors from point to vertices
    ##  polypoints[i-1], polypoints[i], polypoints[i+1]
    sᵢ₋₁ = polypoints[end] - point
    sᵢ   = polypoints[begin] - point
    sᵢ₊₁ = polypoints[begin+1] - point
    ## radius / Euclidean distance between points.
    rᵢ₋₁ = norm(sᵢ₋₁) 
    rᵢ   = norm(sᵢ  )
    rᵢ₊₁ = norm(sᵢ₊₁)
    ## Now, we set the interpolated value to the first point's value, multiplied
    ## by the weight computed relative to the first point in the polygon. 
    wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
    wₜₒₜ = wᵢ
    interpolated_value = values[begin] * wᵢ
    for i in 2:n_points
        ## Increment counters + set variables
        sᵢ₋₁ = sᵢ
        sᵢ   = sᵢ₊₁
        sᵢ₊₁ = polypoints[mod1(i+1, n_points)] - point
        rᵢ₋₁ = rᵢ
        rᵢ   = rᵢ₊₁
        rᵢ₊₁ = norm(sᵢ₊₁) 
        ## Now, we calculate the weight:
        wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
        ## perform a weighted sum with the interpolated value:
        interpolated_value += values[i] * wᵢ
        ## and add the weight to the total weight accumulator.
        wₜₒₜ += wᵢ
    end
    ## Return the normalized interpolated value.
    return interpolated_value / wₜₒₜ
end

# When you have holes, then you have to be careful 
# about the order you iterate around points.

# Specifically, you have to iterate around each linear ring separately 
# and ensure there are no degenerate/repeated points at the start and end!

function barycentric_interpolate(::MeanValue, exterior::AbstractVector{<: Point{N, T1}}, interiors::AbstractVector{<: AbstractVector{<: Point{N, T1}}}, values::AbstractVector{V}, point::Point{N, T2}) where {N, T1 <: Real, T2 <: Real, V}
    ## @boundscheck @assert length(values) == (length(exterior) + isempty(interiors) ? 0 : sum(length.(interiors)))
    ## @boundscheck @assert length(exterior) >= 3

    current_index = 1
    l_exterior = length(exterior)

    sᵢ₋₁ = exterior[end] - point
    sᵢ   = exterior[begin] - point
    sᵢ₊₁ = exterior[begin+1] - point
    rᵢ₋₁ = norm(sᵢ₋₁) # radius / Euclidean distance between points.
    rᵢ   = norm(sᵢ  ) # radius / Euclidean distance between points.
    rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
    # Now, we set the interpolated value to the first point's value, multiplied
    # by the weight computed relative to the first point in the polygon. 
    wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
    wₜₒₜ = wᵢ
    interpolated_value = values[begin] * wᵢ

    for i in 2:l_exterior
        # Increment counters + set variables
        sᵢ₋₁ = sᵢ
        sᵢ   = sᵢ₊₁
        sᵢ₊₁ = exterior[mod1(i+1, l_exterior)] - point
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
    for hole in interiors
        l_hole = length(hole)
        sᵢ₋₁ = hole[end] - point
        sᵢ   = hole[begin] - point
        sᵢ₊₁ = hole[begin+1] - point
        rᵢ₋₁ = norm(sᵢ₋₁) # radius / Euclidean distance between points.
        rᵢ   = norm(sᵢ  ) # radius / Euclidean distance between points.
        rᵢ₊₁ = norm(sᵢ₊₁) # radius / Euclidean distance between points.
        ## Now, we set the interpolated value to the first point's value, multiplied
        ## by the weight computed relative to the first point in the polygon. 
        wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ

        interpolated_value += values[current_index] * wᵢ

        wₜₒₜ += wᵢ
        current_index += 1
    
        for i in 2:l_hole
            ## Increment counters + set variables
            sᵢ₋₁ = sᵢ
            sᵢ   = sᵢ₊₁
            sᵢ₊₁ = hole[mod1(i+1, l_hole)] - point
            rᵢ₋₁ = rᵢ
            rᵢ   = rᵢ₊₁
            rᵢ₊₁ = norm(sᵢ₊₁) ## radius / Euclidean distance between points.
            wᵢ = (t_value(sᵢ₋₁, sᵢ, rᵢ₋₁, rᵢ) + t_value(sᵢ, sᵢ₊₁, rᵢ, rᵢ₊₁)) / rᵢ 
            interpolated_value += values[current_index] * wᵢ
            wₜₒₜ += wᵢ
            current_index += 1
        end
    end
    return interpolated_value / wₜₒₜ

end

struct Wachspress <: AbstractBarycentricCoordinateMethod
end
