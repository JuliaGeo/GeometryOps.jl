#=
# Spherical Delaunay triangulation using stereographic projections

This file encodes two approaches to spherical Delaunay triangulation.

1. `StereographicDelaunayTriangulation` projects the coordinates to a rotated stereographic projection, then computes the Delaunay triangulation on the plane.
    This approach was taken from d3-geo-voronoi.
    ```
    Copyright 2018-2021 Philippe Rivière

    Permission to use, copy, modify, and/or distribute this software for any purpose
    with or without fee is hereby granted, provided that the above copyright notice
    and this permission notice appear in all copies.

    THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
    REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
    FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
    INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
    OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
    TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
    THIS SOFTWARE.
    ```
2. `SphericalConvexHull` computes the convex hull of the points in 3D Cartesian space, which is by definition the Delaunay triangulation.  This triangulation is currently done using the unreleased Quickhull.jl.
=#

import GeometryOps as GO, GeoInterface as GI, GeoFormatTypes as GFT
import Proj # for easy stereographic projection - TODO implement in Julia
import DelaunayTriangulation as DelTri # Delaunay triangulation on the 2d plane
import Quickhull # convex hulls in d+1 are Delaunay triangulations in dimension d
import CoordinateTransformations, Rotations

using Downloads # does what it says on the tin
using JSON3 # to load data
using CairoMakie, GeoMakie # for plotting
import Makie: Point3d

abstract type SphericalTriangulationAlgorithm end
struct StereographicDelaunayTriangulation <: SphericalTriangulationAlgorithm end
struct SphericalConvexHull <: SphericalTriangulationAlgorithm end 

spherical_triangulation(input_points; kwargs...) = spherical_triangulation(StereographicDelaunayTriangulation(), input_points; kwargs...)

function spherical_triangulation(::StereographicDelaunayTriangulation, input_points; facetype = CairoMakie.GeometryBasics.TriangleFace)
    # @assert GI.crstrait(first(input_points)) isa GI.AbstractGeographicCRSTrait
    points = GO.tuples(input_points)
    # @assert points isa Vector{GI.Point}

    # In 
    pivot_ind = findfirst(x -> all(isfinite, x), points)
    
    pivot_point = points[pivot_ind]
    necessary_rotation = #=Rotations.RotY(-π) *=# Rotations.RotY(-deg2rad(90-pivot_point[2])) * Rotations.RotZ(-deg2rad(pivot_point[1]))

    net_transformation_to_corrected_cartesian = CoordinateTransformations.LinearMap(necessary_rotation) ∘ UnitCartesianFromGeographic()
    stereographic_points = (StereographicFromCartesian() ∘ net_transformation_to_corrected_cartesian).(points)
    triangulation_points = copy(stereographic_points)
    
    # Iterate through the points and find infinite/invalid points
    max_radius = 1
    really_far_idxs = Int[]
    for (i, point) in enumerate(stereographic_points)
        m = hypot(point[1], point[2])
        if !isfinite(m) || m > 1e32
            push!(really_far_idxs, i)
        elseif m > max_radius
            max_radius = m
        end
    end
    @debug max_radius length(really_far_idxs)
    far_value = 1e6 * sqrt(max_radius)

    if !isempty(really_far_idxs)
        triangulation_points[really_far_idxs] .= (Point2(far_value, 0.0),)
    end

    boundary_points = reverse([
        (-far_value, -far_value),
        (-far_value, far_value),
        (far_value, far_value),
        (far_value, -far_value),
        (-far_value, -far_value),
    ])

    # boundary_nodes, pts = DelTri.convert_boundary_points_to_indices(boundary_points; existing_points=triangulation_points)

    # triangulation = DelTri.triangulate(triangulation_points; boundary_nodes)
    triangulation = DelTri.triangulate(triangulation_points)
    #  for diagnostics, try `fig, ax, sc = triplot(triangulation, show_constrained_edges=true, show_convex_hull=true)`

    # First, get all the "solid" faces, ie, faces not attached to boundary nodes
    original_triangles = collect(DelTri.each_solid_triangle(triangulation))
    boundary_face_inds = findall(Base.Fix1(DelTri.is_boundary_triangle, triangulation), original_triangles)
    faces = map(facetype, view(original_triangles, setdiff(axes(original_triangles, 1), boundary_face_inds)))
    
    for boundary_face in view(original_triangles, boundary_face_inds)
        push!(faces, facetype(map(boundary_face) do i; first(DelTri.is_boundary_node(triangulation, i)) ? pivot_ind : i end))
    end
    
    for ghost_face in DelTri.each_ghost_triangle(triangulation)
        push!(faces, facetype(map(ghost_face) do i; DelTri.is_ghost_vertex(i) ? pivot_ind : i end))
    end
    # Remove degenerate triangles
    filter!(faces) do face
        !(DelTri.geti(face) == DelTri.getj(face) || DelTri.getj(face) == DelTri.getk(face) || DelTri.geti(face) == DelTri.getk(face))
    end

    return faces
end

function spherical_triangulation(::SphericalConvexHull, input_points; facetype = CairoMakie.GeometryBasics.TriangleFace{Int32})
    points = GO.tuples(input_points) # we have to decompose the points into tuples, so they work with Quickhull.jl
    # @assert points isa Vector{GI.Point}
    cartesian_points = map(UnitCartesianFromGeographic(), points)
    # The Delaunay triangulation of points on a sphere is simply the convex hull of those points in 3D Cartesian space.
    # We can use e.g Quickhull.jl to get us such a convex hull.
    hull = Quickhull.quickhull(cartesian_points)
    # We return only the faces from these triangulation methods, so we simply map
    # the facetype to the returned values from `Quickhull.facets`.
    return map(facetype, Quickhull.facets(hull))
end


# necessary coordinate transformations

struct StereographicFromCartesian <: CoordinateTransformations.Transformation
end

function (::StereographicFromCartesian)(xyz::AbstractVector)
    @assert length(xyz) == 3 "StereographicFromCartesian expects a 3D Cartesian vector"
    x, y, z = xyz
    # The Wikipedia definition has the north pole at infinity,
    # this implementation has the south pole at infinity.
    return Point2(x/(1-z), y/(1-z))
end

struct CartesianFromStereographic <: CoordinateTransformations.Transformation
end

function (::CartesianFromStereographic)(stereographic_point)
    X, Y = stereographic_point
    x2y2_1 = X^2 + Y^2 + 1
    return Point3(2X/x2y2_1, 2Y/x2y2_1, (x2y2_1 - 2)/x2y2_1)
end

struct UnitCartesianFromGeographic <: CoordinateTransformations.Transformation 
end

function (::UnitCartesianFromGeographic)(geographic_point)
    # Longitude is directly translatable to a spherical coordinate
    # θ (azimuth)
    θ = deg2rad(GI.x(geographic_point))
    # The polar angle is 90 degrees minus the latitude
    # ϕ (polar angle)
    ϕ = deg2rad(90 - GI.y(geographic_point))
    # Since this is the unit sphere, the radius is assumed to be 1,
    # and we don't need to multiply by it.
    return Point3(
        sin(ϕ) * cos(θ),
        sin(ϕ) * sin(θ),
        cos(ϕ)
    )
end

struct GeographicFromUnitCartesian <: CoordinateTransformations.Transformation 
end

function (::GeographicFromUnitCartesian)(xyz::AbstractVector)
    @assert length(xyz) == 3 "GeographicFromUnitCartesian expects a 3D Cartesian vector"
    x, y, z = xyz
    return Point2(
        atan(y, x),
        atan(hypot(x, y), z),
    )
end

function randsphericalangles(n)
    θ = 2π .* rand(n)
    ϕ = acos.(2 .* rand(n) .- 1)
    return Point2.(θ, ϕ)
end

function randsphere(n)
    θϕ = randsphericalangles(n)
    return Point3.(
        sin.(last.(θϕs)) .* cos.(first.(θϕs)),
        sin.(last.(θϕs)) .* sin.(first.(θϕs)),
        cos.(last.(θϕs))
    )
end

# θϕs = randsphericalangles(50)
# θs, ϕs = first.(θϕs), last.(θϕs)
# pts = Point3.(
#     sin.(ϕs) .* cos.(θs),
#     sin.(ϕs) .* sin.(θs),
#     cos.(ϕs)
# )

# f, a, p = scatter(pts; color = [fill(:blue, 49); :red])

# function Makie.rotate!(t::Makie.Transformable, rotation::Rotations.Rotation)
#     quat = Rotations.QuatRotation(rotation)
#     Makie.rotate!(Makie.Absolute, t, Makie.Quaternion(quat.q.v1, quat.q.v2, quat.q.v3, quat.q.s))
# end


# rotate!(p, Rotations.RotX(π/2))
# rotate!(p, Rotations.RotX(0))

# pivot_point = θϕs[end]
# necessary_rotation = Rotations.RotY(pi) * Rotations.RotY(-pivot_point[2]) * Rotations.RotZ(-pivot_point[1])

# rotate!(p, necessary_rotation)

# f