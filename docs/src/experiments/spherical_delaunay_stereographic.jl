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

function spherical_triangulation(::SphericalConvexHull, input_points; facetype = CairoMakie.GeometryBasics.TriangleFace)
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




# These points are known to be good points, i.e., lon, lat, alt
points = Point3{Float64}.(JSON3.read(read(Downloads.download("https://gist.githubusercontent.com/Fil/6bc12c535edc3602813a6ef2d1c73891/raw/3ae88bf307e740ddc020303ea95d7d2ecdec0d19/points.json"), String)))
faces = delaunay_triangulate_spherical(points)
# or
# faces = Quickhull.quickhull(map(UnitCartesianFromGeographic(), points)) |> Quickhull.faces |> collect
# @assert isempty(setdiff(unique!(sort!(reduce(vcat, faces))), 1:length(points)))

# This is the super-cool scrollable 3D globe (though it's a bit deformed...)
f, a, p = Makie.mesh(map(UnitCartesianFromGeographic(), points), faces; color = last.(points), colormap = Reverse(:RdBu), colorrange = (-20, 40), shading = NoShading)

# We can also replicate the observable notebook almost exactly (just missing ExactPredicates):
f, a, p = Makie.mesh(points, faces; axis = (; type = GeoAxis, dest = "+proj=bertin1953 +lon_0=-16.5 +lat_0=-42 +x_0=7.93 +y_0=0.09"), color = last.(points), colormap = Reverse(:RdBu), colorrange = (-20, 40), shading = NoShading)
# Whoops, this doesn't look so good!  Let's try to do this more "manually" instead.
# We'll use NaturalNeighbours.jl for this, but we first have to reconstruct the triangulation
# with the same faces, but in a Bertin projection...
using NaturalNeighbours
lonlat2bertin = Proj.Transformation(GFT.EPSG(4326), GFT.ProjString("+proj=bertin1953 +type=crs +lon_0=-16.5 +lat_0=-42 +x_0=7.93 +y_0=0.09"); always_xy = true)
lons = LinRange(-180, 180, 300)
lats = LinRange(-90, 90, 150)
bertin_points = lonlat2bertin.(lons, lats')

projected_points = GO.reproject(GO.tuples(points), source_crs = GFT.EPSG(4326), target_crs = "+proj=bertin1953 +lon_0=-16.5 +lat_0=-42 +x_0=7.93 +y_0=0.09")

ch = DelTri.convex_hull(projected_points) # assumes each point is in the triangulation
boundary_nodes = DelTri.get_vertices(ch) 
bertin_boundary_poly = GI.Polygon([GI.LineString(DelTri.get_points(ch)[DelTri.get_vertices(ch)])])

tri = DelTri.Triangulation(projected_points, faces, boundary_nodes)
itp = NaturalNeighbours.interpolate(tri, last.(points); derivatives = true)

mat = [
    if GO.contains(bertin_boundary_poly, (x, y))
        itp(x, y; method = Nearest(), project = false)
    else
        NaN
    end
    for (x, y) in bertin_points
]
# TODO: this currently doesn't work, because some points are not inside a triangle and so cannot be located.
# Options are:
# 1. Reject all points outside the convex hull of the projected points.  (Tried but failed)
# 2. ...





# Now, we proceed with the script implementation which has quite a bit of debug information + plots
pivot_ind = findfirst(isfinite, points)

# debug
point_colors = fill(:blue, length(points))
point_colors[pivot_ind] = :red
# end debug

pivot_point = points[pivot_ind]
necessary_rotation = #=Rotations.RotY(-π) *=# Rotations.RotY(-deg2rad(90-pivot_point[2])) * Rotations.RotZ(-deg2rad(pivot_point[1]))
#
net_transformation_to_corrected_cartesian = CoordinateTransformations.LinearMap(necessary_rotation) ∘ UnitCartesianFromGeographic()

scatter(map(net_transformation_to_corrected_cartesian, points); color = point_colors)
#
stereographic_points = (StereographicFromCartesian() ∘ net_transformation_to_corrected_cartesian).(points)
scatter(stereographic_points; color = point_colors)
#
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
@show max_radius length(really_far_idxs)
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

boundary_nodes, pts = DelTri.convert_boundary_points_to_indices(boundary_points; existing_points=triangulation_points)

triangulation = DelTri.triangulate(pts; boundary_nodes)
# triangulation = DelTri.triangulate(triangulation_points)
triplot(triangulation)
DelTri.validate_triangulation(triangulation)
#  for diagnostics, try `fig, ax, sc = triplot(triangulation, show_constrained_edges=true, show_convex_hull=true)`

# First, get all the "solid" faces, ie, faces not attached to boundary nodes
original_triangles = collect(DelTri.each_solid_triangle(triangulation))
boundary_face_inds = findall(Base.Fix1(DelTri.is_boundary_triangle, triangulation), original_triangles)
faces = map(CairoMakie.GeometryBasics.TriangleFace, view(original_triangles, setdiff(axes(original_triangles, 1), boundary_face_inds)))

for boundary_face in view(original_triangles, boundary_face_inds)
    push!(faces, CairoMakie.GeometryBasics.TriangleFace(map(boundary_face) do i; first(DelTri.is_boundary_node(triangulation, i)) ? pivot_ind : i end))
end

for ghost_face in DelTri.each_ghost_triangle(triangulation)
    push!(faces, CairoMakie.GeometryBasics.TriangleFace(map(ghost_face) do i; DelTri.is_ghost_vertex(i) ? pivot_ind : i end))
end
# Remove degenerate triangles
filter!(faces) do face
    !(face[1] == face[2] || face[2] == face[3] || face[1] == face[3])
end

wireframe(CairoMakie.GeometryBasics.normal_mesh(map(UnitCartesianFromGeographic(), points), faces))

grid_lons = LinRange(-180, 180, 10)
grid_lats = LinRange(-90, 90, 10)
lon_grid = [CairoMakie.GeometryBasics.LineString([Point{2, Float64}((grid_lons[j], -90)), Point{2, Float64}((grid_lons[j], 90))]) for j in 1:10] |> vec
lat_grid = [CairoMakie.GeometryBasics.LineString([Point{2, Float64}((-180, grid_lats[i])), Point{2, Float64}((180, grid_lats[i]))]) for i in 1:10] |> vec

sampled_grid = GO.segmentize(lat_grid; max_distance = 0.01)

geographic_grid = GO.transform(UnitCartesianFromGeographic(), sampled_grid) .|> x -> GI.LineString{true, false}(x.geom)
stereographic_grid = GO.transform(sampled_grid) do point
    (StereographicFromCartesian() ∘ UnitCartesianFromGeographic())(point)
end .|> x -> GI.LineString{false, false}(x.geom)



fig = Figure(); ax = LScene(fig[1, 1])
for (i, line) in enumerate(geographic_grid)
    lines!(ax, line; linewidth = 3, color = Makie.resample_cmap(:viridis, length(stereographic_grid))[i])
end
fig

#
fig = Figure(); ax = LScene(fig[1, 1])
for (i, line) in enumerate(stereographic_grid)
    lines!(ax, line; linewidth = 3, color = Makie.resample_cmap(:viridis, length(stereographic_grid))[i])
end
fig

all(reduce(vcat, faces) |> sort! |> unique! .== 1:length(points))

faces[findall(x -> 1 in x, faces)]

# try NaturalNeighbours, fail miserably
using NaturalNeighbours

itp = NaturalNeighbours.interpolate(triangulation, last.(points); derivatives = true)
lons = LinRange(-180, 180, 360)
lats = LinRange(-90, 90, 180)
_x = vec([x for x in lons, _ in lats])
_y = vec([y for _ in lons, y in lats])

sibson_1_vals = itp(_x, _y; method=Sibson(1))

transformed_points = GO.transform(points) do point
    # first, spherical to Cartesian
    longitude, latitude = deg2rad(GI.x(point)), deg2rad(GI.y(point))
    radius = 1 # we will operate on the unit sphere
    xyz = Point3d(
        radius * sin(latitude) * cos(longitude),
        radius * sin(latitude) * sin(longitude),
        radius * cos(latitude)
    )
    # then, rotate Cartesian so the pivot point is at the south pole
    rotated_xyz = Rotations.AngleAxis
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

θϕs = randsphericalangles(50)
θs, ϕs = first.(θϕs), last.(θϕs)
pts = Point3.(
    sin.(ϕs) .* cos.(θs),
    sin.(ϕs) .* sin.(θs),
    cos.(ϕs)
)

f, a, p = scatter(pts; color = [fill(:blue, 49); :red])

function Makie.rotate!(t::Makie.Transformable, rotation::Rotations.Rotation)
    quat = Rotations.QuatRotation(rotation)
    Makie.rotate!(Makie.Absolute, t, Makie.Quaternion(quat.q.v1, quat.q.v2, quat.q.v3, quat.q.s))
end


rotate!(p, Rotations.RotX(π/2))
rotate!(p, Rotations.RotX(0))

pivot_point = θϕs[end]
necessary_rotation = Rotations.RotY(pi) * Rotations.RotY(-pivot_point[2]) * Rotations.RotZ(-pivot_point[1])

rotate!(p, necessary_rotation)

f