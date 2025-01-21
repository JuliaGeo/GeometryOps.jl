
include(joinpath(@__DIR__, "spherical_delaunay.jl"))

# These points are known to be good points, i.e., lon, lat, alt
points = Point3{Float64}.(JSON3.read(read(Downloads.download("https://gist.githubusercontent.com/Fil/6bc12c535edc3602813a6ef2d1c73891/raw/3ae88bf307e740ddc020303ea95d7d2ecdec0d19/points.json"), String)))
faces = spherical_triangulation(points)
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


