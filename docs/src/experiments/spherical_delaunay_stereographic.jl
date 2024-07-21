#=
# Spherical Delaunay triangulation using stereographic projections

This is the approach which d3-geo-voronoi and friends use.  The alternative is STRIPACK which computes in 3D on the sphere.
The 3D approach is basically that the 3d convex hull of a set of points on the sphere is its Delaunay triangulation.
=#

import GeometryOps as GO, GeoInterface as GI
import Proj # for easy stereographic projection - TODO implement in Julia
import DelaunayTriangulation as DelTri # Delaunay triangulation on the 2d plane
import CoordinateTransformations, Rotations

using Downloads # does what it says on the tin
using JSON3 # to load data
using CairoMakie # for plotting
import Makie: Point3d

points = Point3{Float64}.(JSON3.read(read(Downloads.download("https://gist.githubusercontent.com/Fil/6bc12c535edc3602813a6ef2d1c73891/raw/3ae88bf307e740ddc020303ea95d7d2ecdec0d19/points.json"), String)))

pivot_ind = findfirst(isfinite, points)
pivot_point = points[pivot_ind]
necessary_rotation = Rotations.RotY(pi) * Rotations.RotY(-pivot_point[2]) * Rotations.RotZ(-pivot_point[1])

net_transformation_to_corrected_cartesian = CoordinateTransformations.LinearMap(necessary_rotation) ∘ UnitCartesianFromGeographic()

stereographic_points = (StereographicFromCartesian() ∘ net_transformation_to_corrected_cartesian).(points)

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
far_value = 1e6 * sqrt(max_radius)
if !isempty(really_far_idxs)
    stereographic_points[really_far_idxs] .= Point2(far_value, 0)
end

# Add infinite horizon points
triangulation_points = copy(stereographic_points)
push!(triangulation_points, Point2(0, far_value))
push!(triangulation_points, Point2(-far_value, 0))
push!(triangulation_points, Point2(0, -far_value))

triangulation = DelTri.triangulate(triangulation_points)

triangulation.points[1:end-3] .= stereographic_points
triangulation.points[end-2:end] .= (stereographic_points[pivot_ind],)

f, a, p = triplot(triangulation; axis = (; type = Axis3,));
m = p.plots[1].plots[1][1][]

geo_mesh_points = vcat(points, repeat([pivot_point], 3))
cartesian_mesh_points = UnitCartesianFromGeographic().(geo_mesh_points)
f, a, p = mesh(cartesian_mesh_points, getfield(getfield(m, :simplices), :faces))
wireframe(p.mesh[])
CairoMakie.GeometryBasics.Mesh(UnitCartesianFromGeographic().(vcat(points, repeat([pivot_point], 3))), splat(CairoMakie.GeometryBasics.TriangleFace).(collect(triangulation.triangles))) |> wireframe

geo_triangulation = deepcopy(triangulation)
geo_triangulation.points .= geo_mesh_points

itp = NaturalNeighbours.interpolate(geo_triangulation, last.(geo_mesh_points); derivatives = true)
lons = LinRange(-180, 180, 360)
lats = LinRange(-90, 90, 180)
_x = vec([x for x in lons, _ in lats])
_y = vec([y for _ in lons, y in lats])

sibson_1_vals = itp(_x, _y; method=Sibson(1))


struct StereographicFromCartesian <: CoordinateTransformations.Transformation
end

function (::StereographicFromCartesian)(xyz::AbstractVector)
    @assert length(xyz) == 3 "StereographicFromCartesian expects a 3D Cartesian vector"
    x, y, z = xyz
    return Point2(x/(1-z), y/(1-z))
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

pivot = θϕs[end]
rotate!(p, )

