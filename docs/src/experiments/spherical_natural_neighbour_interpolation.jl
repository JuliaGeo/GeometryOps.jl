#=
# Spherical natural neighbour interpolatoin
=#

import GeometryOps as GO, GeoInterface as GI, GeoFormatTypes as GFT
import Proj # for easy stereographic projection - TODO implement in Julia
import DelaunayTriangulation as DelTri # Delaunay triangulation on the 2d plane
import CoordinateTransformations, Rotations

using Downloads # does what it says on the tin
using JSON3 # to load data
using CairoMakie, GeoMakie # for plotting
import Makie: Point3d

include(joinpath(@__DIR__, "spherical_delaunay.jl"))

using LinearAlgebra
using GeometryBasics
struct SphericalCap{T}
    point::Point3{T}
    radius::T
end

function SphericalCap(point::Point3{T1}, radius::T2) where {T1, T2}
    return SphericalCap{promote_type(T1, T2)}(point, radius)
end


function circumcenter_on_unit_sphere(a, b, c)
    LinearAlgebra.normalize(a × b + b × c + c × a)
end

spherical_distance(x::Point3, y::Point3) = acos(clamp(x ⋅ y, -1.0, 1.0))

"Get the circumcenter of the triangle (a, b, c) on the unit sphere.  Returns a normalized 3-vector."
function SphericalCap(a::Point3, b::Point3, c::Point3)
    circumcenter = circumcenter_on_unit_sphere(a, b, c)
    circumradius = spherical_distance(a, circumcenter)
    return SphericalCap(circumcenter, circumradius)
end

function _is_ccw_unit_sphere(v_0,v_c,v_i)
    # checks if the smaller interior angle for the great circles connecting u-v and v-w is CCW
    return(LinearAlgebra.dot(LinearAlgebra.cross(v_c - v_0,v_i - v_c), v_i) < 0)
end

function angle_between(a, b, c)
    ab = b - a
    bc = c - b
    norm_dot = (ab ⋅ bc) / (LinearAlgebra.norm(ab) * LinearAlgebra.norm(bc))
    angle =  acos(clamp(norm_dot, -1.0, 1.0))
    if _is_ccw_unit_sphere(a, b, c)
        return angle
    else
        return 2π - angle
    end
end

function bowyer_watson_envelope!(applicable_points, query_point, points, faces, caps = map(splat(SphericalCap), (view(points, face) for face in faces)); applicable_cap_indices = Int64[])
    # brute force for now, but try the jump and search algorithm later
    # can use e.g GeometryBasics.decompose(Point3{Float64}, GeometryBasics.Sphere(Point3(0.0), 1.0), 5) 
    # to get starting points, or similar
    empty!(applicable_cap_indices)
    for (i, cap) in enumerate(caps)
        if cap.radius ≥ spherical_distance(query_point, cap.point)
            push!(applicable_cap_indices, i)
        end
    end
    # Now that we have the face indices, we need to get the applicable points
    empty!(applicable_points)
    # for face_idx in applicable_cap_indices
    #     append!(applicable_points, faces[face_idx])
    # end
    for i in applicable_cap_indices
        current_face = faces[i]
        edge_reoccurs = false
        for current_edge in ((current_face[1], current_face[2]), (current_face[2], current_face[3]), (current_face[3], current_face[1]))
            for j in applicable_cap_indices
                if j == i
                    continue # can't compare a triangle to itself
                end
                face_to_compare = faces[j]
                for edge_to_compare in ((face_to_compare[1], face_to_compare[2]), (face_to_compare[2], face_to_compare[3]), (face_to_compare[3], face_to_compare[1]))
                    if edge_to_compare == current_edge || reverse(edge_to_compare) == current_edge
                        edge_reoccurs = true
                        break
                    end
                end
            end
            if !edge_reoccurs # edge is unique
                push!(applicable_points, current_edge[1])
                push!(applicable_points, current_edge[2])
            end
            edge_reoccurs = false
        end
    end
    unique!(applicable_points)
    # Pick a random point (in this case, the first point) and 
    # compute the angle subtended by the first point, the query
    # point, and each individual point.  Sorting these angles
    # allows us to "wind" the polygon in a clockwise way.
    three_d_points = view(points, applicable_points)
    angles = [angle_between(three_d_points[1], query_point, point) for point in three_d_points]
    pt_inds = sortperm(angles)
    permute!(applicable_points, pt_inds)
    # push!(applicable_points, applicable_points[begin])
    return applicable_points
end

import NaturalNeighbours: previndex_circular, nextindex_circular
function laplace_ratio(points, envelope, i #= current vertex index =#, interpolation_point)
    u = envelope[i]
    prev_u = envelope[previndex_circular(envelope, i)]
    next_u = envelope[nextindex_circular(envelope, i)]
    p, q, r = points[u], points[prev_u], points[next_u]
    g1 = circumcenter_on_unit_sphere(q, p, interpolation_point)
    g2 = circumcenter_on_unit_sphere(p, r, interpolation_point)
    ℓ = spherical_distance(g1, g2)
    d = spherical_distance(p, interpolation_point)
    w = ℓ / d
    return w, u, prev_u, next_u
end

struct NaturalCoordinates{F, I}
    coordinates::Vector{F}
    indices::Vector{I}
    interpolation_point::Point3{Float64}
end

function laplace_nearest_neighbour_coords(points, faces, interpolation_point, caps = map(splat(SphericalCap), (view(points, face) for face in faces)); envelope = Int64[], cap_idxs = Int64[])
    envelope = bowyer_watson_envelope!(envelope, interpolation_point, points, faces, caps; applicable_cap_indices = cap_idxs)
    coords = NaturalCoordinates(Float64[], Int64[], interpolation_point)
    for i in eachindex(envelope)
        w, u, prev_u, next_u = laplace_ratio(points, envelope, i, interpolation_point)
        push!(coords.coordinates, w)
        push!(coords.indices, u)
    end
    coords.coordinates ./= sum(coords.coordinates)
    return coords
end

function eval_laplace_coordinates(points, faces, zs, interpolation_point, caps = map(splat(SphericalCap), (view(points, face) for face in faces)); envelope = Int64[], cap_idxs = Int64[])
    coords = laplace_nearest_neighbour_coords(points, faces, interpolation_point, caps; envelope, cap_idxs)
    if isempty(coords.coordinates)
        return Inf
    end
    return sum((coord * zs[idx] for (coord, idx) in zip(coords.coordinates, coords.indices)))
end

function get_voronoi_polygon(point_index, points, faces, caps)
    # For a Voronoi polygon, we can simply retrieve all triangles that have
    # point_index as a vertex, and then filter them together.  
end



# These points are known to be good points, i.e., lon, lat, alt
geographic_points = Point3{Float64}.(JSON3.read(read(Downloads.download("https://gist.githubusercontent.com/Fil/6bc12c535edc3602813a6ef2d1c73891/raw/3ae88bf307e740ddc020303ea95d7d2ecdec0d19/points.json"), String)))
cartesian_points = UnitCartesianFromGeographic().(geographic_points)
z_values = last.(geographic_points)

faces = spherical_triangulation(geographic_points)

faces2 = spherical_triangulation(SphericalConvexHull(), geographic_points)

unique!(sort!(reduce(vcat, faces))) # so how am I getting this index?


caps = map(splat(SphericalCap), (view(cartesian_points, face) for face in faces))

lons = -180.0:1.0:180.0
lats = -90.0:1.0:90.0

eval_laplace_coordinates(cartesian_points, faces, z_values, LinearAlgebra.normalize(Point3(1.0, 1.0, 0.0)))
eval_laplace_coordinates(cartesian_points, faces2, z_values, LinearAlgebra.normalize(Point3(1.0, 1.0, 0.0)))

values = map(UnitCartesianFromGeographic().(Point2.(lons, lats'))) do point
    eval_laplace_coordinates(cartesian_points, faces, z_values, point)
end

heatmap(lons, lats, values; axis = (; aspect = DataAspect()))

f = Figure();
a = LScene(f[1, 1])
p = meshimage!(a, lons, lats, rotl90(values); npoints = (720, 360))
p.transformation.transform_func[] = Makie.PointTrans{3}(UnitCartesianFromGeographic())
# scatter!(a, cartesian_points)
f # not entirely sure what's going on here
# diagnostics
# f, a, p = scatter(reduce(vcat, (view(cartesian_points, face) for face in view(faces, neighbour_inds))))
# scatter!(query_point; color = :red, markersize = 40)

query_point = LinearAlgebra.normalize(Point3(1.0, 1.0, 0.5))
pt_inds = bowyer_watson_envelope!(Int64[], query_point, cartesian_points, faces, caps)
three_d_points = cartesian_points[pt_inds]
angles = [angle_between(three_d_points[1], query_point, point) for point in three_d_points]

# pushfirst!(angles, 0)
pt_inds[sortperm(angles)]
f, a, p = scatter([query_point]; markersize = 30, color = :green, axis = (; type = LScene));
scatter!(a, cartesian_points)
scatter!(a, view(cartesian_points, pt_inds[sortperm(angles)]); color = eachindex(pt_inds), colormap = :RdBu, markersize = 20)
wireframe!(a, GeometryBasics.Mesh(cartesian_points, faces); alpha = 0.3)

f


function Makie.convert_arguments(::Type{Makie.Mesh}, cap::SphericalCap)
    offset_point = LinearAlgebra.normalize(cap.point + LinearAlgebra.normalize(Point3(cos(cap.radius), sin(cap.radius), 0.0)))
    points = [cap.point + Rotations.AngleAxis(θ, cap.point...) * offset_point for θ in LinRange(0, 2π, 20)]
    push!(points, cap.point)
    faces = [GeometryBasics.TriangleFace(i, i+1, 21) for i in 1:19]
    return (GeometryBasics.normal_mesh(points, faces),)
end


# Animate the natural neighbours of a loxodromic curve
# To create a loxodromic curve we use a Mercator projection
import GeometryOps as GO, GeoInterface as GI, GeoFormatTypes as GFT
using CairoMakie, GeoMakie
start_point = (-95.7129, 90.0)
end_point = (-130.0, 0.0)

geographic_path = GI.LineString([start_point, end_point]; crs = GFT.EPSG(4326))
mercator_path = GO.reproject(geographic_path; target_crs = GFT.EPSG(3857))
sampled_loxodromic_path = GO.segmentize(GO.LinearSegments(; max_distance = 1e5), mercator_path) # max distance in mercator coordinates!
sampled_loxodromic_cartesian_path = GO.transform(UnitCartesianFromGeographic(), GO.reproject(sampled_loxodromic_path; target_crs = GFT.EPSG(4326))) |> GI.getpoint

lines(sampled_loxodromic_cartesian_path)

f, a, p = lines(sampled_loxodromic_path; axis = (; type = GeoAxis))
lines!(a, GeoMakie.coastlines(); color = (:black, 0.5))
f



# A minimal test case to debug
randpoints = randsphere(50)
push!(randpoints, LinearAlgebra.normalize(Point3d(0.0, 1.0, 1.0)))

zs = [zeros(length(randpoints)-1); 1.0]

faces = spherical_triangulation(SphericalConvexHull(), randpoints)
caps = map(splat(SphericalCap), (view(randpoints, face) for face in faces))
values = map(UnitCartesianFromGeographic().(Point2.(lons, lats'))) do query_point
    eval_laplace_coordinates(randpoints, faces, zs, query_point, caps)
end

scatter(randpoints; color = zs)

heatmap(lons, lats, values)


f = Figure();
a = LScene(f[1, 1])
p = meshimage!(a, lons, lats, rotl90(values); npoints = (720, 360))
p.transformation.transform_func[] = Makie.PointTrans{3}(UnitCartesianFromGeographic())
# scatter!(a, cartesian_points)
f # not entirely sure what's going on here

f, a, p = wireframe(GeometryBasics.Mesh(randpoints, faces))

struct UnitSphereSegments <: GO.SegmentizeMethod
    "The maximum distance between segments, in radians (solid angle units).  A nice value here for the globe is ~0.01." 
    max_distance::Float64
end

function GO._fill_linear_kernel!(p, q, npoints = 10)
    # perform slerp interpolation between p and q
    # return npoints points on the great circle between p and q
    # using spherical linear interpolation
    # return the points as a vector of Point3{Float64}
    ω = acos(clamp(dot(p, q), -1.0, 1.0))
    if isapprox(ω, 0, atol=1e-8)
        return fill(p, (npoints,))
    end
    return [
        (sin((1 - t) * ω) * p + sin(t * ω) * q) / sin(ω)
        for t in range(0, 1, length=npoints)
    ]
end

function _segmentize(method::UnitSphereSegments, geom, T::Union{GI.LineStringTrait, GI.LinearRingTrait})
    first_coord = GI.getpoint(geom, 1)
    x1, y1, z1 = GI.x(first_coord), GI.y(first_coord), GI.z(first_coord)
    new_coords = NTuple{3, Float64}[]
    sizehint!(new_coords, GI.npoint(geom))
    push!(new_coords, (x1, y1, z1   ))
    for coord in Iterators.drop(GI.getpoint(geom), 1)
        x2, y2, z2 = GI.x(coord), GI.y(coord), GI.z(coord)
        _fill_linear_kernel!(method, new_coords, x1, y1, z1, x2, y2, z2)
        x1, y1, z1 = x2, y2, z2
    end 
    return rebuild(geom, new_coords)
end

function get_face_lines(points, faces)
    return GeometryBasics.MultiLineString(map(faces) do face
        return GeometryBasics.LineString(vcat(
            interpolate_line_unit_sphere(points[face[1]], points[face[2]], 10),
            interpolate_line_unit_sphere(points[face[2]], points[face[3]], 10),
            interpolate_line_unit_sphere(points[face[3]], points[face[1]], 10),
        ))
    end)
end

lines(get_face_lines(cartesian_points, faces))