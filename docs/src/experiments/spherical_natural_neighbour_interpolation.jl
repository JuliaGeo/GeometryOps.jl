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

# include(joinpath(@__DIR__, "spherical_delaunay_stereographic.jl"))

using LinearAlgebra

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

function bowyer_watson_envelope!(applicable_cap_indices, query_point, points, faces, caps = map(splat(SphericalCap), (view(cartesian_points, face) for face in faces)))
    # brute force for now, but try the jump and search algorithm later
    # can use e.g GeometryBasics.decompose(Point3{Float64}, GeometryBasics.Sphere(Point3(0.0), 1.0), 5) 
    # to get starting points, or similar
    empty!(applicable_cap_indices)
    for (i, cap) in enumerate(caps)
        if cap.radius > spherical_distance(query_point, cap.point)
            push!(applicable_cap_indices, i)
        end
    end
    # Now that we have the face indices, we need to get the applicable points
    applicable_points = Int64[]
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
        end
    end
    return unique!(applicable_points)
end







# These points are known to be good points, i.e., lon, lat, alt
geographic_points = Point3{Float64}.(JSON3.read(read(Downloads.download("https://gist.githubusercontent.com/Fil/6bc12c535edc3602813a6ef2d1c73891/raw/3ae88bf307e740ddc020303ea95d7d2ecdec0d19/points.json"), String)))
z_values = last.(geographic_points)
faces = spherical_triangulation(geographic_points)
# correct the faces, since the order seems to be off
faces = reverse.(faces)

unique!(sort!(reduce(vcat, faces))) # so how am I getting this index?

cartesian_points = UnitCartesianFromGeographic().(geographic_points)

caps = map(splat(SphericalCap), (view(cartesian_points, face) for face in faces))

lons = -180.0:180.0
lats = -90.0:90.0

eval_laplace_coords(cartesian_points, z_values, faces, Point3(1.0, 0.0, 0.0))

values = map(UnitCartesianFromGeographic().(Point2.(lons, lats'))) do point
    eval_laplace_coords(cartesian_points, z_values, faces, point)
end

heatmap(lons, lats, values)

# diagnostics
# f, a, p = scatter(reduce(vcat, (view(cartesian_points, face) for face in view(faces, neighbour_inds))))
# scatter!(query_point; color = :red, markersize = 40)
import NaturalNeighbours: previndex_circular, nextindex_circular
function laplace_ratio(points, envelope, i #= current vertex index =#, interpolation_point)
    u = envelope[i]
    prev_u = envelope[previndex_circular(envelope, i)]
    next_u = envelope[nextindex_circular(envelope, i)]
    p = points[u]
    q, r = points[prev_u], points[next_u]
    g1 = circumcenter_on_unit_sphere(q, p, interpolation_point)
    g2 = circumcenter_on_unit_sphere(p, r, interpolation_point)
    ℓ = spherical_distance(g1, g2)
    d = spherical_distance(p, interpolation_point)
    w = ℓ / d
    return w, u, prev_u, next_u
end

function eval_laplace_coords(points, zs, faces, interpolation_point; envelope = Int64[])
    envelope = bowyer_watson_envelope!(envelope, interpolation_point, points, faces, caps)
    weighted_sum = 0.0
    weight = 0.0
    for i in eachindex(envelope)
        w, u, prev_u, next_u = laplace_ratio(points, envelope, i, interpolation_point)
        weighted_sum += w * zs[i]
        weight += w
    end
    return weighted_sum / weight
end