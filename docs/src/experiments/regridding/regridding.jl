import GeoInterface as GI, GeometryOps as GO
using SortTileRecursiveTree: STRtree
using SparseArrays: spzeros
using Extents

using CairoMakie, GeoInterfaceMakie

include("sphericalpoints.jl")

function area_of_intersection_operator(grid1, grid2; nodecapacity1 = 10, nodecapacity2 = 10)
    area_of_intersection_operator(GO.Planar(), grid1, grid2; nodecapacity1 = nodecapacity1, nodecapacity2 = nodecapacity2)
end

function area_of_intersection_operator(m::GO.Manifold, grid1, grid2; nodecapacity1 = 10, nodecapacity2 = 10) # grid1 and grid2 are both vectors of polygons
    A = spzeros(Float64, length(grid1), length(grid2))
    # Prepare STRtrees for the two grids, to speed up intersection queries
    # we may want to separately tune nodecapacity if one is much larger than the other.  
    # specifically we may want to tune leaf node capacity via Hilbert packing while still 
    # constraining inner node capacity.  But that can come later.
    tree1 = STRtree(grid1; nodecapacity = nodecapacity1) 
    tree2 = STRtree(grid2; nodecapacity = nodecapacity2)
    # Do the dual query, which is the most efficient way to do this,
    # by iterating down both trees simultaneously, rejecting pairs of nodes that do not intersect.
    # when we find an intersection, we calculate the area of the intersection and add it to the result matrix.
    GO.SpatialTreeInterface.do_dual_query(Extents.intersects, tree1, tree2) do i1, i2
        p1, p2 = grid1[i1], grid2[i2]
        # may want to check if the polygons intersect first, 
        # to avoid antimeridian-crossing multipolygons viewing a scanline.
        intersection_polys = try # can remove this now, got all the errors cleared up in the fix.
            # At some future point, we may want to add the manifold here
            # but for right now, GeometryOps only supports planar polygons anyway.
            GO.intersection(p1, p2; target = GI.PolygonTrait())
        catch e
            @error "Intersection failed!" i1 i2
            rethrow(e)
        end

        area_of_intersection = GO.area(m, intersection_polys)
        if area_of_intersection > 0
            A[i1, i2] += area_of_intersection
        end
    end

    return A
end

grid1 = begin
    gridpoints = [(i, j) for i in 0:2, j in 0:2]
    [GI.Polygon([GI.LinearRing([gridpoints[i, j], gridpoints[i, j+1], gridpoints[i+1, j+1], gridpoints[i+1, j], gridpoints[i, j]])]) for i in 1:size(gridpoints, 1)-1, j in 1:size(gridpoints, 2)-1] |> vec
end

grid2 = begin
    diamondpoly = GI.Polygon([GI.LinearRing([(0, 1), (1, 2), (2, 1), (1, 0), (0, 1)])])
    trianglepolys = GI.Polygon.([
        [GI.LinearRing([(0, 0), (1, 0), (0, 1), (0, 0)])],
        [GI.LinearRing([(0, 1), (0, 2), (1, 2), (0, 1)])],
        [GI.LinearRing([(1, 2), (2, 1), (2, 2), (1, 2)])],
        [GI.LinearRing([(2, 1), (2, 0), (1, 0), (2, 1)])],
    ])
    [diamondpoly, trianglepolys...]
end

A = area_of_intersection_operator(grid1, grid2)

# Now, let's perform some interpolation!
area1 = vec(sum(A, dims=2))
# test: @assert area1 == GO.area.(grid1)
area2 = vec(sum(A, dims=1))
# test: @assert area2 == GO.area.(grid2)

values_on_grid2 = [0, 0, 5, 0, 0]
poly(grid2; color = values_on_grid2, strokewidth = 2, strokecolor = :red)

values_on_grid1 = A * values_on_grid2 ./ area1
@assert sum(values_on_grid1 .* area1) == sum(values_on_grid2 .* area2)
poly(grid1; color = values_on_grid1, strokewidth = 2, strokecolor = :blue)

values_back_on_grid2 = A' * values_on_grid1 ./ area2
@assert sum(values_back_on_grid2 .* area2) == sum(values_on_grid2 .* area2)
poly(grid2; color = values_back_on_grid2, strokewidth = 2, strokecolor = :green)
# We can see here that some data has diffused into the central diamond cell of grid2,
# since it was overlapped by the top left cell of grid1.


using SpeedyWeather
using GeoMakie

SpeedyWeatherGeoMakieExt = Base.get_extension(SpeedyWeather, :SpeedyWeatherGeoMakieExt)

grid1 = rand(OctaHEALPixGrid, 5 + 100)
grid2 = rand(FullGaussianGrid, 4 + 100)

faces1 = SpeedyWeatherGeoMakieExt.get_faces(grid1)
faces2 = SpeedyWeatherGeoMakieExt.get_faces(grid2)

polys1 = GI.Polygon.(GI.LinearRing.(eachcol(faces1))) .|> GO.CutAtAntimeridianAndPoles() .|> GO.fix
polys2 = GI.Polygon.(GI.LinearRing.(eachcol(faces2))) .|> GO.CutAtAntimeridianAndPoles() .|> GO.fix

A = @time area_of_intersection_operator(polys1, polys2)

p1 = polys1[93]
p2 = polys2[105]

f, a, p = poly(p1)
poly!(a, p2)
f

# bug found in Foster Hormann tracing
# but geos also does the same thing
boxpoly = GI.Polygon([GI.LinearRing([(0, 0), (2, 0), (2, 2), (0, 2), (0, 0)])])
diamondpoly = GI.Polygon([GI.LinearRing([(0, 1), (1, 2), (2, 1), (1, 0), (0, 1)])])

diffpoly = GO.difference(boxpoly, diamondpoly; target = GI.PolygonTrait()) |> only
cutpolys = GO.cut(diffpoly, GI.Line([(0, 0), (4, 0)])) # even cut misbehaves!




