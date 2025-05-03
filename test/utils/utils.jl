using Test

import GeometryOps as GO, GeoInterface as GI
import Extents

point = GI.Point(1.0, 1.0)
linestring = GI.LineString([(1.0, 1.0), (2.0, 2.0)])
multilinestring = GI.MultiLineString([[(1.0, 1.0), (2.0, 2.0)], [(3.0, 3.0), (4.0, 4.0)]])

polygon = GI.Polygon([[(1.0, 1.0), (2.0, 2.0), (2.0, 1.0), (1.0, 1.0)]])
multipolygon = GI.MultiPolygon([[[(1.0, 1.0), (2.0, 2.0), (2.0, 1.0), (1.0, 1.0)]], [[(3.0, 3.0), (4.0, 4.0), (4.0, 3.0), (3.0, 3.0)]]])

# Test eachedge
@testset "eachedge" begin
    # Test LineString
    edges = collect(GO.eachedge(linestring))
    @test length(edges) == 1
    @test edges[1] == ((1.0, 1.0), (2.0, 2.0))

    # Test MultiLineString
    edges = collect(GO.eachedge(multilinestring))
    @test length(edges) == 2
    @test edges[1] == ((1.0, 1.0), (2.0, 2.0))
    @test edges[2] == ((3.0, 3.0), (4.0, 4.0))

    # Test Polygon
    edges = collect(GO.eachedge(polygon))
    @test length(edges) == 3
    @test edges[1] == ((1.0, 1.0), (2.0, 2.0))
    @test edges[2] == ((2.0, 2.0), (2.0, 1.0))
    @test edges[3] == ((2.0, 1.0), (1.0, 1.0))

    # Test MultiPolygon
    edges = collect(GO.eachedge(multipolygon))
    @test length(edges) == 6  # 3 edges per polygon
    @test edges[1] == ((1.0, 1.0), (2.0, 2.0))
    @test edges[2] == ((2.0, 2.0), (2.0, 1.0))
    @test edges[3] == ((2.0, 1.0), (1.0, 1.0))
    @test edges[4] == ((3.0, 3.0), (4.0, 4.0))
    @test edges[5] == ((4.0, 4.0), (4.0, 3.0))
    @test edges[6] == ((4.0, 3.0), (3.0, 3.0))

    # Test error cases
    @test_throws ArgumentError GO.eachedge(point)
    @test_throws ArgumentError GO.eachedge(GI.MultiPoint([(1.0, 1.0), (2.0, 2.0)]))
end

# Test to_edgelist
@testset "to_edgelist" begin
    # Test LineString
    edges = GO.to_edgelist(linestring)
    @test length(edges) == 1
    @test GI.getpoint(edges[1], 1) == (1.0, 1.0)
    @test GI.getpoint(edges[1], 2) == (2.0, 2.0)

    # Test MultiLineString
    edges = GO.to_edgelist(multilinestring)
    @test length(edges) == 2
    @test GI.getpoint(edges[1], 1) == (1.0, 1.0)
    @test GI.getpoint(edges[1], 2) == (2.0, 2.0)
    @test GI.getpoint(edges[2], 1) == (3.0, 3.0)
    @test GI.getpoint(edges[2], 2) == (4.0, 4.0)

    # Test Polygon
    edges = GO.to_edgelist(polygon)
    @test length(edges) == 3
    @test GI.getpoint(edges[1], 1) == (1.0, 1.0)
    @test GI.getpoint(edges[1], 2) == (2.0, 2.0)
    @test GI.getpoint(edges[2], 1) == (2.0, 2.0)
    @test GI.getpoint(edges[2], 2) == (2.0, 1.0)
    @test GI.getpoint(edges[3], 1) == (2.0, 1.0)
    @test GI.getpoint(edges[3], 2) == (1.0, 1.0)

    # Test MultiPolygon
    edges = GO.to_edgelist(multipolygon)
    @test length(edges) == 6

    # Test with extent filtering
    extent = Extents.Extent(X=(1.5, 2.5), Y=(1.5, 2.5))
    edges, indices = GO.to_edgelist(extent, linestring, Float64)
    @test length(edges) == 1
    @test indices == [1]

    # Test with extent that doesn't intersect
    extent = Extents.Extent(X=(3.0, 4.0), Y=(3.0, 4.0))
    edges, indices = GO.to_edgelist(extent, linestring, Float64)
    @test isempty(edges)
    @test isempty(indices)

    # Test error cases
    @test_throws ArgumentError GO.to_edgelist(point)
    @test_throws ArgumentError GO.to_edgelist(GI.MultiPoint([(1.0, 1.0), (2.0, 2.0)]))
end

# Test lazy_edgelist
@testset "lazy_edgelist" begin
    # Test LineString
    edges = collect(GO.lazy_edgelist(linestring))
    @test length(edges) == 1
    @test GI.getpoint(first(edges), 1) == (1.0, 1.0)
    @test GI.getpoint(first(edges), 2) == (2.0, 2.0)

    # Test MultiLineString
    edges = collect(GO.lazy_edgelist(multilinestring))
    @test length(edges) == 2
    @test edges == GO.to_edgelist(multilinestring)

    # Test Polygon
    edges = collect(GO.lazy_edgelist(polygon))
    @test length(edges) == 3
    @test edges == GO.to_edgelist(polygon)

    # Test MultiPolygon
    edges = collect(GO.lazy_edgelist(multipolygon))
    @test length(edges) == 6
    @test edges == GO.to_edgelist(multipolygon)

    # Test error cases
    @test_throws ArgumentError collect(GO.lazy_edgelist(point))
    @test_throws ArgumentError collect(GO.lazy_edgelist(GI.MultiPoint([(1.0, 1.0), (2.0, 2.0)])))
end

# Test edge_extents
@testset "edge_extents" begin
    # Test LineString
    extents = GO.edge_extents(linestring)
    @test length(extents) == 1
    @test extents[1].X == (1.0, 2.0)
    @test extents[1].Y == (1.0, 2.0)

    # Test MultiLineString
    extents = GO.edge_extents(multilinestring)
    @test length(extents) == 2
    @test extents[1].X == (1.0, 2.0)
    @test extents[1].Y == (1.0, 2.0)
    @test extents[2].X == (3.0, 4.0)
    @test extents[2].Y == (3.0, 4.0)

    # Test Polygon
    extents = GO.edge_extents(polygon)
    @test length(extents) == 3
    @test extents[1].X == (1.0, 2.0)
    @test extents[1].Y == (1.0, 2.0)
    @test extents[2].X == (2.0, 2.0)
    @test extents[2].Y == (1.0, 2.0)
    @test extents[3].X == (1.0, 2.0)
    @test extents[3].Y == (1.0, 1.0)

    # Test MultiPolygon
    extents = GO.edge_extents(multipolygon)
    @test length(extents) == 6
    @test extents[1].X == (1.0, 2.0)
    @test extents[1].Y == (1.0, 2.0)
    @test extents[4].X == (3.0, 4.0)
    @test extents[4].Y == (3.0, 4.0)

    # Test error cases
    @test_throws ArgumentError GO.edge_extents(point)
    @test_throws ArgumentError GO.edge_extents(GI.MultiPoint([(1.0, 1.0), (2.0, 2.0)]))
end

# Test lazy_edge_extents
@testset "lazy_edge_extents" begin
    # Test LineString
    extents = collect(GO.lazy_edge_extents(linestring))
    @test length(extents) == 1
    @test extents[1].X == (1.0, 2.0)
    @test extents[1].Y == (1.0, 2.0)

    # Test MultiLineString
    extents = collect(GO.lazy_edge_extents(multilinestring))
    @test length(extents) == 2
    @test extents[1].X == (1.0, 2.0)
    @test extents[1].Y == (1.0, 2.0)
    @test extents[2].X == (3.0, 4.0)
    @test extents[2].Y == (3.0, 4.0)

    # Test Polygon
    extents = collect(GO.lazy_edge_extents(polygon))
    @test length(extents) == 3
    @test extents == GO.edge_extents(polygon)

    # Test MultiPolygon
    extents = collect(GO.lazy_edge_extents(multipolygon))
    @test length(extents) == 6
    @test extents == GO.edge_extents(multipolygon)

    # Test error cases
    @test_throws ArgumentError collect(GO.lazy_edge_extents(point))
    @test_throws ArgumentError collect(GO.lazy_edge_extents(GI.MultiPoint([(1.0, 1.0), (2.0, 2.0)])))
end

