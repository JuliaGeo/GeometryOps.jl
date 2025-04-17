using Test

using GeometryOps.SpatialTreeInterface
using GeometryOps.SpatialTreeInterface: isspatialtree, isleaf, getchild, nchild, child_indices_extents, node_extent
using GeometryOps.SpatialTreeInterface: query, depth_first_search, dual_depth_first_search
using GeometryOps.SpatialTreeInterface: FlatNoTree
using Extents
using GeoInterface: GeoInterface as GI
using SortTileRecursiveTree: STRtree

# Generic test functions for spatial trees
function test_basic_interface(TreeType)
    @testset "Basic interface" begin
        # Create a simple tree with one extent
        extents = [Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0))]
        tree = TreeType(extents)

        @test isspatialtree(tree)
        @test isspatialtree(typeof(tree))
    end
end

function test_child_indices_extents(TreeType)
    @testset "child_indices_extents" begin
        # Create a tree with three extents
        extents = [
            Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
            Extents.Extent(X=(1.0, 2.0), Y=(1.0, 2.0)),
            Extents.Extent(X=(2.0, 3.0), Y=(2.0, 3.0))
        ]
        tree = TreeType(extents)
        
        # Test that we get the correct indices and extents
        indices_extents = collect(child_indices_extents(tree))
        @test length(indices_extents) == 3
        @test indices_extents[1] == (1, extents[1])
        @test indices_extents[2] == (2, extents[2])
        @test indices_extents[3] == (3, extents[3])
    end
end

function test_query_functionality(TreeType)
    @testset "Query functionality" begin
        # Create a tree with three extents
        extents = [
            Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
            Extents.Extent(X=(1.0, 2.0), Y=(1.0, 2.0)),
            Extents.Extent(X=(2.0, 3.0), Y=(2.0, 3.0))
        ]
        tree = TreeType(extents)
        
        # Test query with a predicate that matches all
        all_pred = x -> true
        results = query(tree, all_pred)
        @test sort(results) == [1, 2, 3]

        # Test query with a predicate that matches none
        none_pred = x -> false
        results = query(tree, none_pred)
        @test isempty(results)

        # Test query with a specific extent predicate
        search_extent = Extents.Extent(X=(0.5, 1.5), Y=(0.5, 1.5))
        results = query(tree, Base.Fix1(Extents.intersects, search_extent))
        @test sort(results) == [1, 2]  # Should match first two extents
    end
end

function test_dual_query_functionality(TreeType)
    @testset "Dual query functionality" begin
        # Create two trees with overlapping extents
        tree1 = TreeType([
            Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
            Extents.Extent(X=(1.0, 2.0), Y=(1.0, 2.0))
        ])
        tree2 = TreeType([
            Extents.Extent(X=(0.5, 1.5), Y=(0.5, 1.5)),
            Extents.Extent(X=(1.5, 2.5), Y=(1.5, 2.5))
        ])

        # Test dual query with a predicate that matches all
        all_pred = (x, y) -> true
        results = Tuple{Int, Int}[]
        dual_depth_first_search((i, j) -> push!(results, (i, j)), all_pred, tree1, tree2)
        @test length(results) == 4  # 2 points in tree1 * 2 points in tree2

        # Test dual query with a specific predicate
        intersects_pred = (x, y) -> Extents.intersects(x, y)
        results = Tuple{Int, Int}[]
        dual_depth_first_search((i, j) -> push!(results, (i, j)), intersects_pred, tree1, tree2)
        @test sort(results) == [(1,1), (2,1), (2,2)]
    end
end

function test_geometry_support(TreeType)
    @testset "Geometry support" begin
        # Create a tree with 100 points
        points = tuple.(1:100, 1:100)
        tree = TreeType(points)
        
        # Test basic interface
        @test isspatialtree(tree)
        @test isspatialtree(typeof(tree))
        
        # Test query functionality
        all_pred = x -> true
        results = query(tree, all_pred)
        @test sort(results) == collect(1:100)

        none_pred = x -> false
        results = query(tree, none_pred)
        @test isempty(results)

        search_extent = Extents.Extent(X=(45.0, 55.0), Y=(45.0, 55.0))
        results = query(tree, Base.Fix1(Extents.intersects, search_extent))
        @test sort(results) == collect(45:55)
    end
end

# Test FlatNoTree implementation
@testset "FlatNoTree" begin
    test_basic_interface(FlatNoTree)
    test_child_indices_extents(FlatNoTree)
    test_query_functionality(FlatNoTree)
    test_dual_query_functionality(FlatNoTree)
    test_geometry_support(FlatNoTree)
end

# Test STRtree implementation
@testset "STRtree" begin
    test_basic_interface(STRtree)
    test_child_indices_extents(STRtree)
    test_query_functionality(STRtree)
    test_dual_query_functionality(STRtree)
    test_geometry_support(STRtree)
end



