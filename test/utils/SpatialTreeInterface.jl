using Test

using GeometryOps.SpatialTreeInterface
using GeometryOps.SpatialTreeInterface: isspatialtree, isleaf, getchild, nchild, child_indices_extents, node_extent
using GeometryOps.SpatialTreeInterface: query, do_query, do_dual_query
using Extents
using GeoInterface: GeoInterface as GI

using GeometryOps.SpatialTreeInterface: FlatNoTree

# Test FlatNoTree implementation
@testset "FlatNoTree" begin
    # Test with a simple vector of extents
    extents = [
        Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
        Extents.Extent(X=(1.0, 2.0), Y=(1.0, 2.0)),
        Extents.Extent(X=(2.0, 3.0), Y=(2.0, 3.0))
    ]
    tree = FlatNoTree(extents)

    @testset "Basic interface" begin
        @test isleaf(tree)
        @test isspatialtree(tree)
        @test isspatialtree(typeof(tree))
    end

    @testset "child_indices_extents" begin
        # Test that we get the correct indices and extents
        indices_extents = collect(child_indices_extents(tree))
        @test length(indices_extents) == 3
        @test indices_extents[1] == (1, extents[1])
        @test indices_extents[2] == (2, extents[2])
        @test indices_extents[3] == (3, extents[3])
    end

    @testset "Query functionality" begin
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

    @testset "Dual query functionality" begin
        # Create two trees for dual query testing
        tree1 = FlatNoTree([
            Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
            Extents.Extent(X=(1.0, 2.0), Y=(1.0, 2.0))
        ])
        tree2 = FlatNoTree([
            Extents.Extent(X=(0.5, 1.5), Y=(0.5, 1.5)),
            Extents.Extent(X=(1.5, 2.5), Y=(1.5, 2.5))
        ])

        # Test dual query with a predicate that matches all
        all_pred = (x, y) -> true
        results = Tuple{Int, Int}[]
        do_dual_query((i, j) -> push!(results, (i, j)), all_pred, tree1, tree2)
        @test sort(results) == [(1,1), (1,2), (2,1), (2,2)]

        # Test dual query with a specific predicate
        intersects_pred = (x, y) -> Extents.intersects(x, y)
        results = Tuple{Int, Int}[]
        do_dual_query((i, j) -> push!(results, (i, j)), intersects_pred, tree1, tree2)
        @test sort(results) == [(1,1), (2,1), (2,2)]
    end
end

