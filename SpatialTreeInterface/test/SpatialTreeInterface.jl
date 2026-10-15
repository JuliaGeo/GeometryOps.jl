using Test
import GeoInterface as GI
using SpatialTreeInterface
using SpatialTreeInterface: isspatialtree, isleaf, getchild, nchild, child_indices_extents, node_extent
using SpatialTreeInterface: query, depth_first_search, dual_depth_first_search
using SpatialTreeInterface: FlatNoTree, spatialtree, node_extent_is_expensive
using SpatialTreeInterface.LoopStateMachine: Action
using Extents
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

        # An extent or a geometry is sanitized into an intersection predicate
        @test query(tree, search_extent) == [1, 2]
        @test query(tree, (0.5, 0.5)) == [1]
    end
end

function test_dual_query_functionality(TreeType)
    @testset "Dual query functionality - simple" begin
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
    @testset "Dual tree query with many boundingboxes" begin
        xs = 1:100
        ys = 1:100
        extent_grid = [Extents.Extent(X=(x, x+1), Y=(y, y+1)) for x in xs, y in ys] |> vec
        point_grid = [(x + 0.5, y + 0.5) for x in xs, y in ys] |> vec

        extent_tree = TreeType(extent_grid)
        point_tree = TreeType(point_grid)

        found_everything = falses(length(extent_grid))
        dual_depth_first_search(Extents.intersects, extent_tree, point_tree) do i, j
            if i == j
                found_everything[i] = true
            end
        end
        @test all(found_everything)
    end

    @testset "Imbalanced dual query - tree 1 deeper than tree 2" begin
        xs = 0:0.01:10
        ys = 0:0.01:10
        extent_grid = [Extents.Extent(X=(x, x+0.1), Y=(y, y+0.1)) for x in xs, y in ys] |> vec
        point_grid = [(x + 0.5, y + 0.5) for x in 0:9, y in 0:9] |> vec

        extent_tree = TreeType(extent_grid)
        point_tree = TreeType(point_grid)

        found_everything = falses(length(point_grid))
        dual_depth_first_search(Extents.intersects, extent_tree, point_tree) do i, j
            if Extents.intersects(extent_grid[i], GI.extent(point_grid[j]))
                found_everything[j] = true
            end
        end
        @test all(found_everything)
    end

    @testset "Imbalanced dual query - tree 2 deeper than tree 1" begin
        xs = 0:0.01:10
        ys = 0:0.01:10
        extent_grid = [Extents.Extent(X=(x, x+0.1), Y=(y, y+0.1)) for x in xs, y in ys] |> vec
        point_grid = [(x + 0.5, y + 0.5) for x in 0:9, y in 0:9] |> vec

        extent_tree = TreeType(extent_grid)
        point_tree = TreeType(point_grid)

        found_everything = falses(length(point_grid))
        dual_depth_first_search(Extents.intersects, point_tree, extent_tree) do i, j
            if Extents.intersects(GI.extent(point_grid[i]), extent_grid[j])
                found_everything[i] = true
            end
        end
        @test all(found_everything)
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

@testset "Interface fallbacks" begin
    @test !isspatialtree(1)
    @test !isspatialtree(Int)
    @test !node_extent_is_expensive(1)
    @test_throws ErrorException isleaf(1)
    @test_throws ErrorException nchild(1)
    @test_throws ErrorException getchild(1, 1)
    @test length(methods(spatialtree)) == 0
end

# Test FlatNoTree implementation
@testset "FlatNoTree" begin
    test_basic_interface(FlatNoTree)
    test_child_indices_extents(FlatNoTree)
    test_query_functionality(FlatNoTree)
    test_dual_query_functionality(FlatNoTree)
    test_geometry_support(FlatNoTree)
end

# Test STRtree implementation (via the SortTileRecursiveTree extension)
@testset "STRtree" begin
    @test !isnothing(Base.get_extension(SpatialTreeInterface, :SpatialTreeInterfaceSortTileRecursiveTreeExt))
    test_basic_interface(STRtree)
    test_child_indices_extents(STRtree)
    test_query_functionality(STRtree)
    test_dual_query_functionality(STRtree)
    test_geometry_support(STRtree)
end

# Wraps another tree, derives its extents on access, and counts the derivations -
# which is what `node_extent_is_expensive` is there to reduce.
const NE_CALLS = Ref(0)

struct CountedExtentNode{N}
    node::N
end
SpatialTreeInterface.isspatialtree(::Type{<:CountedExtentNode}) = true
SpatialTreeInterface.node_extent_is_expensive(::Type{<:CountedExtentNode}) = true
SpatialTreeInterface.isleaf(n::CountedExtentNode) = isleaf(n.node)
SpatialTreeInterface.nchild(n::CountedExtentNode) = nchild(n.node)
SpatialTreeInterface.getchild(n::CountedExtentNode) = (CountedExtentNode(c) for c in getchild(n.node))
SpatialTreeInterface.getchild(n::CountedExtentNode, i) = CountedExtentNode(getchild(n.node, i))
SpatialTreeInterface.child_indices_extents(n::CountedExtentNode) = child_indices_extents(n.node)
function SpatialTreeInterface.node_extent(n::CountedExtentNode)
    NE_CALLS[] += 1
    return node_extent(n.node)
end

# Same, but not opted in, so the two can be compared directly.
struct PlainCountedNode{N}
    node::N
end
SpatialTreeInterface.isspatialtree(::Type{<:PlainCountedNode}) = true
SpatialTreeInterface.isleaf(n::PlainCountedNode) = isleaf(n.node)
SpatialTreeInterface.nchild(n::PlainCountedNode) = nchild(n.node)
SpatialTreeInterface.getchild(n::PlainCountedNode) = (PlainCountedNode(c) for c in getchild(n.node))
SpatialTreeInterface.getchild(n::PlainCountedNode, i) = PlainCountedNode(getchild(n.node, i))
SpatialTreeInterface.child_indices_extents(n::PlainCountedNode) = child_indices_extents(n.node)
function SpatialTreeInterface.node_extent(n::PlainCountedNode)
    NE_CALLS[] += 1
    return node_extent(n.node)
end

const XYExtent = Extents.Extent{(:X,:Y),Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}}

# A tree whose nodes have differing numbers of children, held in a tuple, so
# `getchild` is inferred as `Tuple{Vararg{VaryingFanoutNode}}`.  A lazy
# `(child, extent)` iterator over that is a `Generator{I} where I`, which boxes;
# the traversal must not build one.
struct VaryingFanoutNode
    children::Tuple{Vararg{VaryingFanoutNode}}
    items::Vector{Tuple{Int,XYExtent}}
    extent::XYExtent
end
SpatialTreeInterface.isspatialtree(::Type{VaryingFanoutNode}) = true
SpatialTreeInterface.isleaf(n::VaryingFanoutNode) = isempty(n.children)
SpatialTreeInterface.nchild(n::VaryingFanoutNode) = length(n.children)
SpatialTreeInterface.getchild(n::VaryingFanoutNode) = n.children
SpatialTreeInterface.getchild(n::VaryingFanoutNode, i) = n.children[i]
SpatialTreeInterface.node_extent(n::VaryingFanoutNode) = n.extent
SpatialTreeInterface.child_indices_extents(n::VaryingFanoutNode) = n.items

function varying_fanout_tree(items::Vector{Tuple{Int,XYExtent}}, depth = 0)
    length(items) <= 4 && return VaryingFanoutNode((), items, reduce(Extents.union, last.(items)))
    step = cld(length(items), 2 + mod(depth, 4))  # fanout cycles 2, 3, 4, 5 by level
    kids = Tuple(varying_fanout_tree(items[i:min(i+step-1, end)], depth + 1)
                 for i in 1:step:length(items))
    return VaryingFanoutNode(kids, Tuple{Int,XYExtent}[],
                             reduce(Extents.union, [c.extent for c in kids]))
end
varying_fanout_tree(extents::Vector{XYExtent}) = varying_fanout_tree(collect(enumerate(extents)))

collect_pairs(args...) = (out = Tuple{Int,Int}[];
                          dual_depth_first_search((i, j) -> push!(out, (i, j)), args...);
                          sort!(out))

@testset "dual_depth_first_search extent carrying" begin
    # deliberately mismatched depths - that is where extents are carried through
    # a one-sided descent, and where the 4-arg/6-arg forms could diverge
    fine = [Extents.Extent(X=(x, x+0.1), Y=(y, y+0.1)) for x in 0:0.05:5, y in 0:0.05:5] |> vec
    coarse = [Extents.Extent(X=(x, x+1.0), Y=(y, y+1.0)) for x in 0:5, y in 0:5] |> vec

    @testset "4-arg and 6-arg forms agree ($(nameof(TreeType)))" for TreeType in
            (STRtree, FlatNoTree, varying_fanout_tree)
        t1 = TreeType(fine)
        t2 = TreeType(coarse)
        four = collect_pairs(Extents.intersects, t1, t2)
        six = (out = Tuple{Int,Int}[];
               dual_depth_first_search((i, j) -> push!(out, (i, j)), Extents.intersects,
                                       t1, node_extent(t1), t2, node_extent(t2));
               sort!(out))
        @test four == six
        @test !isempty(four)
        # and the pairs really are the intersecting ones
        @test all(((i, j),) -> Extents.intersects(fine[i], coarse[j]), four)
    end

    @testset "node_extent_is_expensive changes cost, not results" begin
        base = STRtree(fine)
        other = STRtree(coarse)

        NE_CALLS[] = 0
        plain = collect_pairs(Extents.intersects, PlainCountedNode(base), PlainCountedNode(other))
        plain_calls = NE_CALLS[]

        NE_CALLS[] = 0
        opted = collect_pairs(Extents.intersects, CountedExtentNode(base), CountedExtentNode(other))
        opted_calls = NE_CALLS[]

        # identical answers...
        @test plain == opted
        @test plain == collect_pairs(Extents.intersects, base, other)
        # ...but the opted-in tree derives each extent far fewer times
        @test opted_calls < plain_calls
    end

    @testset "default trait allocates nothing (varying fanout)" begin
        t1 = varying_fanout_tree(fine)
        t2 = varying_fanout_tree(coarse)
        @test node_extent_is_expensive(t1) == false
        noop(i, j) = nothing
        dual_depth_first_search(noop, Extents.intersects, t1, t2)  # compile
        @test (@allocated dual_depth_first_search(noop, Extents.intersects, t1, t2)) == 0
    end

    # `_child_extents` picks the caching strategy; it has to fold to a concrete
    # type from the node type alone, or the traversal boxes its loop variables.
    @testset "the caching strategy is a compile-time choice" begin
        cheap = varying_fanout_tree(fine)
        costly = CountedExtentNode(STRtree(fine))
        @test (@inferred SpatialTreeInterface._child_extents(cheap)) === nothing
        @test (@inferred SpatialTreeInterface._child_extents(STRtree(fine))) === nothing
        @test (@inferred SpatialTreeInterface._child_extents(costly)) isa Vector{<:Extents.Extent}
    end
end

@testset "dual_depth_first_search full_return" begin
    # Each of the four branches of the recursion has to propagate an
    # `Action(:full_return, x)` back out through every enclosing level.
    fine = [Extents.Extent(X=(x, x+0.1), Y=(y, y+0.1)) for x in 0:0.05:5, y in 0:0.05:5] |> vec
    coarse = [Extents.Extent(X=(x, x+1.0), Y=(y, y+1.0)) for x in 0:5, y in 0:5] |> vec
    single = [Extents.Extent(X=(0.0, 6.0), Y=(0.0, 6.0))]

    # deep-vs-deep (both-internal branch), deep-vs-shallow and shallow-vs-deep
    # (the one-leaf branches), and leaf-vs-leaf
    @testset "$label" for (label, a, b) in (
        ("both internal", STRtree(fine), STRtree(coarse)),
        ("leaf vs internal", STRtree(single), STRtree(fine)),
        ("internal vs leaf", STRtree(fine), STRtree(single)),
        ("leaf vs leaf", STRtree(single), STRtree(single)),
    )
        n = Ref(0)
        result = dual_depth_first_search(Extents.intersects, a, b) do i, j
            n[] += 1
            return Action(:full_return, (i, j))
        end
        # returns the Action itself, and stops at the very first pair
        @test result isa Action
        @test result.name == :full_return
        @test n[] == 1

        # the same thing through the 6-arg form
        m = Ref(0)
        result6 = dual_depth_first_search(Extents.intersects, a, node_extent(a), b, node_extent(b)) do i, j
            m[] += 1
            return Action(:full_return, (i, j))
        end
        @test result6 isa Action
        @test result6.name == :full_return
        @test m[] == 1
        @test result6.x == result.x
    end

    @testset "full_return with an opted-in tree" begin
        a = CountedExtentNode(STRtree(fine))
        b = CountedExtentNode(STRtree(coarse))
        n = Ref(0)
        result = dual_depth_first_search(Extents.intersects, a, b) do i, j
            n[] += 1
            return Action(:full_return, i)
        end
        @test result isa Action
        @test result.name == :full_return
        @test n[] == 1
    end
end

@testset "depth_first_search full_return" begin
    fine = [Extents.Extent(X=(x, x+0.1), Y=(y, y+0.1)) for x in 0:0.05:5, y in 0:0.05:5] |> vec
    tree = STRtree(fine)
    n = Ref(0)
    result = depth_first_search(Base.Fix1(Extents.intersects, Extents.Extent(X=(0, 6), Y=(0, 6))), tree) do i
        n[] += 1
        return Action(:full_return, i)
    end
    @test result isa Action
    @test n[] == 1
    # the two-argument form collects the matching indices
    @test sort!(depth_first_search(Base.Fix1(Extents.intersects, Extents.Extent(X=(0, 6), Y=(0, 6))), tree)) == 1:length(fine)
end
