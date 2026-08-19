using Test
import GeometryOps as GO, GeoInterface as GI
using GeometryOps.SpatialTreeInterface
using GeometryOps.SpatialTreeInterface: isspatialtree, isleaf, getchild, nchild, child_indices_extents, node_extent
using GeometryOps.SpatialTreeInterface: query, depth_first_search, dual_depth_first_search
using GeometryOps.SpatialTreeInterface: FlatNoTree, spatialtree, node_extent_is_expensive, children_extent_type
using GeometryOps.LoopStateMachine: Action
using GeometryOps.FlexibleRTrees: RTree, STR
using GeometryOps.NaturalIndexing: NaturalIndex
using Extents
using SortTileRecursiveTree: STRtree
using NaturalEarth
using Polylabel

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

function test_find_point_in_all_countries(TreeType)
    all_countries = NaturalEarth.naturalearth("admin_0_countries", 10)
    tree = TreeType(all_countries.geometry)

    ber = (13.4050, 52.5200)   # Berlin
    nyc = (-74.0060, 40.7128)  # New York City
    sin = (103.8198, 1.3521)   # Singapore

    @testset "locate points using query" begin
        @testset let point = ber, name = "Berlin"
            # Test Berlin (should be in Germany)
            results = query(tree, point)
            @test any(i -> all_countries.ADM0_A3[i] == "DEU", results)
        end
        @testset let point = nyc, name = "New York City"
            # Test NYC (should be in USA)
            results = query(tree, point)
            @test any(i -> all_countries.ADM0_A3[i] == "USA", results)
        end
        @testset let point = sin, name = "Singapore"
            # Test Singapore
            results = query(tree, point)
            @test any(i -> all_countries.ADM0_A3[i] == "SGP", results)
        end
    end
end

# Test FlatNoTree implementation
@testset "FlatNoTree" begin
    test_basic_interface(FlatNoTree)
    test_child_indices_extents(FlatNoTree)
    test_query_functionality(FlatNoTree)
    test_dual_query_functionality(FlatNoTree)
    test_geometry_support(FlatNoTree)
    test_find_point_in_all_countries(FlatNoTree)
end

# Test STRtree implementation
@testset "STRtree" begin
    test_basic_interface(STRtree)
    test_child_indices_extents(STRtree)
    test_query_functionality(STRtree)
    test_dual_query_functionality(STRtree)
    test_geometry_support(STRtree)
    test_find_point_in_all_countries(STRtree)
end

@testset "spatialtree" begin
    geometries = [missing, nothing, (0.0, 0.0), nothing, (1.0, 1.0)]
    tree = spatialtree(geometries)

    @test tree isa RTree{STR}
    @test query(tree, Extents.Extent(X = (-0.1, 1.1), Y = (-0.1, 1.1))) == [3, 5]
    @test isnothing(spatialtree([missing, nothing, missing]))

    # the default is `Planar()`
    @test query(spatialtree(GO.Planar(), geometries),
        Extents.Extent(X = (-0.1, 1.1), Y = (-0.1, 1.1))) == [3, 5]

    # `Spherical()` indexes the 3D regions on the unit sphere: the cap encloses
    # the pole none of its vertices reach, and the skipped entries still do not
    # shift the reported indices
    cap = GI.Polygon([[(lon, 80.0) for lon in 0.0:30.0:360.0]])
    spherical = spatialtree(GO.Spherical(), [missing, nothing, cap])
    @test Extents.extent(spherical) isa Extents.Extent{(:X, :Y, :Z)}
    @test query(spherical, Extents.Extent(X = (-0.01, 0.01), Y = (-0.01, 0.01), Z = (0.99, 1.0))) == [3]
end

# Test NaturalIndex implementation
@testset "STRtree" begin
    test_basic_interface(NaturalIndex)
    test_child_indices_extents(NaturalIndex)
    test_query_functionality(NaturalIndex)
    test_dual_query_functionality(NaturalIndex)
    test_geometry_support(NaturalIndex)
    test_find_point_in_all_countries(NaturalIndex)
end

# Wraps another tree, derives its extents on access, and counts the derivations -
# which is what `node_extent_is_expensive` is there to reduce.
const NE_CALLS = Ref(0)

struct CountedExtentNode{N}
    node::N
end
GO.SpatialTreeInterface.isspatialtree(::Type{<:CountedExtentNode}) = true
GO.SpatialTreeInterface.node_extent_is_expensive(::Type{<:CountedExtentNode}) = true
GO.SpatialTreeInterface.isleaf(n::CountedExtentNode) = isleaf(n.node)
GO.SpatialTreeInterface.nchild(n::CountedExtentNode) = nchild(n.node)
GO.SpatialTreeInterface.getchild(n::CountedExtentNode) = (CountedExtentNode(c) for c in getchild(n.node))
GO.SpatialTreeInterface.getchild(n::CountedExtentNode, i) = CountedExtentNode(getchild(n.node, i))
GO.SpatialTreeInterface.child_indices_extents(n::CountedExtentNode) = child_indices_extents(n.node)
function GO.SpatialTreeInterface.node_extent(n::CountedExtentNode)
    NE_CALLS[] += 1
    return node_extent(n.node)
end

# Same, but not opted in, so the two can be compared directly.
struct PlainCountedNode{N}
    node::N
end
GO.SpatialTreeInterface.isspatialtree(::Type{<:PlainCountedNode}) = true
GO.SpatialTreeInterface.isleaf(n::PlainCountedNode) = isleaf(n.node)
GO.SpatialTreeInterface.nchild(n::PlainCountedNode) = nchild(n.node)
GO.SpatialTreeInterface.getchild(n::PlainCountedNode) = (PlainCountedNode(c) for c in getchild(n.node))
GO.SpatialTreeInterface.getchild(n::PlainCountedNode, i) = PlainCountedNode(getchild(n.node, i))
GO.SpatialTreeInterface.child_indices_extents(n::PlainCountedNode) = child_indices_extents(n.node)
function GO.SpatialTreeInterface.node_extent(n::PlainCountedNode)
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
GO.SpatialTreeInterface.isspatialtree(::Type{VaryingFanoutNode}) = true
GO.SpatialTreeInterface.isleaf(n::VaryingFanoutNode) = isempty(n.children)
GO.SpatialTreeInterface.nchild(n::VaryingFanoutNode) = length(n.children)
GO.SpatialTreeInterface.getchild(n::VaryingFanoutNode) = n.children
GO.SpatialTreeInterface.getchild(n::VaryingFanoutNode, i) = n.children[i]
GO.SpatialTreeInterface.node_extent(n::VaryingFanoutNode) = n.extent
GO.SpatialTreeInterface.child_indices_extents(n::VaryingFanoutNode) = n.items

function varying_fanout_tree(items::Vector{Tuple{Int,XYExtent}}, depth = 0)
    length(items) <= 4 && return VaryingFanoutNode((), items, reduce(Extents.union, last.(items)))
    step = cld(length(items), 2 + mod(depth, 4))  # fanout cycles 2, 3, 4, 5 by level
    kids = Tuple(varying_fanout_tree(items[i:min(i+step-1, end)], depth + 1)
                 for i in 1:step:length(items))
    return VaryingFanoutNode(kids, Tuple{Int,XYExtent}[],
                             reduce(Extents.union, [c.extent for c in kids]))
end
varying_fanout_tree(extents::Vector{XYExtent}) = varying_fanout_tree(collect(enumerate(extents)))

# An extent-like type that is not an `Extents.Extent`, so a tree using both has
# an extent type that genuinely changes with depth.
struct BoxExtent
    X::Tuple{Float64,Float64}
    Y::Tuple{Float64,Float64}
end
_asbox(e) = BoxExtent(e.X, e.Y)
_boxy_intersects(a, b) = Extents.intersects(
    Extents.Extent(X = a.X, Y = a.Y), Extents.Extent(X = b.X, Y = b.Y))

# Odd levels report an `Extent`, even levels a `BoxExtent`.  Siblings agree,
# levels do not - the case one tree-wide stack cannot serve.
struct AlternatingExtentNode{N}
    node::N
    level::Int
end
AlternatingExtentNode(node) = AlternatingExtentNode(node, 0)
GO.SpatialTreeInterface.isspatialtree(::Type{<:AlternatingExtentNode}) = true
GO.SpatialTreeInterface.node_extent_is_expensive(::Type{<:AlternatingExtentNode}) = true
GO.SpatialTreeInterface.children_extent_type(::Type{<:AlternatingExtentNode}) = Union{XYExtent,BoxExtent}
GO.SpatialTreeInterface.isleaf(n::AlternatingExtentNode) = isleaf(n.node)
GO.SpatialTreeInterface.nchild(n::AlternatingExtentNode) = nchild(n.node)
GO.SpatialTreeInterface.getchild(n::AlternatingExtentNode) =
    (AlternatingExtentNode(c, n.level + 1) for c in getchild(n.node))
GO.SpatialTreeInterface.getchild(n::AlternatingExtentNode, i) =
    AlternatingExtentNode(getchild(n.node, i), n.level + 1)
GO.SpatialTreeInterface.child_indices_extents(n::AlternatingExtentNode) = child_indices_extents(n.node)
GO.SpatialTreeInterface.node_extent(n::AlternatingExtentNode) =
    iseven(n.level) ? node_extent(n.node) : _asbox(node_extent(n.node))

# Same idea, but the trait is spelled on the node rather than on its type; both
# spellings have to work, and both have to stay inferrable.
struct InstanceTraitNode{N}
    node::N
end
GO.SpatialTreeInterface.isspatialtree(::Type{<:InstanceTraitNode}) = true
GO.SpatialTreeInterface.node_extent_is_expensive(::Type{<:InstanceTraitNode}) = true
GO.SpatialTreeInterface.children_extent_type(n::InstanceTraitNode) = XYExtent
GO.SpatialTreeInterface.isleaf(n::InstanceTraitNode) = isleaf(n.node)
GO.SpatialTreeInterface.nchild(n::InstanceTraitNode) = nchild(n.node)
GO.SpatialTreeInterface.getchild(n::InstanceTraitNode) = (InstanceTraitNode(c) for c in getchild(n.node))
GO.SpatialTreeInterface.getchild(n::InstanceTraitNode, i) = InstanceTraitNode(getchild(n.node, i))
GO.SpatialTreeInterface.child_indices_extents(n::InstanceTraitNode) = child_indices_extents(n.node)
GO.SpatialTreeInterface.node_extent(n::InstanceTraitNode) = node_extent(n.node)

collect_pairs(args...) = (out = Tuple{Int,Int}[];
                          dual_depth_first_search((i, j) -> push!(out, (i, j)), args...);
                          sort!(out))

@testset "dual_depth_first_search extent carrying" begin
    # deliberately mismatched depths - that is where extents are carried through
    # a one-sided descent, and where the 4-arg/6-arg forms could diverge
    fine = [Extents.Extent(X=(x, x+0.1), Y=(y, y+0.1)) for x in 0:0.05:5, y in 0:0.05:5] |> vec
    coarse = [Extents.Extent(X=(x, x+1.0), Y=(y, y+1.0)) for x in 0:5, y in 0:5] |> vec

    @testset "4-arg and 6-arg forms agree ($(nameof(TreeType)))" for TreeType in
            (STRtree, NaturalIndex, FlatNoTree, varying_fanout_tree)
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
        base = NaturalIndex(fine)
        other = NaturalIndex(coarse)

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

    @testset "default trait allocates nothing ($label)" for (label, t1, t2) in (
        ("fixed fanout", NaturalIndex(fine), NaturalIndex(coarse)),
        ("varying fanout", varying_fanout_tree(fine), varying_fanout_tree(coarse)),
    )
        @test node_extent_is_expensive(t1) == false
        noop(i, j) = nothing
        dual_depth_first_search(noop, Extents.intersects, t1, t2)  # compile
        @test (@allocated dual_depth_first_search(noop, Extents.intersects, t1, t2)) == 0
    end

    # `_extent_stack` picks the caching strategy; it has to fold to a concrete
    # type from the node type alone, or the traversal boxes its loop variables.
    @testset "the caching strategy is a compile-time choice" begin
        cheap = varying_fanout_tree(fine)
        costly = CountedExtentNode(NaturalIndex(fine))
        ext = node_extent(NaturalIndex(fine))
        @test (@inferred GO.SpatialTreeInterface._extent_stack(nothing, cheap, ext)) === nothing
        @test (@inferred GO.SpatialTreeInterface._extent_stack(nothing, NaturalIndex(fine), ext)) === nothing
        @test (@inferred GO.SpatialTreeInterface._extent_stack(nothing, costly, ext)) isa Vector{<:Extents.Extent}
        # an existing stack is reused rather than replaced
        stack = typeof(ext)[]
        @test (@inferred GO.SpatialTreeInterface._extent_stack(stack, costly, ext)) === stack
    end

    # the opted-in path used to allocate a vector per visited node pair; now it
    # allocates one scratch stack for the whole descent
    @testset "an opted-in tree allocates a bounded amount" begin
        a = CountedExtentNode(NaturalIndex(fine))
        b = CountedExtentNode(NaturalIndex(coarse))
        noop(i, j) = nothing
        dual_depth_first_search(noop, Extents.intersects, a, b)  # compile
        @test (@allocated dual_depth_first_search(noop, Extents.intersects, a, b)) < 4096
    end

    # the traversal asks the node, so a tree may answer on the node or on its type
    @testset "children_extent_type is asked of the node, in either spelling" begin
        onnode = InstanceTraitNode(NaturalIndex(fine))
        ontype = AlternatingExtentNode(NaturalIndex(fine))
        @test children_extent_type(onnode) === XYExtent                    # defined on the node
        @test children_extent_type(ontype) === Union{XYExtent,BoxExtent}   # defined on the type
        @test children_extent_type(NaturalIndex(fine)) === nothing         # says nothing, so the default
        # both spellings fold, so the stack type is still a compile-time choice
        @test (@inferred GO.SpatialTreeInterface._extent_stack(nothing, onnode, node_extent(onnode))) isa
              Vector{XYExtent}
        @test (@inferred GO.SpatialTreeInterface._extent_stack(nothing, ontype, node_extent(ontype))) isa
              Vector{Union{XYExtent,BoxExtent}}
        # ...and the whole traversal stays inferrable through each
        noop(i, j) = nothing
        @test (@inferred dual_depth_first_search(noop, Extents.intersects, onnode,
                                                 InstanceTraitNode(NaturalIndex(coarse)))) === nothing
        @test (@inferred dual_depth_first_search(noop, _boxy_intersects, ontype,
                                                 AlternatingExtentNode(NaturalIndex(coarse)))) === nothing
        # the node-spelled tree answers exactly as the plain one does
        @test collect_pairs(Extents.intersects, onnode, InstanceTraitNode(NaturalIndex(coarse))) ==
              collect_pairs(Extents.intersects, NaturalIndex(fine), NaturalIndex(coarse))
    end

    # a tree whose extent type changes with depth: siblings agree, levels do not
    @testset "children_extent_type carries a per-level extent type" begin
        a = AlternatingExtentNode(NaturalIndex(fine))
        b = AlternatingExtentNode(NaturalIndex(coarse))
        @test node_extent_is_expensive(a)
        @test children_extent_type(a) === Union{XYExtent,BoxExtent}
        # the answers are the ones the plain tree gives
        @test collect_pairs(_boxy_intersects, a, b) ==
              collect_pairs(Extents.intersects, NaturalIndex(fine), NaturalIndex(coarse))
        # and the stack the traversal picks is the declared child type
        @test (@inferred GO.SpatialTreeInterface._extent_stack(nothing, a, node_extent(a))) isa
              Vector{Union{XYExtent,BoxExtent}}
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
        ("both internal", NaturalIndex(fine), NaturalIndex(coarse)),
        ("leaf vs internal", NaturalIndex(single), NaturalIndex(fine)),
        ("internal vs leaf", NaturalIndex(fine), NaturalIndex(single)),
        ("leaf vs leaf", NaturalIndex(single), NaturalIndex(single)),
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
        a = CountedExtentNode(NaturalIndex(fine))
        b = CountedExtentNode(NaturalIndex(coarse))
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

# This testset is not used because Polylabel.jl has some issues.

#=


    @testset "Dual query functionality - every country's polylabel against every country" begin

        # Note that this is a perfectly balanced tree query - we don't yet have a test for unbalanced
        # trees (but could easily add one, e.g. by getting polylabels of admin-1 or admin-2 regions)
        # from Natural Earth, or by using GADM across many countries.

        all_countries = NaturalEarth.naturalearth("admin_0_countries", 10)
        all_adm0_a3 = all_countries.ADM0_A3
        all_geoms = all_countries.geometry
        # US minor outlying islands - bug in Polylabel.jl
        # A lot of small geoms have this issue, that there will be an error from the queue
        # because the cell exists in the queue already.
        # Not sure what to do about it, I don't want to check containment every time...
        deleteat!(all_adm0_a3, 205)
        deleteat!(all_geoms, 205)

        geom_tree = TreeType(all_geoms)

        polylabels = [Polylabel.polylabel(geom; rtol = 0.019) for geom in all_geoms]
        polylabel_tree = TreeType(polylabels)

        found_countries = falses(length(polylabels))

        dual_depth_first_search(Extents.intersects, geom_tree, polylabel_tree) do i, j
            if i == j
                found_countries[i] = true
            end
        end

        @test all(found_countries)
    end
=#