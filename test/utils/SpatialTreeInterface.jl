using Test
import GeometryOps as GO, GeoInterface as GI
using GeometryOps.SpatialTreeInterface
using GeometryOps.SpatialTreeInterface: isspatialtree, isleaf, getchild, nchild, child_indices_extents, node_extent
using GeometryOps.SpatialTreeInterface: query, depth_first_search, dual_depth_first_search
using GeometryOps.SpatialTreeInterface: FlatNoTree
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

# Test NaturalIndex implementation
@testset "STRtree" begin
    test_basic_interface(NaturalIndex)
    test_child_indices_extents(NaturalIndex)
    test_query_functionality(NaturalIndex)
    test_dual_query_functionality(NaturalIndex)
    test_geometry_support(NaturalIndex)
    test_find_point_in_all_countries(NaturalIndex)
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