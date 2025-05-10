# # Polygon clipping helpers
# This file contains the shared helper functions for the polygon clipping functionalities.

# This file specifically defines helpers for the Foster-Hormann clipping algorithm.


"""
    abstract type IntersectionAccelerator

The abstract supertype for all intersection accelerator types.

The idea is that these speed up the edge-edge intersection checking process,
perhaps at the cost of memory.

The naive case is `NestedLoop`, which is just a nested loop, running in O(n*m) time.

Then we have `SingleSTRtree`, which is a single STRtree, running in O(n*log(m)) time.

Then we have `DoubleSTRtree`, which is a simultaneous double-tree traversal of two STRtrees.

Finally, we have `AutoAccelerator`, which chooses the best
accelerator based on the size of the input polygons.  This gets materialized in `build_a_list` for now.
`AutoAccelerator` should also try to respect existing spatial indices, if they exist.
"""
abstract type IntersectionAccelerator end
struct NestedLoop <: IntersectionAccelerator end
struct SingleSTRtree <: IntersectionAccelerator end
struct DoubleSTRtree <: IntersectionAccelerator end
struct SingleNaturalTree <: IntersectionAccelerator end
struct DoubleNaturalTree <: IntersectionAccelerator end
struct ThinnedDoubleNaturalTree <: IntersectionAccelerator end

"""
    AutoAccelerator()

Let the algorithm choose the best accelerator based on the size of the input polygons.

Once we have prepared geometry, this will also consider the existing preparations on the geoms.
"""
struct AutoAccelerator <: IntersectionAccelerator end

"""
    FosterHormannClipping{M <: Manifold, A <: Union{Nothing, Accelerator}} <: GeometryOpsCore.Algorithm{M} 

Applies the Foster-Hormann clipping algorithm.

# Arguments
- `manifold::M`: The manifold on which the algorithm operates.
- `accelerator::A`: The accelerator to use for the algorithm.  Can be `nothing` for automatic choice, or a custom accelerator.
"""
struct FosterHormannClipping{M <: Manifold, A <: IntersectionAccelerator} <: GeometryOpsCore.Algorithm{M} 
    manifold::M
    accelerator::A
    # TODO: add exact flag
    # TODO: should exact flag be in the type domain?
end
FosterHormannClipping(; manifold::Manifold = Planar(), accelerator = nothing) = FosterHormannClipping(manifold, isnothing(accelerator) ? AutoAccelerator() : accelerator)
FosterHormannClipping(manifold::Manifold, accelerator::Union{Nothing, IntersectionAccelerator} = nothing) = FosterHormannClipping(manifold, isnothing(accelerator) ? AutoAccelerator() : accelerator)
FosterHormannClipping(accelerator::Union{Nothing, IntersectionAccelerator}) = FosterHormannClipping(Planar(), isnothing(accelerator) ? AutoAccelerator() : accelerator)
# special case for spherical / geodesic manifolds
# since they can't use STRtrees (because those don't work on the sphere)
FosterHormannClipping(manifold::Union{Spherical, Geodesic}, accelerator::Union{Nothing, IntersectionAccelerator} = nothing) = FosterHormannClipping(manifold, isnothing(accelerator) ? NestedLoop() : (accelerator isa AutoAccelerator ? NestedLoop() : accelerator))

# This enum defines which side of an edge a point is on
@enum PointEdgeSide left=1 right=2 unknown=3

# Constants assigned for readability
const enter, exit = true, false
const crossing, bouncing = true, false

#= A point can either be the start or end of an overlapping chain of points between two
polygons, or not an endpoint of a chain. =#
@enum EndPointType start_chain=1 end_chain=2 not_endpoint=3

#= This is the struct that makes up a_list and b_list. Many values are only used if point is
an intersection point (ipt). =#
@kwdef struct PolyNode{T <: AbstractFloat}
    point::Tuple{T,T}          # (x, y) values of given point
    inter::Bool = false        # If ipt, true, else 0
    neighbor::Int = 0          # If ipt, index of equivalent point in a_list or b_list, else 0
    idx::Int = 0               # If crossing point, index within sorted a_idx_list
    ent_exit::Bool = false     # If ipt, true if enter and false if exit, else false
    crossing::Bool = false     # If ipt, true if intersection crosses from out/in polygon, else false
    endpoint::EndPointType = not_endpoint # If ipt, denotes if point is the start or end of an overlapping chain
    fracs::Tuple{T,T} = (0., 0.) # If ipt, fractions along edges to ipt (a_frac, b_frac), else (0, 0)
end

#= Create a new node with all of the same field values as the given PolyNode unless
alternative values are provided, in which case those should be used. =#
PolyNode(node::PolyNode{T};
    point = node.point, inter = node.inter, neighbor = node.neighbor, idx = node.idx,
    ent_exit = node.ent_exit, crossing = node.crossing, endpoint = node.endpoint,
    fracs = node.fracs,
) where T = PolyNode{T}(;
    point = point, inter = inter, neighbor = neighbor, idx = idx, ent_exit = ent_exit,
    crossing = crossing, endpoint = endpoint, fracs = fracs)

# Checks equality of two PolyNodes by backing point value, fractional value, and intersection status
equals(pn1::PolyNode, pn2::PolyNode) = pn1.point == pn2.point && pn1.inter == pn2.inter && pn1.fracs == pn2.fracs
Base.:(==)(pn1::PolyNode, pn2::PolyNode) = equals(pn1, pn2)

# Finally, we define a nice error type for when the clipping tracing algorithm hits every point in a polygon.
# This stores the polygons, the a_list, and the b_list, and the a_idx_list.
# allowing the user to understand what happened and why.
"""
    TracingError{T1, T2} <: Exception

An error that is thrown when the clipping tracing algorithm fails somehow.
This is a bug in the algorithm, and should be reported.

The polygons are contained in the exception object, accessible by try-catch or as `err` in the REPL.
"""
struct TracingError{T1, T2, T} <: Exception
    message::String
    poly_a::T1
    poly_b::T2
    a_list::Vector{PolyNode{T}}
    b_list::Vector{PolyNode{T}}
    a_idx_list::Vector{Int}
end

function Base.showerror(io::IO, e::TracingError{T1, T2}) where {T1, T2}
    print(io, "TracingError: ")
    println(io, e.message)
    println(io, "Please open an issue with the polygons contained in this error object.")
    println(io)
    if max(GI.npoint(e.poly_a), GI.npoint(e.poly_b)) < 10
        println(io, "Polygon A:")
        println(io, GI.coordinates(e.poly_a))
        println(io)
        println(io, "Polygon B:")
        println(io, GI.coordinates(e.poly_b))
    else
        println(io, "The polygons are contained in the exception object, accessible by try-catch or as `err` in the REPL.")
    end
end



#=
    _build_ab_list(::Type{T}, poly_a, poly_b, delay_cross_f, delay_bounce_f; exact) ->
        (a_list, b_list, a_idx_list)

This function takes in two polygon rings and calls '_build_a_list', '_build_b_list', and
'_flag_ent_exit' in order to fully form a_list and b_list. The 'a_list' and 'b_list' that it
returns are the fully updated vectors of PolyNodes that represent the rings 'poly_a' and
'poly_b', respectively. This function also returns 'a_idx_list', which at its "ith" index
stores the index in 'a_list' at which the "ith" intersection point lies.
=#
function _build_ab_list(alg::FosterHormannClipping, ::Type{T}, poly_a, poly_b, delay_cross_f::F1, delay_bounce_f::F2; exact) where {T, F1, F2}
    # Make a list for nodes of each polygon
    a_list, a_idx_list, n_b_intrs = _build_a_list(alg, T, poly_a, poly_b; exact)
    b_list = _build_b_list(alg, T, a_idx_list, a_list, n_b_intrs, poly_b)

    # Flag crossings
    _classify_crossing!(alg, T, a_list, b_list; exact)

    # Flag the entry and exits
    _flag_ent_exit!(alg, T, GI.LinearRingTrait(), poly_b, a_list, delay_cross_f, Base.Fix2(delay_bounce_f, true); exact)
    _flag_ent_exit!(alg, T, GI.LinearRingTrait(), poly_a, b_list, delay_cross_f, Base.Fix2(delay_bounce_f, false); exact)

    # Set node indices and filter a_idx_list to just crossing points
    _index_crossing_intrs!(alg, a_list, b_list, a_idx_list)

    return a_list, b_list, a_idx_list
end


"The number of vertices past which we should use a STRtree for edge intersection checking."
const GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS = 32
# Fallback convenience method so we can just pass the algorithm in
function foreach_pair_of_maybe_intersecting_edges_in_order(
    alg::FosterHormannClipping{M, A}, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64
) where {FA, FAAfter, FI, T, M, A}
    return foreach_pair_of_maybe_intersecting_edges_in_order(alg.manifold, alg.accelerator, f_on_each_a, f_after_each_a, f_on_each_maybe_intersect, poly_a, poly_b, T)
end

"""
    foreach_pair_of_maybe_intersecting_edges_in_order(
        manifold::M, accelerator::A,
        f_on_each_a::FA,
        f_after_each_a::FAAfter,
        f_on_each_maybe_intersect::FI,
        geom_a,
        geom_b,
        ::Type{T} = Float64
    ) where {FA, FAAfter, FI, T, M <: Manifold, A <: IntersectionAccelerator}

Decompose `geom_a` and `geom_b` into edge lists (unsorted), and then, logically, 
perform the following iteration:

```julia
for (a_edge, i) in enumerate(eachedge(geom_a))
    f_on_each_a(a_edge, i)
    for (b_edge, j) in enumerate(eachedge(geom_b))
        if may_intersect(a_edge, b_edge)
            f_on_each_maybe_intersect(a_edge, b_edge)
        end
    end
    f_after_each_a(a_edge, i)
end
```

This may not be the exact acceleration that is performed - but it is 
the logical sequence of events.  It also uses the `accelerator`, 
and can automatically choose the best one based on an internal heuristic
if you pass in an [`AutoAccelerator`](@ref).  

For example, the `SingleSTRtree` accelerator is used along
with extent thinning to avoid unnecessary edge intersection 
checks in the inner loop.

"""
function foreach_pair_of_maybe_intersecting_edges_in_order(
    manifold::M, accelerator::A, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64
) where {FA, FAAfter, FI, T, M <: Manifold, A <: IntersectionAccelerator}
    # TODO: dispatch on manifold
    # this is suitable for planar
    # but spherical / geodesic will need s2 support at some point,
    # or -- even now -- just buffering
    na = GI.npoint(poly_a)
    nb = GI.npoint(poly_b)

    accelerator = if accelerator isa AutoAccelerator
        if na < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS && nb < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS
            NestedLoop()
        elseif na < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS || nb < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS
            SingleNaturalTree()
        else
            DoubleNaturalTree()
        end
    else
        accelerator
    end

    if accelerator isa NestedLoop
        # if we don't have enough vertices in either of the polygons to merit a tree,
        # then we can just do a simple nested loop
        # this becomes extremely useful in e.g. regridding, 
        # where we know the polygon will only ever have a few vertices.
        # This is also applicable to any manifold, since the checking is done within
        # the loop.
        # First, loop over "each edge" in poly_a
        for (i, (a1t, a2t)) in enumerate(eachedge(poly_a, T))
            a1t == a2t && continue
            isnothing(f_on_each_a) ||f_on_each_a(a1t, i)
            for (j, (b1t, b2t)) in enumerate(eachedge(poly_b, T))
                b1t == b2t && continue
                LoopStateMachine.@controlflow f_on_each_maybe_intersect(((a1t, a2t), i), ((b1t, b2t), j)) # this should be aware of manifold by construction.
            end
            isnothing(f_after_each_a) || f_after_each_a(a1t, i)
        end
        # And we're done!
    elseif accelerator isa SingleSTRtree
        # This is the "middle ground" case - run only a strtree 
        # on poly_b without doing so on poly_a.
        # This is less complex than running a dual tree traversal,
        # and reduces the overhead of constructing an edge list and tree on poly_a.
        ext_a, ext_b = GI.extent(poly_a), GI.extent(poly_b)
        edges_b, indices_b = to_edgelist(ext_a, poly_b, T)
        if isempty(edges_b) && !isnothing(f_on_each_a) && !isnothing(f_after_each_a)
            # shortcut - nothing can possibly intersect
            # so we just call f_on_each_a for each edge in poly_a
            for i in 1:GI.npoint(poly_a)-1
                pt = _tuple_point(GI.getpoint(poly_a, i), T)
                f_on_each_a(pt, i)
                f_after_each_a(pt, i)
            end
            return nothing
        end

        # This is the STRtree generated from the edges of poly_b
        tree_b = STRtree(edges_b)

        # this is a pre-allocation that will store the resuits of the query into tree_b
        query_result = Int[] 
        
        # Loop over each vertex in poly_a
        for (i, (a1t, a2t)) in enumerate(eachedge(poly_a, T))
            a1t == a2t && continue
            l1 = GI.Line(SVector{2}(a1t, a2t))
            ext_l = GI.extent(l1)
            # l = GI.Line(SVector{2}(a1t, a2t); extent=ext_l) # this seems to be unused - TODO remove
            isnothing(f_on_each_a) || f_on_each_a(a1t, i)
            # Query the STRtree for any edges in b that may intersect this edge
            # This is sorted because we want to pretend we're doing the same thing
            # as the nested loop above, and iterating through poly_b in order.
            if Extents.intersects(ext_l, ext_b)
                empty!(query_result)
                SortTileRecursiveTree.query!(query_result, tree_b.rootnode, ext_l)
                sort!(query_result) # STRTree.jl's query! does not sort!, even though query does...
                # Loop over the edges in b that might intersect the edges in a
                for j in query_result
                    b1t, b2t = edges_b[j].geom
                    b1t == b2t && continue
                    # Manage control flow if the function returns a LoopStateMachine.Action
                    # like Break(), Continue(), or Return()
                    # This allows the function to break out of the loop early if it wants
                    # without being syntactically inside the loop.
                    LoopStateMachine.@controlflow f_on_each_maybe_intersect(((a1t, a2t), i), ((b1t, b2t), indices_b[j])) # note the indices_b[j] here - we are using the index of the edge in the original edge list, not the index of the edge in the STRtree.
                end
            end
            isnothing(f_after_each_a) || f_after_each_a(a1t, i)
        end
    elseif accelerator isa SingleNaturalTree
        ext_a, ext_b = GI.extent(poly_a), GI.extent(poly_b)
        edges_b = to_edgelist(poly_b, T)

        b_tree = NaturalIndexing.NaturalIndex(edges_b)

        for (i, (a1t, a2t)) in enumerate(eachedge(poly_a, T))
            a1t == a2t && continue
            ext_l = Extents.Extent(X = minmax(a1t[1], a2t[1]), Y = minmax(a1t[2], a2t[2]))
            isnothing(f_on_each_a) || f_on_each_a(a1t, i)
            # Query the STRtree for any edges in b that may intersect this edge
            # This is sorted because we want to pretend we're doing the same thing
            # as the nested loop above, and iterating through poly_b in order.
            if Extents.intersects(ext_l, ext_b)
                # Loop over the edges in b that might intersect the edges in a
                SpatialTreeInterface.depth_first_search(Base.Fix1(Extents.intersects, ext_l), b_tree) do j
                    b1t, b2t = edges_b[j].geom
                    b1t == b2t && return LoopStateMachine.Continue()
                    # LoopStateMachine control is managed outside the loop, by the depth_first_search function.
                    return f_on_each_maybe_intersect(((a1t, a2t), i), ((b1t, b2t), j)) # note the indices_b[j] here - we are using the index of the edge in the original edge list, not the index of the edge in the STRtree.
                end
            end
            isnothing(f_after_each_a) || f_after_each_a(a1t, i)
        end

    elseif accelerator isa DoubleNaturalTree
        edges_a = to_edgelist(poly_a, T)
        edges_b = to_edgelist(poly_b, T)

        tree_a = NaturalIndexing.NaturalIndex(edges_a)
        tree_b = NaturalIndexing.NaturalIndex(edges_b)

        last_a_idx::Int = 0

        SpatialTreeInterface.dual_depth_first_search(Extents.intersects, tree_a, tree_b) do a_edge_idx, b_edge_idx
            a1t, a2t = edges_a[a_edge_idx].geom
            b1t, b2t = edges_b[b_edge_idx].geom

            if last_a_idx < a_edge_idx
                if !isnothing(f_on_each_a)
                    for i in (last_a_idx+1):(a_edge_idx-1)
                        f_on_each_a((edges_a[i].geom[1]), i)
                        !isnothing(f_after_each_a) && f_after_each_a((edges_a[i].geom[1]), i)
                    end
                end
                !isnothing(f_on_each_a) && f_on_each_a(a1t, a_edge_idx)
            end

            f_on_each_maybe_intersect(((a1t, a2t), a_edge_idx), ((b1t, b2t), b_edge_idx))

            if last_a_idx < a_edge_idx
                if !isnothing(f_after_each_a)
                    f_after_each_a(a1t, a_edge_idx)
                end
                last_a_idx = a_edge_idx
            end
        end

        if last_a_idx == 0 # the query did not find any intersections
            if !isnothing(f_on_each_a) && isnothing(f_after_each_a)
                return
            else
                for (i, edge) in enumerate(edges_a)
                    !isnothing(f_on_each_a) && f_on_each_a(edge.geom[1], i)
                    !isnothing(f_after_each_a) && f_after_each_a(edge.geom[1], i)
                end
            end
        elseif last_a_idx < length(edges_a)
            # the query terminated early - this will almost always be the case.
            if !isnothing(f_on_each_a) && isnothing(f_after_each_a)
                return
            else
                for (i, edge) in zip(last_a_idx+1:length(edges_a), view(edges_a, last_a_idx+1:length(edges_a)))
                    !isnothing(f_on_each_a) && f_on_each_a(edge.geom[1], i)
                    !isnothing(f_after_each_a) && f_after_each_a(edge.geom[1], i)
                end
            end
        end
    elseif accelerator isa ThinnedDoubleNaturalTree
        ext_a, ext_b = GI.extent(poly_a), GI.extent(poly_b)
        mutual_extent = Extents.intersection(ext_a, ext_b)

        edges_a, indices_a = to_edgelist(mutual_extent, poly_a, T)
        edges_b, indices_b = to_edgelist(mutual_extent, poly_b, T)

        tree_a = NaturalIndexing.NaturalIndex(edges_a)
        tree_b = NaturalIndexing.NaturalIndex(edges_b)

        last_a_idx::Int = 1

        SpatialTreeInterface.dual_depth_first_search(Extents.intersects, tree_a, tree_b) do a_thinned_idx, b_thinned_idx
            a_edge_idx = indices_a[a_thinned_idx]
            b_edge_idx = indices_b[b_thinned_idx]

            a1t, a2t = edges_a[a_thinned_idx].geom
            b1t, b2t = edges_b[b_thinned_idx].geom

            if last_a_idx < a_edge_idx
                if !isnothing(f_on_each_a)
                    for i in last_a_idx:(a_edge_idx-1)
                        f_on_each_a(a1t, a_edge_idx)
                        !isnothing(f_after_each_a) && f_after_each_a(a1t, a_edge_idx)
                    end
                end
                !isnothing(f_on_each_a) && f_on_each_a(a1t, a_edge_idx)
            end

            f_on_each_maybe_intersect(((a1t, a2t), a_edge_idx), ((b1t, b2t), b_edge_idx))

            if last_a_idx < a_edge_idx
                if !isnothing(f_after_each_a)
                    f_after_each_a(a1t, a_edge_idx)
                end
                last_a_idx = a_edge_idx
            end
        end
    else
        error("Unsupported accelerator type: $accelerator.  FosterHormannClipping only supports NestedLoop() or SingleSTRtree().")
    end

    return nothing

end

#=
    _build_a_list(::Type{T}, poly_a, poly_b) -> (a_list, a_idx_list)

This function take in two polygon rings and creates a vector of PolyNodes to represent
poly_a, including its intersection points with poly_b. The information stored in each
PolyNode is needed for clipping using the Greiner-Hormann clipping algorithm.
    
Note: After calling this function, a_list is not fully formed because the neighboring
indices of the intersection points in b_list still need to be updated. Also we still have
not update the entry and exit flags for a_list.
    
The a_idx_list is a list of the indices of intersection points in a_list. The value at
index i of a_idx_list is the location in a_list where the ith intersection point lies.
=#
function _build_a_list(alg::FosterHormannClipping{M, A}, ::Type{T}, poly_a, poly_b; exact) where {T, M, A}
    n_a_edges = _nedge(poly_a)
    a_list = PolyNode{T}[]  # list of points in poly_a
    sizehint!(a_list, n_a_edges)
    a_idx_list = Vector{Int}()  # finds indices of intersection points in a_list
    local a_count::Int = 0  # number of points added to a_list
    local n_b_intrs::Int = 0
    local prev_counter::Int = 0

    function on_each_a(a_pt, i)
        new_point = PolyNode{T}(;point = a_pt)
        a_count += 1
        push!(a_list, new_point)
        prev_counter = a_count
        return nothing
    end

    function after_each_a(a_pt, i)
        # Order intersection points by placement along edge using fracs value
        if prev_counter < a_count
            Δintrs = a_count - prev_counter
            inter_points = @view a_list[(a_count - Δintrs + 1):a_count]
            sort!(inter_points, by = x -> x.fracs[1])
        end
        return nothing
    end

    function on_each_maybe_intersect(((a_pt1, a_pt2), i), ((b_pt1, b_pt2), j))
        if (b_pt1 == b_pt2)  # don't repeat points
            b_pt1 = b_pt2
            return
        end
        # Determine if edges intersect and how they intersect
        line_orient, intr1, intr2 = _intersection_point(alg.manifold, T, (a_pt1, a_pt2), (b_pt1, b_pt2); exact)
        if line_orient != line_out  # edges intersect
            if line_orient == line_cross  # Intersection point that isn't a vertex
                int_pt, fracs = intr1
                new_intr = PolyNode{T}(;
                    point = int_pt, inter = true, neighbor = j, # j is now equivalent to old j-1
                    crossing = true, fracs = fracs,
                )
                a_count += 1
                n_b_intrs += 1
                push!(a_list, new_intr)
                push!(a_idx_list, a_count)
            else
                (_, (α1, β1)) = intr1
                # Determine if a1 or b1 should be added to a_list
                add_a1 = α1 == 0 && 0 ≤ β1 < 1
                a1_β = add_a1 ? β1 : zero(T)
                add_b1 = β1 == 0 && 0 < α1 < 1
                b1_α = add_b1 ? α1 : zero(T)
                # If lines are collinear and overlapping, a second intersection exists
                if line_orient == line_over
                    (_, (α2, β2)) = intr2
                    if α2 == 0 && 0 ≤ β2 < 1
                        add_a1, a1_β = true, β2
                    end
                    if β2 == 0 && 0 < α2 < 1
                        add_b1, b1_α = true, α2
                    end
                end
                # Add intersection points determined above
                if add_a1
                    n_b_intrs += a1_β == 0 ? 0 : 1
                    a_list[prev_counter] = PolyNode{T}(;
                        point = a_pt1, inter = true, neighbor = j,
                        fracs = (zero(T), a1_β),
                    )
                    push!(a_idx_list, prev_counter)
                end
                if add_b1
                    new_intr = PolyNode{T}(;
                        point = b_pt1, inter = true, neighbor = j,
                        fracs = (b1_α, zero(T)),
                    )
                    a_count += 1
                    push!(a_list, new_intr)
                    push!(a_idx_list, a_count)
                end
            end
        end
        return nothing
    end

    # do the iteration but in an accelerated way
    # this is equivalent to (but faster than)
    #=
    ```julia
    for ((a1, a2), i) in eachedge(poly_a)
        on_each_a(a1, i)
        for ((b1, b2), j) in eachedge(poly_b)
            on_each_maybe_intersect(((a1, a2), i), ((b1, b2), j))
        end
        after_each_a(a1, i)
    end
    ```
    =#
    foreach_pair_of_maybe_intersecting_edges_in_order(alg, on_each_a, after_each_a, on_each_maybe_intersect, poly_a, poly_b, T)

    return a_list, a_idx_list, n_b_intrs
end

#=
    _build_b_list(::Type{T}, a_idx_list, a_list, poly_b) -> b_list

This function takes in the a_list and a_idx_list build in _build_a_list and poly_b and
creates a vector of PolyNodes to represent poly_b. The information stored in each PolyNode
is needed for clipping using the Greiner-Hormann clipping algorithm.
    
Note: after calling this function, b_list is not fully updated. The entry/exit flags still
need to be updated. However, the neighbor value in a_list is now updated.
=#
function _build_b_list(alg::FosterHormannClipping{M, A}, ::Type{T}, a_idx_list, a_list, n_b_intrs, poly_b) where {T, M, A} 
    # Sort intersection points by insertion order in b_list
    sort!(a_idx_list, by = x-> a_list[x].neighbor + a_list[x].fracs[2])
    # Initialize needed values and lists
    n_b_edges = _nedge(poly_b)
    n_intr_pts = length(a_idx_list)
    b_list = PolyNode{T}[]
    sizehint!(b_list, n_b_edges + n_b_intrs)
    intr_curr = 1
    b_count = 0
    # Loop over points in poly_b and add each point and intersection point
    local b_pt1
    for (i, b_p2) in enumerate(GI.getpoint(poly_b))
        b_pt2 = _tuple_point(b_p2, T)
        if i ≤ 1 || (b_pt1 == b_pt2)  # don't repeat points
            b_pt1 = b_pt2
            continue
        end
        b_count += 1
        push!(b_list, PolyNode{T}(; point = b_pt1))
        if intr_curr ≤ n_intr_pts
            curr_idx = a_idx_list[intr_curr]
            curr_node = a_list[curr_idx]
            prev_counter = b_count
            while curr_node.neighbor == i - 1  # Add all intersection points on current edge
                b_idx = 0
                new_intr = PolyNode(curr_node; neighbor = curr_idx)
                if curr_node.fracs[2] == 0  # if curr_node is segment start point
                    # intersection point is vertex of b
                    b_idx = prev_counter
                    b_list[b_idx] = new_intr
                else
                    b_count += 1
                    b_idx = b_count
                    push!(b_list, new_intr)
                end
                a_list[curr_idx] = PolyNode(curr_node; neighbor = b_idx)
                intr_curr += 1
                intr_curr > n_intr_pts && break
                curr_idx = a_idx_list[intr_curr]
                curr_node = a_list[curr_idx]
            end
        end
        b_pt1 = b_pt2
    end
    sort!(a_idx_list)  # return a_idx_list to order of points in a_list
    return b_list
end

#=
    _classify_crossing!(T, poly_b, a_list; exact)

This function marks all intersection points as either bouncing or crossing points. "Delayed"
crossing or bouncing intersections (a chain of edges where the central edges overlap and
thus only the first and last edge of the chain determine if the chain is bounding or
crossing) are marked as follows: the first and the last points are marked as crossing if the
chain is crossing and delayed otherwise and all middle points are marked as bouncing.
Additionally, the start and end points of the chain are marked as endpoints using the
endpoints field. 
=#
function _classify_crossing!(alg::FosterHormannClipping{M, A}, ::Type{T}, a_list, b_list; exact) where {T, M, A}
    napts = length(a_list)
    nbpts = length(b_list)
    # start centered on last point
    a_prev = a_list[end - 1]
    curr_pt = a_list[end]
    i = napts
    # keep track of unmatched bouncing chains
    start_chain_edge, start_chain_idx = unknown, 0
    unmatched_end_chain_edge, unmatched_end_chain_idx = unknown, 0
    same_winding = true
    # loop over list points
    for next_idx in 1:napts
        a_next = a_list[next_idx]
        if curr_pt.inter && !curr_pt.crossing
            j = curr_pt.neighbor
            b_prev = j == 1 ? b_list[end] : b_list[j-1]
            b_next = j == nbpts ? b_list[1] : b_list[j+1]
            # determine if any segments are on top of one another
            a_prev_is_b_prev = a_prev.inter && equals(a_prev, b_prev)
            a_prev_is_b_next = a_prev.inter && equals(a_prev, b_next)
            a_next_is_b_prev = a_next.inter && equals(a_next, b_prev)
            a_next_is_b_next = a_next.inter && equals(a_next, b_next)
            # determine which side of a segments the p points are on
            b_prev_side, b_next_side = _get_sides(#=TODO: alg.manifold, =#b_prev, b_next, a_prev, curr_pt, a_next,
                i, j, a_list, b_list; exact)
            # no sides overlap
            if !a_prev_is_b_prev && !a_prev_is_b_next && !a_next_is_b_prev && !a_next_is_b_next
                if b_prev_side != b_next_side  # lines cross 
                    a_list[i] = PolyNode(curr_pt; crossing = true)
                    b_list[j] = PolyNode(b_list[j]; crossing = true)
                end
            # end of overlapping chain
            elseif !a_next_is_b_prev && !a_next_is_b_next 
                b_side = a_prev_is_b_prev ? b_next_side : b_prev_side
                if start_chain_edge == unknown  # start loop on overlapping chain
                    unmatched_end_chain_edge = b_side
                    unmatched_end_chain_idx = i
                    same_winding = a_prev_is_b_prev
                else  # close overlapping chain
                    # update end of chain with endpoint and crossing / bouncing tags
                    crossing = b_side != start_chain_edge
                    a_list[i] = PolyNode(curr_pt;
                        crossing = crossing,
                        endpoint = end_chain,
                    )
                    b_list[j] = PolyNode(b_list[j];
                        crossing = crossing,
                        endpoint = same_winding ? end_chain : start_chain,
                    )
                    # update start of chain with endpoint and crossing / bouncing tags
                    start_pt = a_list[start_chain_idx]
                    a_list[start_chain_idx] = PolyNode(start_pt;
                        crossing = crossing,
                        endpoint = start_chain,
                    )
                    b_list[start_pt.neighbor] = PolyNode(b_list[start_pt.neighbor];
                        crossing = crossing,
                        endpoint = same_winding ? start_chain : end_chain,
                    )
                end
            # start of overlapping chain
            elseif !a_prev_is_b_prev && !a_prev_is_b_next
                b_side = a_next_is_b_prev ? b_next_side : b_prev_side
                start_chain_edge = b_side
                start_chain_idx = i
                same_winding = a_next_is_b_next
            end
        end
        a_prev = curr_pt
        curr_pt = a_next
        i = next_idx
    end
    # if we started in the middle of overlapping chain, close chain
    if unmatched_end_chain_edge != unknown
        crossing = unmatched_end_chain_edge != start_chain_edge
        # update end of chain with endpoint and crossing / bouncing tags
        end_chain_pt = a_list[unmatched_end_chain_idx]
        a_list[unmatched_end_chain_idx] = PolyNode(end_chain_pt;
            crossing = crossing,
            endpoint = end_chain,
        )
        b_list[end_chain_pt.neighbor] = PolyNode(b_list[end_chain_pt.neighbor];
            crossing = crossing,
            endpoint = same_winding ? end_chain : start_chain,
        )
        # update start of chain with endpoint and crossing / bouncing tags
        start_pt = a_list[start_chain_idx]
        a_list[start_chain_idx] = PolyNode(start_pt;
            crossing = crossing,
            endpoint = start_chain,
        )
        b_list[start_pt.neighbor] = PolyNode(b_list[start_pt.neighbor];
            crossing = crossing,
            endpoint = same_winding ? start_chain : end_chain,
        )
    end
end

# Check if PolyNode is a vertex of original polygon
_is_vertex(pt) = !pt.inter || pt.fracs[1] == 0 || pt.fracs[1] == 1 || pt.fracs[2] == 0 || pt.fracs[2] == 1

#= Determines which side (right or left) of the segment a_prev-curr_pt-a_next the points
b_prev and b_next are on. Given this is only called when curr_pt is an intersection point
that wasn't initially classified as crossing, we know that curr_pt is either from a hinge or
overlapping intersection and thus is an original vertex of either poly_a or poly_b. Due to
floating point error when calculating new intersection points, we only want to use original 
vertices to determine orientation. Thus, for other points, find nearest point that is a
vertex. Given other intersection points will be collinear along existing segments, this
won't change the orientation. =#
function _get_sides(b_prev, b_next, a_prev, curr_pt, a_next, i, j, a_list, b_list; exact)
    b_prev_pt = if _is_vertex(b_prev)
        b_prev.point
    else  # Find original start point of segment formed by b_prev and curr_pt
        prev_idx = findprev(_is_vertex, b_list, j - 1)
        prev_idx = isnothing(prev_idx) ? findlast(_is_vertex, b_list) : prev_idx
        b_list[prev_idx].point
    end
    b_next_pt = if _is_vertex(b_next)
        b_next.point
    else  # Find original end point of segment formed by curr_pt and b_next
        next_idx = findnext(_is_vertex, b_list, j + 1)
        next_idx = isnothing(next_idx) ? findfirst(_is_vertex, b_list) : next_idx
        b_list[next_idx].point
    end
    a_prev_pt = if _is_vertex(a_prev)
        a_prev.point
    else   # Find original start point of segment formed by a_prev and curr_pt
        prev_idx = findprev(_is_vertex, a_list, i - 1)
        prev_idx = isnothing(prev_idx) ? findlast(_is_vertex, a_list) : prev_idx
        a_list[prev_idx].point
    end
    a_next_pt = if _is_vertex(a_next)
        a_next.point
    else  # Find original end point of segment formed by curr_pt and a_next
        next_idx = findnext(_is_vertex, a_list, i + 1)
        next_idx = isnothing(next_idx) ? findfirst(_is_vertex, a_list) : next_idx
        a_list[next_idx].point
    end
    # Determine side orientation of b_prev and b_next
    b_prev_side = _get_side(b_prev_pt, a_prev_pt, curr_pt.point, a_next_pt; exact)
    b_next_side = _get_side(b_next_pt, a_prev_pt, curr_pt.point, a_next_pt; exact)
    return b_prev_side, b_next_side
end

# Determines if Q lies to the left or right of the line formed by P1-P2-P3
function _get_side(Q, P1, P2, P3; exact)
    s1 = Predicates.orient(Q, P1, P2; exact)
    s2 = Predicates.orient(Q, P2, P3; exact)
    s3 = Predicates.orient(P1, P2, P3; exact)

    side = if s3 ≥ 0
        (s1 < 0) || (s2 < 0) ? right : left
    else #  s3 < 0
        (s1 > 0) || (s2 > 0) ? left : right
    end
    return side
end

#= Given a list of PolyNodes, find the first element that isn't an intersection point. Then,
test if this element is in or out of the given polygon. Return the next index, as well as
the enter/exit status of the next intersection point (the opposite of the in/out check). If 
all points are intersection points, find the first element that either is the end of a chain
or a crossing point that isn't in a chain. Then take the midpoint of this point and the next
point in the list and perform the in/out check. If none of these points exist, return
a `next_idx` of `nothing`. =#
function _pt_off_edge_status(alg::FosterHormannClipping{M, A}, ::Type{T}, pt_list, poly, npts; exact) where {T, M, A}
    start_idx, is_non_intr_pt = findfirst(_is_not_intr, pt_list), true
    if isnothing(start_idx)
        start_idx, is_non_intr_pt = findfirst(_next_edge_off, pt_list), false
        isnothing(start_idx) && return (start_idx, false)
    end
    next_idx = start_idx < npts ? (start_idx + 1) : 1
    start_pt = if is_non_intr_pt
        pt_list[start_idx].point
    else
        (pt_list[start_idx].point .+ pt_list[next_idx].point) ./ 2
    end
    start_status = !_point_filled_curve_orientation(alg.manifold, start_pt, poly; in = true, on = false, out = false, exact)
    return next_idx, start_status
end
# Check if a PolyNode is an intersection point
_is_not_intr(pt) = !pt.inter
#= Check if a PolyNode is the last point of a chain or a non-overlapping crossing point.
The next midpoint of one of these points and the next point within a polygon must not be on
the polygon edge. =#
_next_edge_off(pt) = (pt.endpoint == end_chain) || (pt.crossing && pt.endpoint == not_endpoint)

#=
    _flag_ent_exit!(::Type{T}, ::GI.LinearRingTrait, poly, pt_list, delay_cross_f, delay_bounce_f; exact)

This function flags all the intersection points as either an 'entry' or 'exit' point in
relation to the given polygon. For non-delayed crossings we simply alternate the enter/exit
status. This also holds true for the first and last points of a delayed bouncing, where they
both have an opposite entry/exit flag. Conversely, the first and last point of a delayed
crossing have the same entry/exit status. Furthermore, the crossing/bouncing flag of delayed
crossings and bouncings may be updated. This depends on function specific rules that
determine which of the start or end points (if any) should be marked as crossing for used
during polygon tracing. A consistent rule is that the start and end points of a delayed
crossing will have different crossing/bouncing flags, while a the endpoints of a delayed
bounce will be the same.

Used for clipping polygons by other polygons.
=#
function _flag_ent_exit!(alg::FosterHormannClipping{M, A}, ::Type{T}, ::GI.LinearRingTrait, poly, pt_list, delay_cross_f, delay_bounce_f; exact) where {T, M, A}
    npts = length(pt_list)
    # Find starting index if there is one
    next_idx, status = _pt_off_edge_status(alg, T, pt_list, poly, npts; exact)
    isnothing(next_idx) && return
    start_idx = next_idx - 1 
    # Loop over points and mark entry and exit status
    start_chain_idx = 0
    for ii in Iterators.flatten((next_idx:npts, 1:start_idx))
        curr_pt = pt_list[ii]
        if curr_pt.endpoint == start_chain
            start_chain_idx = ii
        elseif curr_pt.crossing || curr_pt.endpoint == end_chain
            start_crossing, end_crossing = curr_pt.crossing, curr_pt.crossing
            if curr_pt.endpoint == end_chain  # ending overlapping chain
                start_pt = pt_list[start_chain_idx]
                if curr_pt.crossing  # delayed crossing
                    #= start and end crossing status are different and depend on current
                    entry/exit status =#
                    start_crossing, end_crossing = delay_cross_f(status)
                else  # delayed bouncing
                    next_idx = ii < npts ? (ii + 1) : 1
                    next_val = (curr_pt.point .+ pt_list[next_idx].point) ./ 2
                    pt_in_poly = _point_filled_curve_orientation(alg.manifold, next_val, poly; in = true, on = false, out = false, exact)
                    #= start and end crossing status are the same and depend on if adjacent
                    edges of pt_list are within poly =#
                    start_crossing = delay_bounce_f(pt_in_poly)
                    end_crossing = start_crossing
                end
                # update start of chain point
                pt_list[start_chain_idx] = PolyNode(start_pt; ent_exit = status, crossing = start_crossing)
                if !curr_pt.crossing
                    status = !status
                end
            end
            pt_list[ii] = PolyNode(curr_pt; ent_exit = status, crossing = end_crossing)
            status = !status
        end
    end
    return
end

#=
    _flag_ent_exit!(::GI.LineTrait, line, pt_list; exact)

This function flags all the intersection points as either an 'entry' or 'exit' point in
relation to the given line. Returns true if there are crossing points to classify, else
returns false. Used for cutting polygons by lines.

Assumes that the first point is outside of the polygon and not on an edge.
=#
function _flag_ent_exit!(alg::FosterHormannClipping{M, A}, ::GI.LineTrait, poly, pt_list; exact) where {M, A}
    status = !_point_filled_curve_orientation(#=TODO: alg.manifold=#pt_list[1].point, poly; in = true, on = false, out = false, exact)
    # Loop over points and mark entry and exit status
    for (ii, curr_pt) in enumerate(pt_list)
        if curr_pt.crossing
            pt_list[ii] = PolyNode(curr_pt; ent_exit = status)
            status = !status
        end
    end
    return
end

#= Filters a_idx_list to just include crossing points and sets the index of all crossing
points (which element they correspond to within a_idx_list). =#
function _index_crossing_intrs!(alg::FosterHormannClipping{M, A}, a_list, b_list, a_idx_list) where {M, A}
    filter!(x -> a_list[x].crossing, a_idx_list)
    for (i, a_idx) in enumerate(a_idx_list)
        curr_node = a_list[a_idx]
        neighbor_node = b_list[curr_node.neighbor]
        a_list[a_idx] = PolyNode(curr_node; idx = i)
        b_list[curr_node.neighbor] = PolyNode(neighbor_node; idx = i)
    end
    return
end

#=
    _trace_polynodes(::Type{T}, a_list, b_list, a_idx_list, f_step)::Vector{GI.Polygon}

This function takes the outputs of _build_ab_list and traces the lists to determine which
polygons are formed as described in Greiner and Hormann. The function f_step determines in
which direction the lists are traced.  This function is different for intersection,
difference, and union. f_step must take in two arguments: the most recent intersection
node's entry/exit status and a boolean that is true if we are currently tracing a_list and
false if we are tracing b_list. The functions used for each clipping operation are follows:
    - Intersection: (x, y) -> x ? 1 : (-1)
    - Difference: (x, y) -> (x ⊻ y) ? 1 : (-1)
    - Union: (x, y) -> x ? (-1) : 1

A list of GeoInterface polygons is returned from this function. 

Note: `poly_a` and `poly_b` are temporary inputs used for debugging and can be removed
eventually.
=#
function _trace_polynodes(alg::FosterHormannClipping{M, A}, ::Type{T}, a_list, b_list, a_idx_list, f_step, poly_a, poly_b) where {T, M, A}
    n_a_pts, n_b_pts = length(a_list), length(b_list)
    total_pts = n_a_pts + n_b_pts
    n_cross_pts = length(a_idx_list)
    return_polys = Vector{_get_poly_type(T)}(undef, 0)
    # Keep track of number of processed intersection points
    visited_pts = 0
    processed_pts = 0
    first_idx = 1
    while processed_pts < n_cross_pts
        curr_list, curr_npoints = a_list, n_a_pts
        on_a_list = true
        # Find first unprocessed intersecting point in subject polygon
        visited_pts += 1
        processed_pts += 1
        first_idx = findnext(x -> x != 0, a_idx_list, first_idx)
        idx = a_idx_list[first_idx]
        a_idx_list[first_idx] = 0
        start_pt = a_list[idx]

        # Set first point in polygon
        curr = curr_list[idx]
        pt_list = [curr.point]

        curr_not_start = true
        while curr_not_start
            step = f_step(curr.ent_exit, on_a_list)
            # changed curr_not_intr to curr_not_same_ent_flag
            same_status, prev_status = true, curr.ent_exit
            while same_status
                if visited_pts >= total_pts
                    throw(TracingError("Clipping tracing hit every point - clipping error.", poly_a, poly_b, a_list, b_list, a_idx_list))
                end
                # Traverse polygon either forwards or backwards
                idx += step
                idx = (idx > curr_npoints) ? mod(idx, curr_npoints) : idx
                idx = (idx == 0) ? curr_npoints : idx

                # Get current node and add to pt_list
                curr = curr_list[idx]
                push!(pt_list, curr.point)
                if (curr.crossing || curr.endpoint != not_endpoint)
                    # Keep track of processed intersection points
                    same_status = curr.ent_exit == prev_status
                    curr_not_start = curr != start_pt && curr != b_list[start_pt.neighbor]
                    !curr_not_start && break
                    if (on_a_list && curr.crossing) || (!on_a_list && a_list[curr.neighbor].crossing)
                        processed_pts += 1
                        a_idx_list[curr.idx] = 0
                    end
                end
                visited_pts += 1
            end
            # Switch to next list and next point
            curr_list, curr_npoints = on_a_list ? (b_list, n_b_pts) : (a_list, n_a_pts)
            on_a_list = !on_a_list
            idx = curr.neighbor
            curr = curr_list[idx]
        end
        push!(return_polys, GI.Polygon([pt_list]))
    end
    return return_polys
end

# Get type of polygons that will be made
# TODO: Increase type options
_get_poly_type(::Type{T}) where T =
    GI.Polygon{false, false, Vector{GI.LinearRing{false, false, Vector{Tuple{T, T}}, Nothing, Nothing}}, Nothing, Nothing}

#=
    _find_non_cross_orientation(a_list, b_list, a_poly, b_poly; exact)

For polygons with no crossing intersection points, either one polygon is inside of another,
or they are separate polygons with no intersection (other than an edge or point).

Return two booleans that represent if a is inside b (potentially with shared edges / points)
and visa versa if b is inside of a.
=#
function _find_non_cross_orientation(m::M, a_list, b_list, a_poly, b_poly; exact) where {M <: Manifold}
    non_intr_a_idx = findfirst(x -> !x.inter, a_list)
    non_intr_b_idx = findfirst(x -> !x.inter, b_list)
    #= Determine if non-intersection point is in or outside of polygon - if there isn't A
    non-intersection point, then all points are on the polygon edge =#
    a_pt_orient = isnothing(non_intr_a_idx) ? point_on :
        _point_filled_curve_orientation(a_list[non_intr_a_idx].point, b_poly; exact)
    b_pt_orient = isnothing(non_intr_b_idx) ? point_on :
        _point_filled_curve_orientation(b_list[non_intr_b_idx].point, a_poly; exact)
    a_in_b = a_pt_orient != point_out && b_pt_orient != point_in
    b_in_a = b_pt_orient != point_out && a_pt_orient != point_in
    return a_in_b, b_in_a
end

_find_non_cross_orientation(alg::FosterHormannClipping{M}, a_list, b_list, a_poly, b_poly; exact) where {M <: Manifold} =
    _find_non_cross_orientation(alg.manifold, a_list, b_list, a_poly, b_poly; exact)

#=
    _add_holes_to_polys!(::Type{T}, return_polys, hole_iterator, remove_poly_idx; exact)

The holes specified by the hole iterator are added to the polygons in the return_polys list.
If this creates more polygons, they are added to the end of the list. If this removes
polygons, they are removed from the list
=#
function _add_holes_to_polys!(alg::FosterHormannClipping{M, A}, ::Type{T}, return_polys, hole_iterator, remove_poly_idx; exact) where {T, M, A}
    n_polys = length(return_polys)
    remove_hole_idx = Int[]
    # Remove set of holes from all polygons
    for i in 1:n_polys
        n_new_per_poly = 0
        for curr_hole in Iterators.map(tuples, hole_iterator) # loop through all holes
            curr_hole = _linearring(curr_hole)
            # loop through all pieces of original polygon (new pieces added to end of list)
            for j in Iterators.flatten((i:i, (n_polys + 1):(n_polys + n_new_per_poly)))
                curr_poly = return_polys[j]
                remove_poly_idx[j] && continue
                curr_poly_ext = GI.nhole(curr_poly) > 0 ? GI.Polygon(StaticArrays.SVector(GI.getexterior(curr_poly))) : curr_poly
                in_ext, on_ext, out_ext = _line_polygon_interactions(#=TODO: alg.manifold=#curr_hole, curr_poly_ext; exact, closed_line = true)
                if in_ext  # hole is at least partially within the polygon's exterior
                    new_hole, new_hole_poly, n_new_pieces = _combine_holes!(alg, T, curr_hole, curr_poly, return_polys, remove_hole_idx)
                    if n_new_pieces > 0
                        append!(remove_poly_idx, falses(n_new_pieces))
                        n_new_per_poly += n_new_pieces
                    end
                    if !on_ext && !out_ext  # hole is completely within exterior
                        push!(curr_poly.geom, new_hole)
                    else  # hole is partially within and outside of polygon's exterior
                        new_polys = difference(alg, curr_poly_ext, new_hole_poly, T; target=GI.PolygonTrait())
                        n_new_polys = length(new_polys) - 1
                        # replace original
                        curr_poly.geom[1] = GI.getexterior(new_polys[1])
                        append!(curr_poly.geom, GI.gethole(new_polys[1]))
                        if n_new_polys > 0  # add any extra pieces
                            append!(return_polys, @view new_polys[2:end])
                            append!(remove_poly_idx, falses(n_new_polys))
                            n_new_per_poly += n_new_polys
                        end
                    end
                # polygon is completely within hole
                elseif coveredby(#=TODO: alg.manifold=#curr_poly_ext, GI.Polygon(StaticArrays.SVector(curr_hole)))
                    remove_poly_idx[j] = true
                end
            end
        end
        n_polys += n_new_per_poly
    end
    # Remove all polygon that were marked for removal
    deleteat!(return_polys, remove_poly_idx)
    return
end

#=
    _combine_holes!(::Type{T}, new_hole, curr_poly, return_polys)

The new hole is combined with any existing holes in curr_poly. The holes can be combined
into a larger hole if they are intersecting. If this happens, then the new, combined hole is
returned with the original holes making up the new hole removed from curr_poly. Additionally,
if the combined holes form a ring, the interior is added to the return_polys as a new
polygon piece. Additionally, holes leftover after combination will be checked for it they
are in the "main" polygon or in one of these new pieces and moved accordingly. 

If the holes don't touch or curr_poly has no holes, then new_hole is returned without any
changes.
=#
function _combine_holes!(alg::FosterHormannClipping{M, A}, ::Type{T}, new_hole, curr_poly, return_polys, remove_hole_idx) where {T, M, A}
    n_new_polys = 0
    empty!(remove_hole_idx)
    new_hole_poly = GI.Polygon(StaticArrays.SVector(new_hole))
    # Combine any existing holes in curr_poly with new hole
    for (k, old_hole) in enumerate(GI.gethole(curr_poly))
        old_hole_poly = GI.Polygon(StaticArrays.SVector(old_hole))
        if intersects(#=TODO: alg.manifold=#new_hole_poly, old_hole_poly)
            # If the holes intersect, combine them into a bigger hole
            hole_union = union(alg, new_hole_poly, old_hole_poly, T; target = GI.PolygonTrait())[1]
            push!(remove_hole_idx, k + 1)
            new_hole = GI.getexterior(hole_union)
            new_hole_poly = GI.Polygon(StaticArrays.SVector(new_hole))
            n_pieces = GI.nhole(hole_union)
            if n_pieces > 0  # if the hole has a hole, then this is a new polygon piece! 
                append!(return_polys, [GI.Polygon([h]) for h in GI.gethole(hole_union)])
                n_new_polys += n_pieces
            end
        end
    end
    # Remove redundant holes
    deleteat!(curr_poly.geom, remove_hole_idx)
    empty!(remove_hole_idx)
    # If new polygon pieces created, make sure remaining holes are in the correct piece
    @views for piece in return_polys[end - n_new_polys + 1:end]
        for (k, old_hole) in enumerate(GI.gethole(curr_poly))
            if !(k in remove_hole_idx) && within(old_hole, piece)
                push!(remove_hole_idx, k + 1)
                push!(piece.geom, old_hole)
            end
        end
    end
    deleteat!(curr_poly.geom, remove_hole_idx)
    return new_hole, new_hole_poly, n_new_polys
end

#= Remove collinear edge points, other than the first and last edge vertex, to simplify
polygon - including both the exterior ring and any holes=#
function _remove_collinear_points!(alg::FosterHormannClipping{M, A}, polys, remove_idx, poly_a, poly_b) where {M, A}
    for (i, poly) in Iterators.reverse(enumerate(polys))
        for (j, ring) in Iterators.reverse(enumerate(GI.getring(poly)))
            n = length(ring.geom)
            # resize and reset removing index buffer
            resize!(remove_idx, n)
            fill!(remove_idx, false)
            local p1, p2
            for (i, p) in enumerate(ring.geom)
                if i == 1
                    p1 = p
                    continue
                elseif i == 2
                    p2 = p
                    continue
                else
                    p3 = p
                    # check if p2 is approximately on the edge formed by p1 and p3 - remove if so
                    # TODO: make this manifold aware
                    if Predicates.orient(p1, p2, p3; exact = False()) == 0
                        remove_idx[i - 1] = true
                    end
                end
                p1, p2 = p2, p3
            end
            # Check if the first point (which is repeated as the last point) is needed 
            if Predicates.orient(ring.geom[end - 1], ring.geom[1], ring.geom[2]; exact = False()) == 0
                remove_idx[1], remove_idx[end] = true, true
            end
            # Remove unneeded collinear points
            deleteat!(ring.geom, remove_idx)
            # Check if enough points are left to form a polygon
            if length(ring.geom) ≤ (remove_idx[1] ? 2 : 3)
                if j == 1
                    deleteat!(polys, i)
                    break
                else
                    deleteat!(poly.geom, j)
                    continue
                end
            end
            if remove_idx[1]  # make sure the last point is repeated
                push!(ring.geom, ring.geom[1])
            end
        end
    end
    return
end
