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
- `manifold::M`: The manifold on which the algorithm operates. `Geodesic` is not supported —
  the algorithm has no geodesic implementation of its intersection primitives — and the
  constructor throws `ArgumentError` for it; use [`Spherical`](@ref) instead.
- `accelerator::A`: The accelerator to use for the algorithm.  Can be `nothing` for automatic choice, or a custom accelerator.
"""
struct FosterHormannClipping{M <: Manifold, A <: IntersectionAccelerator} <: GeometryOpsCore.Algorithm{M}
    manifold::M
    accelerator::A
    # TODO: add exact flag
    # TODO: should exact flag be in the type domain?
    #= Foster-Hormann has no geodesic implementation of `_get_side` and the other clipping
    primitives. Without this check, `FosterHormannClipping(Geodesic())` constructs
    successfully and then throws a bare `MethodError` deep inside the first clip it
    attempts — far from the actual cause, and only on inputs that reach that code path.
    This is the single point every construction path funnels through (see the outer
    constructor immediately below), including direct parametric construction like
    `FosterHormannClipping{Geodesic{Float64}, NestedLoop}(...)`, so the check cannot be
    bypassed. =#
    function FosterHormannClipping{M, A}(manifold::M, accelerator::A) where {M <: Manifold, A <: IntersectionAccelerator}
        manifold isa Geodesic && throw(ArgumentError(
            "FosterHormannClipping does not support the Geodesic manifold ($manifold): Foster-Hormann clipping has no geodesic implementation of its intersection primitives. Use Spherical() instead."
        ))
        return new{M, A}(manifold, accelerator)
    end
end
#= Defining the inner constructor above suppresses Julia's automatically-generated
type-inferring outer constructor (`FosterHormannClipping(manifold::M, accelerator::A)
where {M, A}`), so it has to be written back explicitly — every other constructor below
relies on it to actually build the struct. =#
FosterHormannClipping(manifold::M, accelerator::A) where {M <: Manifold, A <: IntersectionAccelerator} = FosterHormannClipping{M, A}(manifold, accelerator)
FosterHormannClipping(; manifold::Manifold = Planar(), accelerator = nothing) = FosterHormannClipping(manifold, isnothing(accelerator) ? NestedLoop() : accelerator)
FosterHormannClipping(manifold::Manifold, accelerator::Union{Nothing, IntersectionAccelerator} = nothing) = FosterHormannClipping(manifold, isnothing(accelerator) ? NestedLoop() : accelerator)
FosterHormannClipping(accelerator::Union{Nothing, IntersectionAccelerator}) = FosterHormannClipping(Planar(), isnothing(accelerator) ? NestedLoop() : accelerator)
#= Spherical and geodesic manifolds cannot use any of the tree accelerators: an STRtree
indexes planar rectangles, and neither the antimeridian nor the poles survive that. The
automatic choice therefore resolves to `NestedLoop` on the sphere.

This method is deliberately narrow in *both* arguments. Writing it as
`(::Union{Spherical, Geodesic}, ::Union{Nothing, IntersectionAccelerator})` — narrower in
the manifold but wider in the accelerator than the struct's own outer constructor — makes
every two-argument spherical call ambiguous, and with it every one-argument one, since
that forwards through it. =#
FosterHormannClipping(manifold::Union{Spherical, Geodesic}, ::AutoAccelerator) = FosterHormannClipping(manifold, NestedLoop())

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
@kwdef struct PolyNode{T <: AbstractFloat, P}
    point::P                   # the vertex, in whatever representation the manifold computes in
    inter::Bool = false        # If ipt, true, else 0
    neighbor::Int = 0          # If ipt, index of equivalent point in a_list or b_list, else 0
    idx::Int = 0               # If crossing point, index within sorted a_idx_list
    ent_exit::Bool = false     # If ipt, true if enter and false if exit, else false
    crossing::Bool = false     # If ipt, true if intersection crosses from out/in polygon, else false
    endpoint::EndPointType = not_endpoint # If ipt, denotes if point is the start or end of an overlapping chain
    fracs::Tuple{T,T} = (0., 0.) # If ipt, fractions along edges to ipt (a_frac, b_frac), else (0, 0)
    #= 1-based position of this vertex in the ring it was ingested from, or 0 for a computed
    intersection. `Spherical` converts lon/lat input to xyz once at ingress, so the original
    coordinates are no longer recoverable from `point` by anything but a lossy round trip;
    this is what lets egress hand a passthrough vertex back exactly as it arrived. =#
    srcidx::Int32 = Int32(0)
end

#= `PolyNode{T}(; point = p)` picks the point type up from `p`, so the many call sites that
know the float type but not the representation keep their spelling. =#
(::Type{PolyNode{T}})(; point, kwargs...) where {T <: AbstractFloat} =
    PolyNode{T, typeof(point)}(; point, kwargs...)

#= Create a new node with all of the same field values as the given PolyNode unless
alternative values are provided, in which case those should be used. =#
PolyNode(node::PolyNode{T, P};
    point = node.point, inter = node.inter, neighbor = node.neighbor, idx = node.idx,
    ent_exit = node.ent_exit, crossing = node.crossing, endpoint = node.endpoint,
    fracs = node.fracs, srcidx = node.srcidx,
) where {T, P} = PolyNode{T, typeof(point)}(;
    point = point, inter = inter, neighbor = neighbor, idx = idx, ent_exit = ent_exit,
    crossing = crossing, endpoint = endpoint, fracs = fracs, srcidx = srcidx)

# Checks equality of two PolyNodes by backing point value, fractional value, and intersection status
equals(pn1::PolyNode, pn2::PolyNode) = pn1.point == pn2.point && pn1.inter == pn2.inter && pn1.fracs == pn2.fracs
Base.:(==)(pn1::PolyNode, pn2::PolyNode) = equals(pn1, pn2)

"""
    FosterHormannCache{T}()
    FosterHormannCache([T = Float64])
    FosterHormannCache(alg::FosterHormannClipping, [T = Float64])

Preallocated buffers for [`FosterHormannClipping`](@ref).

Pass this as the `cache` keyword argument to [`intersection_area`](@ref) to reuse the
algorithm's working set instead of allocating it per call — worth it when measuring many
polygon pairs in a hot loop, which is what conservative regridding between two discrete
global grids does for every cell pair.

The buffers are the vertex lists of the two rings and the index of intersections within the
first. They are scratch: nothing the call returns points into them, so a result stays valid
after the cache is reused.

`T` must match the float type of the clip. `FosterHormannCache(alg, T)` spells that out.

!!! warning "Thread safety"
    A cache must not be shared across concurrent tasks. Create one per task. The default
    (`cache = nothing`) allocates per call and is always safe.

# Example

```julia
import GeometryOps as GO

alg = GO.FosterHormannClipping(GO.Spherical())
cache = GO.FosterHormannCache(alg)
for (a, b) in cell_pairs
    frac = GO.intersection_area(alg, a, b; cache)
end
```
"""
struct FosterHormannCache{T, P}
    a_list::Vector{PolyNode{T, P}}
    b_list::Vector{PolyNode{T, P}}
    a_idx_list::Vector{Int}
end
FosterHormannCache{T, P}() where {T, P} =
    FosterHormannCache{T, P}(PolyNode{T, P}[], PolyNode{T, P}[], Int[])

#= The representation the manifold computes in, mirroring `SutherlandHodgmanCache`: planar
clipping works in the chart, spherical clipping works on the unit sphere. =#
_fh_point_type(::Planar, ::Type{T}) where {T} = Tuple{T, T}
_fh_point_type(::Spherical, ::Type{T}) where {T} = UnitSpherical.UnitSphericalPoint{T}
FosterHormannCache(m::Manifold, ::Type{T} = Float64) where {T <: AbstractFloat} =
    FosterHormannCache{T, _fh_point_type(m, T)}()
FosterHormannCache(alg::FosterHormannClipping, ::Type{T} = Float64) where {T <: AbstractFloat} =
    FosterHormannCache(alg.manifold, T)
FosterHormannCache(::Type{T} = Float64) where {T <: AbstractFloat} =
    FosterHormannCache{T, Tuple{T, T}}()

#= Put an input vertex into the representation the manifold computes in. Both legs are
already identities where they can be: `_tuple_point` passes a `UnitSphericalPoint` through
untouched, and `_spherical_kernel_point` is the identity on one, so 3D input reaches the
node lists without arithmetic and comes back out bit-identical. =#
_fh_ingest(::Planar, p, ::Type{T}) where {T} = _tuple_point(p, T)
_fh_ingest(::Spherical, p, ::Type{T}) where {T} = _spherical_kernel_point(_tuple_point(p, T))

function _fh_check_cache(cache::FosterHormannCache{C, Q}, ::Type{T}, ::Type{P}) where {C, Q, T, P}
    (C === T && Q === P) || throw(ArgumentError(
        "FosterHormannCache type mismatch: this clip requires " *
        "FosterHormannCache{$T, $P}, got FosterHormannCache{$C, $Q}. Construct the cache " *
        "with `FosterHormannCache(alg, T)` to match the algorithm."))
    return cache
end

#-- Hand back a cleared buffer from the cache, or a fresh one when there is no cache. Both
#-- branches return the same type, so the caller stays inferable either way.
_fh_buffer(::Nothing, ::Type{V}) where {V} = V()
_fh_buffer(v::Vector, ::Type{V}) where {V} = (empty!(v); v)

# Finally, we define a nice error type for when the clipping tracing algorithm hits every point in a polygon.
# This stores the polygons, the a_list, and the b_list, and the a_idx_list.
# allowing the user to understand what happened and why.
"""
    TracingError{T1, T2} <: Exception

An error that is thrown when the clipping tracing algorithm fails somehow.
This is a bug in the algorithm, and should be reported.

The polygons are contained in the exception object, accessible by try-catch or as `err` in the REPL.
"""
struct TracingError{T1, T2, L1 <: AbstractVector{<:PolyNode}, L2 <: AbstractVector{<:PolyNode}} <: Exception
    message::String
    poly_a::T1
    poly_b::T2
    a_list::L1
    b_list::L2
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
function _build_ab_list(alg::FosterHormannClipping, ::Type{T}, poly_a, poly_b, delay_cross_f::F1, delay_bounce_f::F2; exact, cache = nothing) where {T, F1, F2}
    # Make a list for nodes of each polygon
    a_list, a_idx_list, n_b_intrs = _build_a_list(alg, T, poly_a, poly_b; exact, cache)
    b_list = _build_b_list(alg, T, a_idx_list, a_list, n_b_intrs, poly_b; cache)

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
    manifold::M, accelerator::AutoAccelerator, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64
) where {FA, FAAfter, FI, T, M <: Manifold}
    # this is suitable for planar
    # but spherical / geodesic will need s2 support at some point,
    # or -- even now -- just buffering
    na = GI.npoint(poly_a)
    nb = GI.npoint(poly_b)
    # Switching behaviour is turned off in the patch release
    # This should be turned on in a GO v0.2.x
    if na < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS && nb < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS
        return foreach_pair_of_maybe_intersecting_edges_in_order(manifold, NestedLoop(), f_on_each_a, f_after_each_a, f_on_each_maybe_intersect, poly_a, poly_b, T)
    elseif na < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS || nb < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS
        return foreach_pair_of_maybe_intersecting_edges_in_order(manifold, SingleNaturalTree(), f_on_each_a, f_after_each_a, f_on_each_maybe_intersect, poly_a, poly_b, T)
    else
        return foreach_pair_of_maybe_intersecting_edges_in_order(manifold, DoubleNaturalTree(), f_on_each_a, f_after_each_a, f_on_each_maybe_intersect, poly_a, poly_b, T)
    end
end

function foreach_pair_of_maybe_intersecting_edges_in_order(
    manifold::M, accelerator::NestedLoop, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64
) where {FA, FAAfter, FI, T, M <: Manifold}
    # this is suitable for planar
    # but spherical / geodesic will need s2 support at some point,
    # or -- even now -- just buffering
    na = GI.npoint(poly_a)
    nb = GI.npoint(poly_b)
    # if we don't have enough vertices in either of the polygons to merit a tree,
    # then we can just do a simple nested loop
    # this becomes extremely useful in e.g. regridding, 
    # where we know the polygon will only ever have a few vertices.
    # This is also applicable to any manifold, since the checking is done within
    # the loop.
    # First, loop over "each edge" in poly_a
    for (i, (a1t, a2t)) in enumerate(eachedge(manifold, poly_a, T))
        a1t == a2t && continue
        isnothing(f_on_each_a) || f_on_each_a(a1t, i)
        for (j, (b1t, b2t)) in enumerate(eachedge(manifold, poly_b, T))
            b1t == b2t && continue
            LoopStateMachine.@controlflow f_on_each_maybe_intersect(((a1t, a2t), i), ((b1t, b2t), j)) # this should be aware of manifold by construction.
        end
        isnothing(f_after_each_a) || f_after_each_a(a1t, i)
    end
    # And we're done!  This is the super simple implementation.
    return nothing
end

function foreach_pair_of_maybe_intersecting_edges_in_order(
    manifold::M, accelerator::SingleSTRtree, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64
) where {FA, FAAfter, FI, T, M <: Manifold}
    na = GI.npoint(poly_a)
    nb = GI.npoint(poly_b)
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
    return nothing
end

function foreach_pair_of_maybe_intersecting_edges_in_order(
    manifold::M, accelerator::SingleNaturalTree, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64
) where {FA, FAAfter, FI, T, M <: Manifold}
    na = GI.npoint(poly_a)
    nb = GI.npoint(poly_b)
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
    return nothing
end

function foreach_pair_of_maybe_intersecting_edges_in_order(
    manifold::M, accelerator::DoubleNaturalTree, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64
) where {FA, FAAfter, FI, T, M <: Manifold}
    na = GI.npoint(poly_a)
    nb = GI.npoint(poly_b)
    edges_a = to_edgelist(poly_a, T)
    edges_b = to_edgelist(poly_b, T)

    tree_a = NaturalIndexing.NaturalIndex(edges_a)
    tree_b = NaturalIndexing.NaturalIndex(edges_b)

    last_a_idx = 0

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
    return nothing
end
    
function foreach_pair_of_maybe_intersecting_edges_in_order(
    manifold::M, accelerator::ThinnedDoubleNaturalTree, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64
) where {FA, FAAfter, FI, T, M <: Manifold}
    na = GI.npoint(poly_a)
    nb = GI.npoint(poly_b)
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
function _build_a_list(alg::FosterHormannClipping{M, A}, ::Type{T}, poly_a, poly_b; exact, cache = nothing) where {T, M, A}
    n_a_edges = _nedge(poly_a)
    # list of points in poly_a
    P = _fh_point_type(alg.manifold, T)
    a_list = _fh_buffer(cache === nothing ? nothing : cache.a_list, Vector{PolyNode{T, P}})
    #-- A cached buffer already carries the capacity its last call grew it to, and
    #-- `sizehint!` is free to *shrink* to the hint, which would hand back the storage this
    #-- cache exists to keep and realloc it again on the next push.
    cache === nothing && sizehint!(a_list, n_a_edges)
    # finds indices of intersection points in a_list
    a_idx_list = _fh_buffer(cache === nothing ? nothing : cache.a_idx_list, Vector{Int})
    local a_count::Int = 0  # number of points added to a_list
    local n_b_intrs::Int = 0
    local prev_counter::Int = 0

    function on_each_a(a_pt, i)
        #-- edge `i` starts at vertex `i`, so this is the slot egress reads back
        new_point = PolyNode{T}(;point = a_pt, srcidx = Int32(i))
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
                    #= This promotes the vertex already sitting at `prev_counter` -- same
                    point -- so copy it rather than rebuild it, keeping its source slot. =#
                    a_list[prev_counter] = PolyNode(a_list[prev_counter];
                        inter = true, neighbor = j, fracs = (zero(T), a1_β),
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
function _build_b_list(alg::FosterHormannClipping{M, A}, ::Type{T}, a_idx_list, a_list, n_b_intrs, poly_b; cache = nothing) where {T, M, A}
    # Sort intersection points by insertion order in b_list
    sort!(a_idx_list, by = x-> a_list[x].neighbor + a_list[x].fracs[2])
    # Initialize needed values and lists
    n_b_edges = _nedge(poly_b)
    n_intr_pts = length(a_idx_list)
    P = _fh_point_type(alg.manifold, T)
    b_list = _fh_buffer(cache === nothing ? nothing : cache.b_list, Vector{PolyNode{T, P}})
    cache === nothing && sizehint!(b_list, n_b_edges + n_b_intrs)
    intr_curr = 1
    b_count = 0
    # Loop over points in poly_b and add each point and intersection point
    local b_pt1
    for (i, b_p2) in enumerate(GI.getpoint(poly_b))
        b_pt2 = _fh_ingest(alg.manifold, b_p2, T)
        if i ≤ 1 || (b_pt1 == b_pt2)  # don't repeat points
            b_pt1 = b_pt2
            continue
        end
        b_count += 1
        #-- `b_pt1` is the vertex from the previous step, i.e. slot `i - 1`
        push!(b_list, PolyNode{T}(; point = b_pt1, srcidx = Int32(i - 1)))
        if intr_curr ≤ n_intr_pts
            curr_idx = a_idx_list[intr_curr]
            curr_node = a_list[curr_idx]
            prev_counter = b_count
            while curr_node.neighbor == i - 1  # Add all intersection points on current edge
                b_idx = 0
                #-- `curr_node`'s slot indexes poly_a; in b_list it would resolve
                #-- against the wrong ring, so this node converts at egress instead
                new_intr = PolyNode(curr_node; neighbor = curr_idx, srcidx = Int32(0))
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
            b_prev_side, b_next_side = _get_sides(alg.manifold, b_prev, b_next, a_prev, curr_pt, a_next,
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
function _get_sides(m::Manifold, b_prev, b_next, a_prev, curr_pt, a_next, i, j, a_list, b_list; exact)
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
    b_prev_side = _get_side(m, b_prev_pt, a_prev_pt, curr_pt.point, a_next_pt; exact)
    b_next_side = _get_side(m, b_next_pt, a_prev_pt, curr_pt.point, a_next_pt; exact)
    return b_prev_side, b_next_side
end

# Determines if Q lies to the left or right of the line formed by P1-P2-P3
function _get_side(::Planar, Q, P1, P2, P3; exact)
    s1 = Predicates.orient(Q, P1, P2; exact)
    s2 = Predicates.orient(Q, P2, P3; exact)
    s3 = Predicates.orient(P1, P2, P3; exact)

    return _side_from_orientations(s1, s2, s3)
end

#= The same question on the sphere, over the same three orientations.

`spherical_orient(a, b, c)` is `sign((a × b) ⋅ c)`: positive when `c` lies left of the
directed great-circle arc `a → b`. That is the same handedness `Predicates.orient` gives
in the plane, so the three signs combine by exactly the planar rule below and only the
predicate underneath changes.

`P1-P2-P3` here are always original ring vertices — `_get_sides` walks back to real
vertices before calling — and this is reached only for a hinge or overlap at `P2`, which
is where the chart edge of a DGG cell differs most from the great circle through its
endpoints. Classifying that hinge with the planar determinant is what made the crossing /
bouncing decision wrong for non-convex spherical cells. =#
function _get_side(::Spherical, Q, P1, P2, P3; exact)
    q = _spherical_kernel_point(Q)
    p1 = _spherical_kernel_point(P1)
    p2 = _spherical_kernel_point(P2)
    p3 = _spherical_kernel_point(P3)
    orient = _spherical_orient_for(booltype(exact))
    s1 = orient(q, p1, p2)
    s2 = orient(q, p2, p3)
    s3 = orient(p1, p2, p3)

    return _side_from_orientations(s1, s2, s3)
end

#= Which spherical orientation predicate `exact` selects, mirroring how the planar
`_get_side` threads `exact` into `Predicates.orient`.

`exact_spherical_orient` is a true sign function; `spherical_orient` reports `0` inside an
`eps*16` band, which at cell scale is wider than the determinant it is judging. Since the
hinge classified here is precisely where a DGG cell's chart edge departs furthest from the
great circle through its endpoints, a spurious `0` there is a wrong crossing/bouncing
decision, not a harmless tie. =#
@inline _spherical_orient_for(::True) = UnitSpherical.exact_spherical_orient
@inline _spherical_orient_for(::False) = UnitSpherical.spherical_orient

#= Reads the three orientations as a side. `s3` orients the hinge `P1-P2-P3` itself, and
`s1`/`s2` place `Q` against each of its legs: `Q` is inside the hinge's turn only when it
is on the turn's side of both, so a single disagreement puts it on the other side. =#
function _side_from_orientations(s1, s2, s3)
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
        _clip_midpoint(alg.manifold, pt_list[start_idx].point, pt_list[next_idx].point)
    end
    start_status = !_point_filled_curve_orientation(alg.manifold, start_pt, poly; in = true, on = false, out = false, exact)
    return next_idx, start_status
end

#= Whether `p2` carries no shape and may be dropped from a traced ring — that is, whether it
already lies on the edge joining its neighbours.

Which edge that is, is the whole question. A run of vertices along a parallel — the 49th
between Canada and the United States, lat 22 between Egypt and Sudan — is exactly collinear
in the chart, so the planar test drops every interior vertex of the run. On the sphere those
vertices are not redundant at all: the edge joining the ends of the run is a great-circle arc
that bulges poleward of the parallel, by 0.8° over a 28° span at latitude 49. Dropping them
therefore does not simplify the ring, it moves its boundary, and the sliver between the
polyline and the arc is lost from the result.

Asking `spherical_orient` instead asks whether `p2` lies on the great circle through `p1` and
`p3`, which is the edge the spherical clipper actually draws. Vertices along a parallel fail
that test and are kept; vertices genuinely on a shared great circle still go. The same
applies to the chart edges of a DGG cell, which are not great circles either. =#
_is_removable_collinear(::Planar, p1, p2, p3) =
    Predicates.orient(p1, p2, p3; exact = False()) == 0
_is_removable_collinear(::Spherical, p1, p2, p3) =
    UnitSpherical.spherical_orient(_spherical_kernel_point(p1),
        _spherical_kernel_point(p2), _spherical_kernel_point(p3)) == 0

#= A point strictly between two adjacent points of a traced ring, used to ask which side of
the other polygon the piece of boundary between them runs.

The question is only meaningful if the probe lies *on* the boundary piece it is standing in
for. In the plane the chart midpoint does. On the sphere it does not: the boundary is the
great-circle arc, and the chart midpoint sits off it, pulled toward the chord by the arc's
sagitta. Where the two polygons share a border — every interior edge of a tiling, and every
land border in a country dataset — that displacement is perpendicular to the very edge being
classified, so the in/out answer is decided by the sagitta rather than by the geometry, and
the entry/exit alternation it feeds stops alternating.

The spherical midpoint is the normalized sum of the two unit vectors, which is the
great-circle midpoint and needs no angle. `p + q` vanishing means the two are antipodal,
where no midpoint is defined and either of the two equidistant candidates would be a guess;
the chart midpoint is returned there so the caller still gets a point, and the antipodal
edge itself is what `antipodal_edge_split.jl` exists to remove upstream. =#
_clip_midpoint(::Planar, p, q) = (p .+ q) ./ 2
function _clip_midpoint(::Spherical, p, q)
    u = _spherical_kernel_point(p) + _spherical_kernel_point(q)
    n = norm(u)
    n == 0 && return _sph_mid_degenerate(p, q)
    return _sph_mid_as(UnitSphericalPoint(u ./ n), p)
end

#= The midpoint is fed straight back to a predicate alongside the nodes it came from, so it
has to speak their representation, not be round-tripped into the chart. =#
_sph_mid_as(mid, ::UnitSphericalPoint) = mid
_sph_mid_as(mid, _) = _usp_to_lonlat(mid)
#-- exactly antipodal: every great circle through the pair is a bisector, so no midpoint is
#-- more correct than another; take an endpoint and keep the return type stable.
_sph_mid_degenerate(p::UnitSphericalPoint, q) = p
_sph_mid_degenerate(p, q) = (p .+ q) ./ 2

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
                    next_val = _clip_midpoint(alg.manifold, curr_pt.point, pt_list[next_idx].point)
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
    status = !_point_filled_curve_orientation(alg.manifold, pt_list[1].point, poly; in = true, on = false, out = false, exact)
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

# Get type of polygons that will be made
# TODO: Increase type options
_get_poly_type(::Type{T}) where T = _get_poly_type(T, Tuple{T, T})
#-- the wrapper's `Z` flag has to agree with the point type, or the polygon the collector
#-- declares and the one `GI.Polygon` infers from 3D points are different types
_fh_pt_is3d(::Type{<:UnitSpherical.UnitSphericalPoint}) = true
_fh_pt_is3d(::Type) = false
_get_poly_type(::Type{T}, ::Type{P}) where {T, P} =
    GI.Polygon{_fh_pt_is3d(P), false,
        Vector{GI.LinearRing{_fh_pt_is3d(P), false, Vector{P}, Nothing, Nothing}}, Nothing, Nothing}

#= Egress mirrors ingress.

Planar hands back what it was given. Spherical computes in xyz, so the representation the
caller gets is the one it supplied: 3D input is returned untouched -- bit-exact, and already
the representation a DGG caller works in -- while lon/lat input converts back. A vertex that
passed through the clip unchanged is not converted at all on that path either: `srcidx`
names its slot in the input ring, so it is returned as the very value that came in, and only
genuinely computed intersections pay a conversion. =#
_fh_out_point_type(::Planar, poly, ::Type{T}) where {T} = Tuple{T, T}
_fh_out_point_type(::Spherical, poly, ::Type{T}) where {T} =
    GI.is3d(poly) ? UnitSpherical.UnitSphericalPoint{T} : Tuple{T, T}

#-- `nothing` means "the stored representation is already what the caller wants"
_fh_egress_ring(::Planar, poly) = nothing
_fh_egress_ring(::Spherical, poly) = GI.is3d(poly) ? nothing : _fh_source_ring(poly)

#-- The node lists are built from a ring, but the tracer is handed whatever the caller had,
#-- which for the polygon entry points is the polygon. `srcidx` indexes the ring.
_fh_source_ring(g) = _fh_source_ring(GI.trait(g), g)
_fh_source_ring(::GI.PolygonTrait, g) = GI.getexterior(g)
_fh_source_ring(::Any, g) = g

_fh_egress(node, ::Nothing, ::Type{T}) where {T} = node.point
_fh_egress(node, ring, ::Type{T}) where {T} =
    node.srcidx == 0 ? _sph_lonlat(T, node.point) :
                       _tuple_point(GI.getpoint(ring, Int(node.srcidx)), T)

#=
    abstract type _RingSink

What `_trace_polynodes!` does with the vertices it walks. The traversal is the same
whether the result rings are being built or only measured, so that difference is the one
thing the tracer takes as a parameter.

## Interface

A sink is fed one ring at a time, in traversal order, and carries whatever per-ring
working value it likes as `state` — the tracer creates it, threads it through, hands it
back, and never inspects it. Three methods, all required:

| method | returns | contract |
|:-------|:--------|:---------|
| `_ring_start(sink, pt)` | `state` | open a ring whose first vertex is `pt` |
| `_ring_step(sink, state, pt)` | `state` | extend the ring by `pt` |
| `_ring_close!(sink, state)` | `nothing` | fold the finished ring into `sink` |

The ring closes on the vertex it opened with: `_ring_step` is called with the first vertex
again before `_ring_close!`, so a sink that walks edges gets the closing edge for free and
one that collects points gets a closed ring.

A sink accumulates across rings and is read afterwards, so it is the mutable half of the
pair; `state` is per-ring and may be immutable. Two implementations ship:
`_RingCollector` (the result polygons) and `_RingMeasurer` (their total area).
=#
abstract type _RingSink end

#-- Interface fallbacks. Without them a sink missing a method fails as a `MethodError` on
#-- an internal call several frames into the tracer, which says nothing about what is
#-- actually missing.
_ring_start(sink::_RingSink, pt) = _ring_sink_incomplete(sink, "_ring_start(sink, pt)")
_ring_step(sink::_RingSink, state, pt) = _ring_sink_incomplete(sink, "_ring_step(sink, state, pt)")
_ring_close!(sink::_RingSink, state) = _ring_sink_incomplete(sink, "_ring_close!(sink, state)")

_ring_sink_incomplete(sink, sig) = throw(ArgumentError(
    "$(typeof(sink)) is a `_RingSink` but does not implement `$sig`. A ring sink must " *
    "implement `_ring_start`, `_ring_step` and `_ring_close!` — see the interface note " *
    "above `_RingSink` in clipping_processor.jl."))

# The result polygons: a point vector per ring, wrapped as a polygon. The original — and
# still the only — behaviour of `_trace_polynodes`.
struct _RingCollector{P} <: _RingSink
    polys::Vector{P}
end
_RingCollector(::Type{T}, ::Type{P} = Tuple{T, T}) where {T, P} =
    _RingCollector(Vector{_get_poly_type(T, P)}(undef, 0))

_ring_start(::_RingCollector, pt) = [pt]
_ring_step(::_RingCollector, pts, pt) = (push!(pts, pt); pts)
_ring_close!(sink::_RingCollector, pts) = (push!(sink.polys, GI.Polygon([pts])); nothing)

# The total area of those same rings, accumulated as they are walked. This is what lets
# `intersection_area` trace without materializing a ring at all.
mutable struct _RingMeasurer{M <: Manifold, T} <: _RingSink
    manifold::M
    area::T
    nrings::Int   # `isempty(polys)` for the collector: whether the trace found anything
end
_RingMeasurer(m::M, ::Type{T}) where {M, T} = _RingMeasurer{M, T}(m, zero(T), 0)

#-- state is (first vertex, previous vertex, running sum): both formulas below are
#-- two-point recurrences, so the ring never has to exist all at once
_ring_start(sink::_RingMeasurer{M, T}, pt) where {M, T} = (pt, pt, zero(T))
_ring_step(sink::_RingMeasurer, (first_pt, prev, acc), pt) =
    (first_pt, pt, acc + _ring_term(sink.manifold, first_pt, prev, pt))
function _ring_close!(sink::_RingMeasurer, (first_pt, prev, acc))
    #-- the tracer already closed the ring on `first_pt`, so `acc` is complete
    sink.area += abs(_ring_total(sink.manifold, acc))
    sink.nrings += 1
    return nothing
end

#-- the per-vertex terms of `_ring_area`'s two formulas (methods/area.jl), taken one
#-- vertex at a time. Summed in the same order, they give the same answer. The spherical
#-- one is untested: `FosterHormannClipping(Spherical())` is an ambiguous constructor call
#-- today, so no spherical FH algorithm can be built to reach it.
_ring_term(::Planar, first_pt, prev, pt) = _area_component(prev, pt)
#-- `_spherical_kernel_point` is the identity on a `UnitSphericalPoint`, which is what the
#-- tracer now carries, so this wrap costs nothing on the spherical path and still accepts
#-- lon/lat from any other caller.
_ring_term(::Spherical, first_pt, prev, pt) = _spherical_triangle_area(Eriksson(),
    _spherical_kernel_point(first_pt), _spherical_kernel_point(prev),
    _spherical_kernel_point(pt))
_ring_total(::Planar, acc) = acc / 2
_ring_total(::Spherical, acc) = acc

_trace_polynodes(alg::FosterHormannClipping, ::Type{T}, a_list, b_list, a_idx_list, f_step, poly_a, poly_b) where {T} =
    _trace_polynodes!(_RingCollector(T, _fh_out_point_type(alg.manifold, poly_a, T)),
        alg, T, a_list, b_list, a_idx_list, f_step, poly_a, poly_b).polys

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
function _trace_polynodes!(sink::_RingSink, alg::FosterHormannClipping{M, A}, ::Type{T}, a_list, b_list, a_idx_list, f_step, poly_a, poly_b) where {T, M, A}
    ring_a = _fh_egress_ring(alg.manifold, poly_a)
    ring_b = _fh_egress_ring(alg.manifold, poly_b)
    n_a_pts, n_b_pts = length(a_list), length(b_list)
    total_pts = n_a_pts + n_b_pts
    n_cross_pts = length(a_idx_list)
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
        ring = _ring_start(sink, _fh_egress(curr, ring_a, T))

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

                # Get current node and add to the ring
                curr = curr_list[idx]
                ring = _ring_step(sink, ring, _fh_egress(curr, on_a_list ? ring_a : ring_b, T))
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
        _ring_close!(sink, ring)
    end
    return sink
end

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
        _point_filled_curve_orientation(m, a_list[non_intr_a_idx].point, b_poly; exact)
    b_pt_orient = isnothing(non_intr_b_idx) ? point_on :
        _point_filled_curve_orientation(m, b_list[non_intr_b_idx].point, a_poly; exact)
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
                    if _is_removable_collinear(alg.manifold, p1, p2, p3)
                        remove_idx[i - 1] = true
                    end
                end
                p1, p2 = p2, p3
            end
            # Check if the first point (which is repeated as the last point) is needed 
            if _is_removable_collinear(alg.manifold, ring.geom[end - 1], ring.geom[1], ring.geom[2])
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
