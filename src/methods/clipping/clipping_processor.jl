# # Polygon clipping helpers
# This file contains the shared helper functions for the polygon clipping functionalities.

# This file specifically defines helpers for the Foster-Hormann clipping algorithm.


"""
    abstract type IntersectionAccelerator

Supertype for accelerators that reduce edge-pair intersection checks, with optional extra
memory.

`NestedLoop` takes O(n*m) time. `SingleSTRtree` indexes one ring in O(n*log(m)) query time.
`DoubleSTRtree` traverses two trees together.

`AutoAccelerator` selects an accelerator from the input polygon sizes in `build_a_list`.
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

Choose an accelerator from the input polygon sizes.
"""
struct AutoAccelerator <: IntersectionAccelerator end

"""
    FosterHormannClipping{M <: Manifold, A <: Union{Nothing, Accelerator}} <: GeometryOpsCore.Algorithm{M} 

Applies the Foster-Hormann clipping algorithm.

# Arguments
- `manifold::M`: The manifold on which the algorithm operates. `Geodesic` is not supported
  (the constructor throws); use [`Spherical`](@ref) instead.
- `accelerator::A`: The accelerator to use. `NestedLoop()` and `AutoAccelerator()` support
  spherical input; explicit tree accelerators currently require `Planar()`.

Spherical clipping preserves the first input's coordinate representation and requested numeric
type. Input arc predicates are exact; computed intersections use rounded coordinates.

Repeated clipping against the same boundary can produce inconsistent topology and raise
`TracingError`. Keep the default `fix_multipoly` correction for overlapping components.

Robust arbitrary chaining requires intersection provenance or exact constructions. This
implementation does not snap nearby vertices.
"""
struct FosterHormannClipping{M <: Manifold, A <: IntersectionAccelerator} <: GeometryOpsCore.Algorithm{M}
    manifold::M
    accelerator::A
    # TODO: add exact flag
    # TODO: should exact flag be in the type domain?
    #= There is no geodesic implementation of `_get_side` and the other clipping primitives,
    so a `Geodesic` algorithm would surface as a bare `MethodError` deep inside the first
    clip. Rejecting it here catches every construction path, parametric ones included. =#
    function FosterHormannClipping{M, A}(manifold::M, accelerator::A) where {M <: Manifold, A <: IntersectionAccelerator}
        manifold isa Geodesic && throw(ArgumentError(
            "FosterHormannClipping does not support the Geodesic manifold ($manifold): Foster-Hormann clipping has no geodesic implementation of its intersection primitives. Use Spherical() instead."
        ))
        return new{M, A}(manifold, accelerator)
    end
end
#= The inner constructor above suppresses Julia's automatic outer constructor, which every
other constructor below builds through, so it has to be written back explicitly. =#
FosterHormannClipping(manifold::M, accelerator::A) where {M <: Manifold, A <: IntersectionAccelerator} = FosterHormannClipping{M, A}(manifold, accelerator)
FosterHormannClipping(; manifold::Manifold = Planar(), accelerator = nothing) = FosterHormannClipping(manifold, isnothing(accelerator) ? NestedLoop() : accelerator)
FosterHormannClipping(manifold::Manifold, accelerator::Union{Nothing, IntersectionAccelerator} = nothing) = FosterHormannClipping(manifold, isnothing(accelerator) ? NestedLoop() : accelerator)
FosterHormannClipping(accelerator::Union{Nothing, IntersectionAccelerator}) = FosterHormannClipping(Planar(), isnothing(accelerator) ? NestedLoop() : accelerator)
#= Spherical clipping uses `NestedLoop` because tree bounds are planar rectangles.
Keep both argument types narrow: `Union{Nothing, IntersectionAccelerator}` makes
spherical constructor dispatch ambiguous. =#
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
    #= 1-based position in the ring this vertex was ingested from, or 0 for a computed
    intersection. `Spherical` converts to xyz at ingress, so this is what lets egress hand a
    passthrough vertex back exactly as it arrived instead of round-tripping it. =#
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
    FosterHormannCache(alg::FosterHormannClipping, [T = Float64])
    FosterHormannCache(m::Manifold, [T = Float64])
    FosterHormannCache([T = Float64])

Buffers for [`FosterHormannClipping`](@ref): ring vertices, intersection indices, and
converted spherical inner edges.

Pass `cache` to [`intersection_area`](@ref) to reuse these buffers. Results do not reference
the buffers and remain valid after cache reuse.

The cache must match the numeric type and manifold point representation. Use
`FosterHormannCache(alg, T)` to match both; `FosterHormannCache(T)` is planar.

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
    b_edges::Vector{Tuple{P,P}}
end
FosterHormannCache{T, P}() where {T, P} =
    FosterHormannCache{T, P}(PolyNode{T, P}[], PolyNode{T, P}[], Int[], Tuple{P,P}[])

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
identities on points already in that representation, so 3D input reaches the node lists
without arithmetic and comes back out bit-identical. =#
_fh_ingest(::Planar, p, ::Type{T}) where {T} = _tuple_point(p, T)
_fh_ingest(::Spherical, p, ::Type{T}) where {T} = _spherical_edge_point(p, T)

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

# Store the polygons and node lists when tracing fails.
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

Build both rings as `PolyNode` vectors and set their entry/exit flags. Return `(a_list,
b_list, a_idx_list)`, where `a_idx_list[i]` locates intersection `i` in `a_list`.
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

_reusable_inner_edges(m::Manifold, geom, ::Type{T}) where {T} = eachedge(m, geom, T)
_reusable_inner_edges(m::Spherical, geom, ::Type{T}) where {T} =
    GI.is3d(geom) ? eachedge(m, geom, T) : collect(eachedge(m, geom, T))

_check_planar_edge_accelerator(::Planar, accelerator) = nothing
_check_planar_edge_accelerator(m::Manifold, accelerator) = throw(ArgumentError(
    "$(typeof(accelerator)) indexes planar edge extents and does not support $m. Use NestedLoop() or AutoAccelerator() for spherical clipping."))

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

`accelerator` reduces candidate edge pairs while preserving this callback order.
[`AutoAccelerator`](@ref) selects a method by an internal heuristic; `SingleSTRtree` uses a
tree and extent filtering.

"""
function foreach_pair_of_maybe_intersecting_edges_in_order(
    manifold::M, accelerator::AutoAccelerator, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64
) where {FA, FAAfter, FI, T, M <: Manifold}
    if manifold isa Spherical
        return foreach_pair_of_maybe_intersecting_edges_in_order(manifold, NestedLoop(), f_on_each_a, f_after_each_a, f_on_each_maybe_intersect, poly_a, poly_b, T)
    end
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
    manifold::M, accelerator::NestedLoop, f_on_each_a::FA, f_after_each_a::FAAfter, f_on_each_maybe_intersect::FI, poly_a, poly_b, _t::Type{T} = Float64; inner_edges = nothing
) where {FA, FAAfter, FI, T, M <: Manifold}
    # this is suitable for planar
    # but spherical / geodesic will need s2 support at some point,
    # or -- even now -- just buffering
    na = GI.npoint(poly_a)
    nb = GI.npoint(poly_b)
    # Use nested loops for small polygons on any manifold. Convert spherical inner edges once
    # to avoid repeated longitude/latitude conversions.
    edges_b = if inner_edges === nothing || GI.is3d(poly_b)
        _reusable_inner_edges(manifold, poly_b, T)
    else
        empty!(inner_edges)
        append!(inner_edges, eachedge(manifold, poly_b, T))
    end
    # First, loop over "each edge" in poly_a
    for (i, (a1t, a2t)) in enumerate(eachedge(manifold, poly_a, T))
        a1t == a2t && continue
        isnothing(f_on_each_a) || f_on_each_a(a1t, i)
        for (j, (b1t, b2t)) in enumerate(edges_b)
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
    _check_planar_edge_accelerator(manifold, accelerator)
    na = GI.npoint(poly_a)
    nb = GI.npoint(poly_b)
    # Index only `poly_b` to avoid constructing an edge list and tree for `poly_a`.
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
                # Handle `LoopStateMachine.Action` results so callbacks can control the loop.
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
    _check_planar_edge_accelerator(manifold, accelerator)
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
    _check_planar_edge_accelerator(manifold, accelerator)
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
    _check_planar_edge_accelerator(manifold, accelerator)
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

Build `a_list` from `poly_a` vertices and intersections with `poly_b`. Neighbor indices into
`b_list` and entry/exit flags remain unset.

`a_idx_list[i]` is the index of intersection `i` in `a_list`.
=#
function _build_a_list(alg::FosterHormannClipping{M, A}, ::Type{T}, poly_a, poly_b; exact, cache = nothing) where {T, M, A}
    n_a_edges = _nedge(poly_a)
    # list of points in poly_a
    P = _fh_point_type(alg.manifold, T)
    a_list = _fh_buffer(cache === nothing ? nothing : cache.a_list, Vector{PolyNode{T, P}})
    #-- Skip `sizehint!` on cached buffers because it can shrink retained capacity.
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
    if cache !== nothing && alg.manifold isa Spherical && alg.accelerator isa NestedLoop
        foreach_pair_of_maybe_intersecting_edges_in_order(
            alg.manifold, alg.accelerator, on_each_a, after_each_a, on_each_maybe_intersect,
            poly_a, poly_b, T; inner_edges = cache.b_edges,
        )
    else
        foreach_pair_of_maybe_intersecting_edges_in_order(alg, on_each_a, after_each_a, on_each_maybe_intersect, poly_a, poly_b, T)
    end

    return a_list, a_idx_list, n_b_intrs
end

#=
    _build_b_list(::Type{T}, a_idx_list, a_list, poly_b) -> b_list

Build `b_list` from `poly_b` and the intersections in `a_list`. Update neighbor indices in
`a_list`; entry/exit flags remain unset.
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

Classify intersections as crossing or bouncing. For overlapping chains, the outer edges
determine the chain classification.

Mark the first and last points as crossing for a crossing chain, or delayed otherwise. Mark
middle points as bouncing and both ends in the `endpoints` field.
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

#= Classify `b_prev` and `b_next` against the hinge `a_prev-curr_pt-a_next`.
For hinges and overlaps, `curr_pt` is an original vertex. Use the nearest original
vertices for orientation to avoid errors from computed intersection coordinates. =#
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

#= Classify spherical hinges with `sign((a × b) ⋅ c)`. Its handedness matches
`Predicates.orient`, so the three signs combine with the planar side rule. =#
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

#= Select exact signs or the `eps*16` tolerance band from `exact`. The band can
classify cell-scale crossings as collinear. =#
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

#= Use a non-intersection vertex to determine the next intersection's entry/exit flag.
If none exists, probe after a chain end or an unchained crossing. Return `nothing`
as the next index if no such point exists. =#
function _pt_off_edge_status(m::Manifold, pt_list, poly, npts; exact)
    start_idx, is_non_intr_pt = findfirst(_is_not_intr, pt_list), true
    if isnothing(start_idx)
        start_idx, is_non_intr_pt = findfirst(_next_edge_off, pt_list), false
        isnothing(start_idx) && return (start_idx, false)
    end
    next_idx = start_idx < npts ? (start_idx + 1) : 1
    start_pt = if is_non_intr_pt
        pt_list[start_idx].point
    else
        _clip_midpoint(m, pt_list[start_idx].point, pt_list[next_idx].point)
    end
    start_status = !_point_filled_curve_orientation(m, start_pt, poly; in = true, on = false, out = false, exact)
    return next_idx, start_status
end

#= Drop `p2` only if it lies on the manifold edge joining its neighbors.
Vertices along a latitude parallel generally do not lie on that great-circle arc. =#
_is_removable_collinear(::Planar, p1, p2, p3) =
    Predicates.orient(p1, p2, p3; exact = False()) == 0
_is_removable_collinear(::Spherical, p1, p2, p3) =
    UnitSpherical.spherical_orient(_spherical_kernel_point(p1),
        _spherical_kernel_point(p2), _spherical_kernel_point(p3)) == 0

#= Probe the traced boundary at its great-circle midpoint, the normalized endpoint sum.
A chart midpoint lies off the arc and can misclassify shared boundaries.

Antipodal endpoints have a zero sum; upstream `antipodal_edge_split.jl` removes them.
The degenerate branch preserves the return type. =#
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
#-- exactly antipodal: no candidate midpoint is more correct than another, so pick one and
#-- keep the return type stable.
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

Flag intersections as entry or exit. Ordinary crossings alternate; delayed bounces have
opposite endpoint flags, and delayed crossings have equal endpoint flags.

Operation-specific callbacks update crossing/bouncing classifications. Delayed crossings have
different endpoint classifications; delayed bounces have equal classifications.
=#
function _flag_ent_exit!(alg::FosterHormannClipping{M, A}, ::Type{T}, ::GI.LinearRingTrait, poly, pt_list, delay_cross_f, delay_bounce_f; exact) where {T, M, A}
    npts = length(pt_list)
    # Find starting index if there is one
    next_idx, status = _pt_off_edge_status(alg.manifold, pt_list, poly, npts; exact)
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

# GeoInterface's convenience constructor inspects the first polygon even when Z/M
# are explicit. Supply the concrete wrapper type so empty clipping results are valid.
function _fh_multipolygon(polys::Vector{P}; crs = nothing) where {P}
    Z = _fh_pt_is3d(_fh_poly_point_type(P))
    return GI.MultiPolygon{Z, false, typeof(polys), Nothing, typeof(crs)}(polys, nothing, crs)
end

#= Return vertices in the input representation. Fetch unchanged vertices by `srcidx`
to preserve their exact values; convert only computed intersections. =#
_fh_out_point_type(::Planar, poly, ::Type{T}) where {T} = Tuple{T, T}
_fh_out_point_type(::Spherical, poly, ::Type{T}) where {T} =
    GI.is3d(poly) ? UnitSpherical.UnitSphericalPoint{T} : Tuple{T, T}

#= No-crossing results must use the tracer's output representation. These helpers
preserve matching points exactly and convert other representations. =#
_fh_as_point(::Type{<:Tuple}, p, ::Type{T}) where {T} = _fh_tuple_point(p, T)
_fh_as_point(::Type{<:UnitSpherical.UnitSphericalPoint}, p, ::Type{T}) where {T} =
    _spherical_edge_point(p, T)

#-- `_tuple_point` is the identity on a `UnitSphericalPoint`, which is right for the node
#-- lists and wrong here, where a tuple is what was asked for.
_fh_tuple_point(p::UnitSpherical.UnitSphericalPoint, ::Type{T}) where {T} = _sph_lonlat(T, p)
_fh_tuple_point(p, ::Type{T}) where {T} = _tuple_point(p, T)

_fh_as_ring(::Type{P}, ring, ::Type{T}) where {P, T} =
    GI.LinearRing([_fh_as_point(P, p, T) for p in GI.getpoint(ring)])
_fh_as_poly(::Type{P}, poly, ::Type{T}) where {P, T} =
    GI.Polygon([_fh_as_ring(P, r, T) for r in GI.getring(poly)])

#-- for the callers that know the point type but do not carry `T`
_fh_float_type(::Type{<:Tuple{T, T}}) where {T} = T
_fh_float_type(::Type{<:UnitSpherical.UnitSphericalPoint{T}}) where {T} = T
_fh_as_ring(::Type{P}, ring) where {P} = _fh_as_ring(P, ring, _fh_float_type(P))
_fh_as_poly(::Type{P}, poly) where {P} = _fh_as_poly(P, poly, _fh_float_type(P))

#-- Recover the point representation an already-built output vector is committed to, for the
#-- helpers that are handed `polys` rather than the original geometries.
_fh_ring_point_type(::Type{<:GI.LinearRing{Z, M, V}}) where {Z, M, V} = eltype(V)
_fh_poly_point_type(::Type{<:GI.Polygon{Z, M, R}}) where {Z, M, R} = _fh_ring_point_type(eltype(R))

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

Consume traced ring vertices to construct polygons or measure their area.

## Interface

The tracer passes each ring in order and carries an opaque per-ring `state`. All three methods
are required:

| method | returns | contract |
|:-------|:--------|:---------|
| `_ring_start(sink, pt)` | `state` | open a ring whose first vertex is `pt` |
| `_ring_step(sink, state, pt)` | `state` | extend the ring by `pt` |
| `_ring_close!(sink, state)` | `nothing` | fold the finished ring into `sink` |

Before closing, `_ring_step` receives the first vertex again, supplying the closing edge.

The mutable sink accumulates across rings; per-ring `state` may be immutable. `_RingCollector`
builds polygons, and `_RingMeasurer` sums their areas.
=#
abstract type _RingSink end

#-- Report missing sink methods at the interface boundary.
_ring_start(sink::_RingSink, pt) = _ring_sink_incomplete(sink, "_ring_start(sink, pt)")
_ring_step(sink::_RingSink, state, pt) = _ring_sink_incomplete(sink, "_ring_step(sink, state, pt)")
_ring_close!(sink::_RingSink, state) = _ring_sink_incomplete(sink, "_ring_close!(sink, state)")

_ring_sink_incomplete(sink, sig) = throw(ArgumentError(
    "$(typeof(sink)) is a `_RingSink` but does not implement `$sig`. A ring sink must " *
    "implement `_ring_start`, `_ring_step` and `_ring_close!` — see the interface note " *
    "above `_RingSink` in clipping_processor.jl."))

# Collect each traced ring into a point vector and wrap it as a polygon.
struct _RingCollector{P} <: _RingSink
    polys::Vector{P}
end
_RingCollector(::Type{T}, ::Type{P} = Tuple{T, T}) where {T, P} =
    _RingCollector(Vector{_get_poly_type(T, P)}(undef, 0))

function _ring_start(::_RingCollector{Poly}, pt) where {Poly}
    P = _fh_poly_point_type(Poly)
    return P[_fh_as_point(P, pt, _fh_float_type(P))]
end
_ring_step(::_RingCollector, pts::Vector{P}, pt) where {P} =
    (push!(pts, _fh_as_point(P, pt, _fh_float_type(P))); pts)
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
#-- vertex at a time. Summed in the same order, they give the same answer.
_ring_term(::Planar, first_pt, prev, pt) = _area_component(prev, pt)
#-- Use the same canonical unit-sphere points as the public spherical area calculation.
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

Trace `_build_ab_list` outputs into GeoInterface polygons using Greiner-Hormann traversal.
`f_step(entry, in_a)` selects the direction from entry/exit status and the active list:
    - Intersection: (x, y) -> x ? 1 : (-1)
    - Difference: (x, y) -> (x ⊻ y) ? 1 : (-1)
    - Union: (x, y) -> x ? (-1) : 1
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

Return whether each polygon lies inside the other when no intersections cross. Shared edges
and points are allowed; edge probes distinguish containment from disjoint interiors.
=#
function _find_non_cross_orientation(m::M, a_list, b_list, a_poly, b_poly; exact) where {M <: Manifold}
    # Shared vertices do not imply shared edges: an enclave can leave the other
    # boundary between two shared vertices. Probe an edge leaving a shared chain
    # when there is no non-intersection vertex, as entry/exit classification does.
    a_idx, a_out = _pt_off_edge_status(m, a_list, b_poly, length(a_list); exact)
    b_idx, b_out = _pt_off_edge_status(m, b_list, a_poly, length(b_list); exact)
    a_pt_orient = isnothing(a_idx) ? point_on : a_out ? point_out : point_in
    b_pt_orient = isnothing(b_idx) ? point_on : b_out ? point_out : point_in
    a_in_b = a_pt_orient != point_out && b_pt_orient != point_in
    b_in_a = b_pt_orient != point_out && a_pt_orient != point_in
    return a_in_b, b_in_a
end

_find_non_cross_orientation(alg::FosterHormannClipping{M}, a_list, b_list, a_poly, b_poly; exact) where {M <: Manifold} =
    _find_non_cross_orientation(alg.manifold, a_list, b_list, a_poly, b_poly; exact)

#=
    _add_holes_to_polys!(::Type{T}, return_polys, hole_iterator, remove_poly_idx; exact)

Subtract `hole_iterator` from `return_polys`. Append split pieces and remove fully covered
polygons.
=#
function _add_holes_to_polys!(alg::FosterHormannClipping{M, A}, ::Type{T}, return_polys, hole_iterator, remove_poly_idx; exact) where {T, M, A}
    n_polys = length(return_polys)
    remove_hole_idx = Int[]
    # Remove set of holes from all polygons
    for i in 1:n_polys
        n_new_per_poly = 0
        for curr_hole in Iterators.map(h -> _fh_as_ring(_fh_poly_point_type(eltype(return_polys)), h, T), hole_iterator) # loop through all holes
            curr_hole = _linearring(curr_hole)
            # loop through all pieces of original polygon (new pieces added to end of list)
            for j in Iterators.flatten((i:i, (n_polys + 1):(n_polys + n_new_per_poly)))
                curr_poly = return_polys[j]
                remove_poly_idx[j] && continue
                curr_poly_ext = GI.nhole(curr_poly) > 0 ? GI.Polygon(StaticArrays.SVector(GI.getexterior(curr_poly))) : curr_poly
                in_ext, on_ext, out_ext = _line_polygon_interactions(alg.manifold, curr_hole, curr_poly_ext; exact, closed_line = true)
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
                elseif coveredby(alg.manifold, curr_poly_ext, GI.Polygon(StaticArrays.SVector(curr_hole)))
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

Merge `new_hole` with intersecting holes in `curr_poly` and remove those holes. If their union
encloses an island, append it to `return_polys` and reassign remaining holes.

Return `new_hole` unchanged when no existing hole touches it.
=#
function _combine_holes!(alg::FosterHormannClipping{M, A}, ::Type{T}, new_hole, curr_poly, return_polys, remove_hole_idx) where {T, M, A}
    n_new_polys = 0
    empty!(remove_hole_idx)
    new_hole_poly = GI.Polygon(StaticArrays.SVector(new_hole))
    # Combine any existing holes in curr_poly with new hole
    for (k, old_hole) in enumerate(GI.gethole(curr_poly))
        old_hole_poly = GI.Polygon(StaticArrays.SVector(old_hole))
        if intersects(alg.manifold, new_hole_poly, old_hole_poly)
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
            if !(k in remove_hole_idx) && within(alg.manifold, old_hole, piece)
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
