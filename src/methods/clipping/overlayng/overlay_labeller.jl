# NOTE: This functionality is experimental and may change at any time.

# # OverlayLabeller — the five-pass edge labelling (port of JTS `OverlayLabeller`)
#
# Phase 2b of the OverlayNG port (design doc `2026-07-16-overlayng-noding-substrate.md`,
# §3). A faithful port of `operation/overlayng/OverlayLabeller.java`, operating on
# the phase-2a index-based `OverlayGraph`. The single shared `OverlayLabel` between a
# symmetric pair means every `set_location_*!` here is seen by both half-edges —
# JTS relies on this, and so does this port.
#
# Also hosts the overlay op-code enum and `_is_result_of_op` (JTS
# `OverlayNG.isResultOfOp`, the boolean core used by both the labeller and the line
# builder), and the topology-error type (JTS `TopologyException`).
#
# Everything here is internal to GeometryOps — nothing is exported.

# ## Overlay operation codes and the op boolean core

# The four set-theoretic overlay operations (port of the `OverlayNG.INTERSECTION`
# / `UNION` / `DIFFERENCE` / `SYMDIFFERENCE` op codes).
@enum _OverlayOpCode::UInt8 begin
    OVERLAY_INTERSECTION
    OVERLAY_UNION
    OVERLAY_DIFFERENCE
    OVERLAY_SYMDIFFERENCE
end

"""
    _OverlayTopologyError(msg)

A robustness/topology error raised by the overlay engine (port of JTS
`TopologyException`). Signals an inconsistency the builder could not resolve
(e.g. a side-location conflict during area propagation, or a ring that cannot
be closed).
"""
struct _OverlayTopologyError <: Exception
    msg::String
end
Base.showerror(io::IO, e::_OverlayTopologyError) = print(io, "OverlayTopologyError: ", e.msg)

# Port of `OverlayNG.isResultOfOp`: whether a point with the given per-input
# locations lies in the result of `op`. `LOC_BOUNDARY` counts as `LOC_INTERIOR`.
@inline function _is_result_of_op(op::_OverlayOpCode, loc0::Integer, loc1::Integer)
    loc0 == LOC_BOUNDARY && (loc0 = LOC_INTERIOR)
    loc1 == LOC_BOUNDARY && (loc1 = LOC_INTERIOR)
    if op == OVERLAY_INTERSECTION
        return loc0 == LOC_INTERIOR && loc1 == LOC_INTERIOR
    elseif op == OVERLAY_UNION
        return loc0 == LOC_INTERIOR || loc1 == LOC_INTERIOR
    elseif op == OVERLAY_DIFFERENCE
        return loc0 == LOC_INTERIOR && loc1 != LOC_INTERIOR
    else # OVERLAY_SYMDIFFERENCE
        return (loc0 == LOC_INTERIOR && loc1 != LOC_INTERIOR) ||
               (loc0 != LOC_INTERIOR && loc1 == LOC_INTERIOR)
    end
end

# The `op` argument was always a predicate over the two per-input locations:
# the four set-op codes select the four hard-wired predicates above, and any
# callable `(loc0, loc1) -> Bool` may stand in their place (`LOC_BOUNDARY` is
# collapsed to `LOC_INTERIOR` before the call, exactly as for the set ops).
# This method plus the loosened `op` signatures below are the whole surface a
# custom result semantics needs on the dissolving (result-area) pipeline;
# non-dissolving consumers use the face enumeration in polygon_builder.jl.
@inline function _is_result_of_op(pred::F, loc0::Integer, loc1::Integer) where {F}
    loc0 == LOC_BOUNDARY && (loc0 = LOC_INTERIOR)
    loc1 == LOC_BOUNDARY && (loc1 = LOC_INTERIOR)
    return pred(loc0, loc1)::Bool
end

# ## computeLabelling — the five passes, in order (port of `computeLabelling`)

#=
Pass 1 is the one pass that reads a node's star as a CLOSED cycle — it walks it
CCW carrying a side location and checks every boundary edge against it. Clip
pruning (split.jl) can leave a node with only part of its star, which would make
that walk report a phantom side-location conflict, so those nodes are skipped.

Skipping is information-losing, never wrong: the edges pass 1 would have labelled
stay unknown and are resolved by passes 2–5, the last of which locates against
the ORIGINAL input by point-in-area and is total. `arr.truncated` is empty on
every unclipped arrangement, so the guard is inert on the default path.

It also cannot cost anything that matters. A truncated node lies strictly outside
the clip box, and every node inside the box keeps its full star (split.jl); a
node outside the box carries edges from ONE input only, because a node shared by
both inputs lies on a segment of each, hence inside both envelopes, hence inside
the box. So the star pass 1 is denied was never going to relate the two inputs to
each other.
=#
function _compute_labelling!(g::OverlayGraph, input)
    edges = g.edges
    truncated = g.arr.truncated
    any_truncated = !isempty(truncated)
    for ne in graph_node_edges(g)
        (any_truncated && truncated[he_origin(edges, ne)]) && continue
        _propagate_area_locations!(edges, input, ne, 0)
        _input_has_edges(input, 1) && _propagate_area_locations!(edges, input, ne, 1)
    end
    _label_connected_linear_edges!(edges, input)
    _label_collapsed_edges!(edges)
    _label_connected_linear_edges!(edges, input)
    _label_disconnected_edges!(g, input)
    return nothing
end

# ### Pass 1: area-node propagation (port of `propagateAreaLocations`)

# Scans a node's star CCW, propagating side labels for one area input to every
# edge whose location for that input is still unknown. A side-location conflict
# between two boundary edges is a topology error (port of the JTS check).
function _propagate_area_locations!(edges, input, node_edge::Integer, gi::Integer)
    _input_is_area(input, gi) || return nothing
    #-- one-edge node: nothing to propagate (dangling edge)
    he_degree(edges, node_edge) == 1 && return nothing

    e_start = _find_propagation_start_edge(edges, node_edge, gi)
    e_start == 0 && return nothing

    curr_loc = oe_get_location(edges, e_start, gi, POS_LEFT)
    e = he_onext(edges, e_start)
    while e != e_start
        label = oe_label(edges, e)
        if !is_boundary(label, gi)
            #-- non-boundary edge: its location relative to this area is now known
            set_location_line!(label, gi, curr_loc)
        else
            loc_right = oe_get_location(edges, e, gi, POS_RIGHT)
            loc_right == curr_loc ||
                throw(_OverlayTopologyError("side location conflict: arg $gi"))
            loc_left = oe_get_location(edges, e, gi, POS_LEFT)
            #-- loc_left == LOC_NONE should never happen for a boundary edge
            curr_loc = loc_left
        end
        e = he_onext(edges, e)
    end
    return nothing
end

# Port of `findPropagationStartEdge`: a boundary edge for `gi` in the node's
# star, or `0` if none.
function _find_propagation_start_edge(edges, node_edge::Integer, gi::Integer)
    e = Int32(node_edge)
    while true
        is_boundary(oe_label(edges, e), gi) && return e
        e = he_onext(edges, e)
        e == node_edge && break
    end
    return Int32(0)
end

# ### Pass 2/4: connected linear propagation (port of `labelConnectedLinearEdges`)

function _label_connected_linear_edges!(edges, input)
    _propagate_linear_locations!(edges, input, 0)
    _input_has_edges(input, 1) && _propagate_linear_locations!(edges, input, 1)
    return nothing
end

# BFS over linear (line or collapse) edges with a known location, propagating
# that location to connected unknown edges (port of `propagateLinearLocations`).
function _propagate_linear_locations!(edges, input, gi::Integer)
    stack = Int32[]
    for i in eachindex(edges)
        lbl = oe_label(edges, i)
        if is_linear(lbl, gi) && !is_line_location_unknown(lbl, gi)
            push!(stack, Int32(i))
        end
    end
    isempty(stack) && return nothing

    is_input_line = _input_is_line(input, gi)
    while !isempty(stack)
        e_node = pop!(stack)
        _propagate_linear_at_node!(edges, e_node, gi, is_input_line, stack)
    end
    return nothing
end

# Port of `propagateLinearLocationAtNode`. Line parents propagate EXTERIOR only.
function _propagate_linear_at_node!(edges, e_node::Integer, gi::Integer,
        is_input_line::Bool, stack::Vector{Int32})
    line_loc = get_line_location(oe_label(edges, e_node), gi)
    #-- a Line parent only propagates EXTERIOR locations
    is_input_line && line_loc != LOC_EXTERIOR && return nothing

    e = he_onext(edges, e_node)
    while e != e_node
        label = oe_label(edges, e)
        if is_line_location_unknown(label, gi)
            set_location_line!(label, gi, line_loc)
            #-- continue the traversal from the far node (don't re-add e itself)
            push!(stack, he_sym(edges, e))
        end
        e = he_onext(edges, e)
    end
    return nothing
end

# ### Pass 3: collapsed-edge ring-role labelling (port of `labelCollapsedEdges`)

function _label_collapsed_edges!(edges)
    for i in eachindex(edges)
        label = oe_label(edges, i)
        is_line_location_unknown(label, 0) && _label_collapsed_edge!(label, 0)
        is_line_location_unknown(label, 1) && _label_collapsed_edge!(label, 1)
    end
    return nothing
end

function _label_collapsed_edge!(label::OverlayLabel, gi::Integer)
    is_collapse(label, gi) || return nothing
    #-- disconnected collapsed edge: label from its parent ring role (shell/hole)
    set_location_collapse!(label, gi)
    return nothing
end

# ### Pass 5: disconnected-edge labelling (port of `labelDisconnectedEdges`,
# ### with the two robustness divergences F1 and F2)

#=
JTS's `labelDisconnectedEdges` asks the point-in-area locator, once per edge and
once per input, where that edge lies relative to the OTHER input's area, and
takes the answer as final. Two things about that are wrong on a near-coincident
arrangement, and this port fixes both. The two fixes are independent — the second
is not a consequence of the first (measured: a case that throws with ZERO
emitted-vs-exact disagreements in the first).

**F1 — the query point is the node's kernel point, not its emitted coordinate.**
JTS has one coordinate per node, so the distinction does not exist there. Here it
does, and design §0 is explicit that no decision may consume a constructed
coordinate. `node_point` is the *emission* — the single lossy step of the whole
substrate — and on `Spherical` it is doubly lossy for this use: the kernel unit
vector is turned into lon/lat by `atan`/`asin` here, and the locator turns it
straight back into a unit vector by `cos`/`sin`. Measured on real reprojected
data, that round trip moves a pass-through vertex by up to 14 ulps, which is an
order of magnitude more than the ±1 ulp separations that produce slivers in the
first place. Worse, the rounding is not merely noisy but *biased towards the
wrong answer*: an emitted coordinate that lands bit-identically on an input
vertex makes the locator report `LOC_BOUNDARY`, which the `!= LOC_EXTERIOR` test
below reads as INTERIOR. The node's own kernel point is exact, is what every
other decision in the engine already uses, and is right there in the node table.

**F2 — a known area location propagates; PIP only seeds where it cannot reach.**
See `_label_disconnected_area!`.
=#

function _label_disconnected_edges!(g::OverlayGraph, input)
    edges = g.edges
    for gi in (0, 1)
        if _input_is_area(input, gi)
            _label_disconnected_area!(g, input, gi)
        else
            #-- non-area target: a disconnected edge must be EXTERIOR (JTS
            #-- `labelDisconnectedEdge`'s `isArea` arm, unchanged)
            for i in eachindex(edges)
                label = oe_label(edges, i)
                is_line_location_unknown(label, gi) &&
                    set_location_all!(label, gi, LOC_EXTERIOR)
            end
        end
    end
    return nothing
end

# ### F1: the kernel point of a node, for the point-in-area query
#
# The node's own exact point in kernel space, never the emitted output
# coordinate. Manifold-generic: dispatched on the node key's point type, which
# IS the manifold's kernel point type.

@inline _node_kernel_point(arr::NodedArrangement, id::Integer) =
    _node_kernel_point(arr, (@inbounds arr.nodes.keys[id]), id)

@inline _node_kernel_point(arr::NodedArrangement, k::NodeKey, id::Integer) =
    k.is_crossing ? _crossing_kernel_point(arr, k, id) : k.pt

# Planar crossing: `RayCrossingCounter` stores its query point as a
# `Tuple{Float64,Float64}`, so the exact `Rational{BigInt}` crossing cannot be
# handed to it at all. The emitted coordinate is the *certified* correctly-
# rounded image of exactly that rational (emit.jl: accepted only when the
# residual plus the dd error bound is below ½ ulp), so it is the best Float64
# representation of the exact point that exists, and this arm is a no-op change.
# Planar vertex nodes take `k.pt` above, which is bit-identical to `node_point`.
# Hence: no planar behaviour changes here at all.
@inline _crossing_kernel_point(arr::NodedArrangement, k::NodeKey{Tuple{Float64, Float64}},
        id::Integer) = node_point(arr, id)

# Spherical crossing: the exact on-arc crossing direction (`_sph_crossing_dir`
# on the `True()` rational path — the same authority `rk_nodes_coincide` and the
# emission fallback use), scaled by its largest component before conversion so a
# large-numerator rational cannot overflow to `Inf`, then normalized. One
# rounding of a direction, against emission's `atan`/`asin` followed by the
# locator's `cos`/`sin`.
function _crossing_kernel_point(::NodedArrangement, k::NodeKey{<:UnitSphericalPoint},
        ::Integer)
    d = _sph_crossing_dir(True(), k)
    scale = max(abs(d[1]), abs(d[2]), abs(d[3]))
    x = Float64(d[1] / scale); y = Float64(d[2] / scale); z = Float64(d[3] / scale)
    s = sqrt(x * x + y * y + z * z)
    return UnitSphericalPoint(x / s, y / s, z / s)
end

# ### F2: propagate a known area location; seed by PIP only where it cannot reach

#=
The invariant the ring builder needs is that at every node the number of
in-result-area half-edges *leaving* equals the number *arriving*. Walk a node's
star and that is equivalent to the sectors between consecutive edges having a
consistent in-result flag: an outgoing half-edge is kept in the result area iff
the sector on its right is in the result and the one on its left is not, its sym
iff the reverse, so around a closed cycle the two counts are the number of
"enters" and the number of "leaves" of the same set of sectors, which are equal
— PROVIDED adjacent edges agree about the sector they share.

Pass 1 establishes that agreement for input `gi` at every node that has a
`gi`-boundary edge: it walks the star as a cycle carrying one side location. At a
node with NO `gi`-boundary edge there is nothing to walk, and JTS falls back to
asking the locator separately for each incident edge. At sub-ulp separations
those independent answers are effectively coin flips, and two of them disagreeing
across one degree-2 node is exactly an unbalanced node — measured on the real
reprojected-watershed reproducer: 28 unbalanced nodes out of 2,936.

The divergence from JTS is to notice that the answers are not independent. A node
carrying no boundary edge of `gi` is a node at which `gi`'s area location cannot
change: crossing it, one stays on the same side of `gi`'s boundary, because
`gi`'s boundary is not there. (A `gi`-*collapse* edge at the node is `gi` linework
with zero depth delta, so it does not separate either — it does not block
propagation. It is not used as a seed, though: pass 3 gives a collapse its
parent's ring role, which is a boundary marker, not the location of a
neighbouring face.) So every edge incident on such a node has the SAME `gi`
location, and one known value determines all of them — exactly, not
approximately. Propagating it is therefore not a heuristic that trades accuracy
for consistency; it is reading a value the arrangement already fixed.

That leaves only the chains that no known value reaches. Those get ONE
point-in-area verdict for the whole chain, so a chain can no longer land on two
sides at once, which is what broke the builder.

The verdict is a majority over the chain's EDGE INTERIORS, and both halves of
that — edge interiors rather than nodes, majority rather than JTS's unanimity —
are answering the same failure.

JTS asks `locateEdgeBothEnds` (interior iff neither endpoint is EXTERIOR) of ONE
EDGE, where it is a reasonable question. Asked of a whole chain it becomes far
stronger and far more fragile: one node one ulp on the wrong side condemns every
edge of the chain. That is not hypothetical. Converting (lon, lat) to xyz rounds
`cos(φ)cos(λ)` and `cos(φ)sin(λ)` independently, so two lat/lon grid cells that
share a meridian do NOT share a great circle — they get two ~1e-17 rad apart. A
cell nested inside a coarser one and sharing its left edge then has two corners a
hair OUTSIDE it and two a whole degree INSIDE, with no crossings anywhere, and
unanimity discarded the cell outright: `intersection` returned empty for a
strictly contained cell and `union` returned two polygons. (Found through
conservative regridding, where whole grid cells vanished from the weight matrix.)

Nodes are the wrong thing to poll. They are exactly where the two geometries
touch, so they are the points a sub-ulp separation is most likely to misplace,
and on a quadrilateral sharing one edge they are also evenly split — two
ambiguous, two decisive — so no tiebreak over nodes can get both the contained
case and the merely-adjacent case right. An arrangement edge, by contrast, has no
crossing in its interior by construction, so its interior lies strictly inside
one face of the other input: an unambiguous probe. Of a nested cell's four edges
only the shared one is ambiguous and three vote INTERIOR; of an adjacent cell's
four only the shared one is ambiguous and three vote EXTERIOR. Both come out
right, and by 3-to-1 rather than on a tiebreak.

It is still one verdict for the chain and still cannot split a chain across two
sides; it just declines to decide on its least representative point.
=#

# A probe strictly inside arrangement edge `(n0, n1)` — its midpoint. Arrangement
# edges carry no crossing in their interior, so this is inside one face of the
# other input, which is what makes it a better witness than either endpoint.
@inline _edge_probe_point(::Planar, arr, n0::Integer, n1::Integer) =
    let a = _node_kernel_point(arr, n0), b = _node_kernel_point(arr, n1)
        (0.5 * (GI.x(a) + GI.x(b)), 0.5 * (GI.y(a) + GI.y(b)))
    end

# On the sphere the midpoint of the CHORD normalizes to the midpoint of the arc,
# and the chord midpoint is what the labeller wants anyway — the two agree in
# direction. (An arrangement edge is never a half great circle, so `a + b` cannot
# vanish: `AntipodalEdgeSplit` is a precondition of noding.)
@inline _edge_probe_point(::Spherical, arr, n0::Integer, n1::Integer) =
    let a = _node_kernel_point(arr, n0), b = _node_kernel_point(arr, n1)
        rk_normalize_usp(UnitSphericalPoint(GI.x(a) + GI.x(b), GI.y(a) + GI.y(b),
                                            GI.z(a) + GI.z(b)))
    end
function _label_disconnected_area!(g::OverlayGraph, input, gi::Integer)
    edges = g.edges
    nnodes = length(g.node_edges)
    #-- per-node transparency cache: 0 = not computed, 1 = transparent, 2 = opaque
    trans = zeros(Int8, nnodes)
    stack = Int32[]

    #-- Phase A: flood a known location out of every NotPart edge that already
    #-- has one (pass 1 labelled it at its far node) across transparent nodes.
    for i in eachindex(edges)
        lbl = oe_label(edges, i)
        (is_not_part(lbl, gi) && !is_line_location_unknown(lbl, gi)) || continue
        push!(stack, Int32(i))
    end
    while !isempty(stack)
        e = pop!(stack)
        _node_is_transparent(g, gi, trans, he_origin(edges, e)) || continue
        loc = get_line_location(oe_label(edges, e), gi)
        f = he_onext(edges, e)
        while f != e
            lbl = oe_label(edges, f)
            if is_line_location_unknown(lbl, gi)
                set_location_all!(lbl, gi, loc)
                push!(stack, he_sym(edges, f))
            end
            f = he_onext(edges, f)
        end
    end

    #-- Phase B: one PIP verdict per remaining chain.
    visited = falses(length(edges))
    comp = Int32[]
    for i0 in eachindex(edges)
        (visited[i0] || !is_line_location_unknown(oe_label(edges, i0), gi)) && continue
        empty!(comp)
        visited[i0] = true
        push!(stack, Int32(i0))
        while !isempty(stack)
            e = pop!(stack)
            push!(comp, e)
            s = he_sym(edges, e)
            if !visited[s]
                visited[s] = true
                push!(stack, s)
            end
            nid = he_origin(edges, e)
            _node_is_transparent(g, gi, trans, nid) || continue
            f = he_onext(edges, e)
            while f != e
                if !visited[f] && is_line_location_unknown(oe_label(edges, f), gi)
                    visited[f] = true
                    push!(stack, f)
                end
                f = he_onext(edges, f)
            end
        end
        #-- one verdict for the chain, from the majority of its EDGE INTERIORS
        #-- (see the discussion above)
        nout = 0; nin = 0
        for e in comp
            p = _edge_probe_point(input.m, g.arr, he_origin(edges, e), he_dest(edges, e))
            if _input_locate_in_area(input, gi, p) == LOC_EXTERIOR
                nout += 1
            else
                nin += 1
            end
        end
        loc = nout > nin ? LOC_EXTERIOR : LOC_INTERIOR
        for e in comp
            set_location_all!(oe_label(edges, e), gi, loc)
        end
    end
    return nothing
end

# Whether input `gi`'s area location is the same on every edge incident on node
# `nid`: true iff the node carries no `gi`-boundary edge and its star is whole.
# A truncated node's star was thinned by clip pruning (split.jl), so a
# `gi`-boundary edge may simply be missing from it — the same reason pass 1 skips
# those nodes. Memoized in `trans` (0 unset / 1 transparent / 2 opaque).
function _node_is_transparent(g::OverlayGraph, gi::Integer, trans::Vector{Int8}, nid::Integer)
    t = @inbounds trans[nid]
    t != 0 && return t == 1
    ne = @inbounds g.node_edges[nid]
    truncated = g.arr.truncated
    ok = ne != 0 && !(!isempty(truncated) && (@inbounds truncated[nid])) &&
         _find_propagation_start_edge(g.edges, ne, gi) == 0
    @inbounds trans[nid] = ok ? Int8(1) : Int8(2)
    return ok
end

# ## Result-area marking (ports of `markResultAreaEdges` / `unmarkDuplicate…`)

function _mark_result_area_edges!(g::OverlayGraph, op)
    edges = g.edges
    for i in eachindex(edges)
        _mark_in_result_area!(edges, i, op)
    end
    return nothing
end

# Port of `markInResultArea`: mark an edge whose right-side (boundary) or line
# location makes it part of the result-area boundary under `op` (a set-op code
# or any location predicate — see `_is_result_of_op`).
@inline function _mark_in_result_area!(edges, i::Integer, op)
    label = oe_label(edges, i)
    if is_boundary_either(label) && _is_result_of_op(op,
            oe_get_location_boundary_or_line(edges, i, 0, POS_RIGHT),
            oe_get_location_boundary_or_line(edges, i, 1, POS_RIGHT))
        oe_mark_in_result_area!(edges, i)
    end
    return nothing
end

# Port of `unmarkDuplicateEdgesFromResultArea`: an edge whose sym is also in the
# result area cancels (merges edge-adjacent result areas per polygon validity).
function _unmark_duplicate_edges_from_result_area!(g::OverlayGraph)
    edges = g.edges
    for i in eachindex(edges)
        oe_in_result_area_both(edges, i) && oe_unmark_from_result_area_both!(edges, i)
    end
    _check_result_area_balance(g)
    return nothing
end

# ## F0: the result-area degree-balance invariant

#=
Once the result area is marked, every node must have as many in-result-area
half-edges leaving it as arriving: the marked half-edges are precisely the
boundary of the result region traversed with the region on the right, and a
boundary that enters a node has to leave it again.

Everything downstream assumes this. `maximal_edge_ring.jl` links each marked
incoming half-edge to a marked outgoing one and throws `"Ring edge missing at
max-ring build"` when it runs out — which is a report of the symptom, several
steps and one data structure away from the node whose labels are inconsistent.
This check is O(#half-edges), runs on every overlay, and names the node.

Truncated nodes are excluded for the same reason pass 1 skips them: clip pruning
(split.jl) removed part of their star, so a boundary genuinely can enter one and
not leave. They lie strictly outside the clip box and contribute nothing to the
result there.
=#
function _check_result_area_balance(g::OverlayGraph)
    edges = g.edges
    truncated = g.arr.truncated
    any_truncated = !isempty(truncated)
    for nid in eachindex(g.node_edges)
        ne = @inbounds g.node_edges[nid]
        ne == 0 && continue
        (any_truncated && (@inbounds truncated[nid])) && continue
        n_out = 0
        n_in = 0
        e = ne
        while true
            oe_in_result_area(edges, e) && (n_out += 1)
            oe_in_result_area(edges, he_sym(edges, e)) && (n_in += 1)
            e = he_onext(edges, e)
            e == ne && break
        end
        n_out == n_in ||
            throw(_OverlayTopologyError(_result_area_balance_message(g, nid, ne, n_out, n_in)))
    end
    return nothing
end

# The diagnostic body: which node, where it emits to, and its whole star with the
# per-half-edge result-area marks, so the inconsistent pair is readable directly
# off the message.
function _result_area_balance_message(g::OverlayGraph, nid::Integer, ne::Integer,
        n_out::Integer, n_in::Integer)
    edges = g.edges
    io = IOBuffer()
    print(io, "result-area degree imbalance at node ", nid,
          " ", node_point(g.arr, nid),
          ": out=", n_out, " in=", n_in, "; star (edge => dest, out/in):")
    e = Int32(ne)
    while true
        print(io, " ", e, "=>", he_dest(edges, e), " ",
              oe_in_result_area(edges, e) ? "1" : "0", "/",
              oe_in_result_area(edges, he_sym(edges, e)) ? "1" : "0")
        e = he_onext(edges, e)
        e == ne && break
    end
    return String(take!(io))
end
