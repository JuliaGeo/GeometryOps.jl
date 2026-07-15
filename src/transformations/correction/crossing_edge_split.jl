# # Crossing Edge Split

export CrossingEdgeSplit

#=
On the sphere a ring that is simple in lon/lat can still self-intersect: two
non-adjacent edges, reinterpreted as great-circle arcs, can cross properly
(a planar needle a few meters wide is enough — the arcs bulge by more than
the needle's width; Natural Earth 110m Sudan is a real instance). Standard
planar validity tooling cannot see this class, and left in place one
crossing turns the ring into a figure-eight whose containment topology is
not the one the data meant. `prepare` on the `Spherical` manifold therefore
rejects such rings by default (see its `validate` docstring), naming this
correction as the remedy.

The repair splits each ring at its proper crossing points and reassembles
the resulting loops as separate rings — even-odd semantics, both lobes
kept, the same resolution as S2Builder's undirected
`split_crossing_edges(true)` repair (which recovers planar even-odd truth
exactly on this class). A polygon whose shell splits therefore becomes a
`MultiPolygon`.

## Example
=#
# ```@example crossingsplit
# import GeometryOps as GO, GeoInterface as GI
# # an explicit figure-eight: edges 1 and 3 cross near (0, 0)
# bowtie = GI.Polygon([GI.LinearRing([(-10., -10.), (10., 10.), (10., -10.), (-10., 10.), (-10., -10.)])])
# GO.fix(bowtie; corrections = [GO.CrossingEdgeSplit()])
# ```
#=
The result is a `MultiPolygon` of the two lobes, sharing the constructed
crossing vertex, after which `GO.prepare(GO.RelateNG(; manifold =
GO.Spherical()), …)` validates clean.

## Implementation
=#

"""
    CrossingEdgeSplit() <: GeometryCorrection

Split every polygon ring at the points where two of its non-adjacent edges
cross properly as great-circle arcs, reassembling the resulting loops as
separate rings (even-odd semantics: both lobes of a figure-eight are kept,
matching S2Builder's undirected `split_crossing_edges` repair). A polygon
whose shell splits becomes a `MultiPolygon`; when it carries holes, each
(likewise repaired) hole loop is assigned to the shell loop that contains
it. This is the remedy for the ring-crossing `ArgumentError` thrown by
`prepare` on the `Spherical` manifold.

Crossing points are constructed in `Float64` (lon/lat of the exact crossing
direction) — corrections construct geometry; they don't decide predicates —
the same standard as [`AntipodalEdgeSplit`](@ref)'s midpoint insertion.

!!! warning "Scope"
    The correction handles *isolated pairwise* crossings — the needle and
    bowtie class, where no edge participates in more than one crossing and
    no two crossings interleave around the ring. Rings with tangled
    multi-crossing topology beyond that throw an `ArgumentError` rather
    than emitting wrong geometry. Vertex touches (rings meeting at a point)
    are valid and are not split.

It can be called on any polygonal geometry as usual
(`CrossingEdgeSplit()(geom)`), or passed to `GeometryOps.fix`. Because a
split changes the geometry type (`Polygon` → `MultiPolygon`), apply it
directly to `MultiPolygon` inputs rather than through `fix`'s per-polygon
traversal.

See also [`GeometryCorrection`](@ref), [`AntipodalEdgeSplit`](@ref).
"""
struct CrossingEdgeSplit <: GeometryCorrection end

application_level(::CrossingEdgeSplit) = GI.PolygonTrait

function (::CrossingEdgeSplit)(::GI.PolygonTrait, polygon)
    shell_loops = _split_ring_at_crossings(GI.getexterior(polygon))
    hole_splits = [_split_ring_at_crossings(h) for h in GI.gethole(polygon)]
    #-- identity, no copy, when nothing crossed
    shell_loops === nothing && all(isnothing, hole_splits) && return polygon

    shells = shell_loops === nothing ?
        [_ring_lonlat_open(GI.getexterior(polygon))] : shell_loops
    holes = Vector{Vector{Tuple{Float64, Float64}}}()
    for (hole, split) in zip(GI.gethole(polygon), hole_splits)
        split === nothing ? push!(holes, _ring_lonlat_open(hole)) : append!(holes, split)
    end

    _close(loop) = GI.LinearRing([loop; [loop[1]]])
    length(shells) == 1 &&
        return GI.Polygon([_close(shells[1]), map(_close, holes)...])

    #-- the shell split: assign each hole loop to the shell loop enclosing it
    m = Spherical()
    shell_krs = [SphericalKernelRing(m, _close(s); exact = True()) for s in shells]
    shell_holes = [Vector{Vector{Tuple{Float64, Float64}}}() for _ in shells]
    for h in holes
        s = _enclosing_shell(m, shell_krs, h)
        s === nothing && throw(ArgumentError(
            "CrossingEdgeSplit: a hole loop lies in no shell loop after the " *
            "split — the ring topology is tangled beyond the isolated-crossing " *
            "class this correction repairs"))
        push!(shell_holes[s], h)
    end
    return GI.MultiPolygon([GI.Polygon([_close(s), map(_close, hs)...])
                            for (s, hs) in zip(shells, shell_holes)])
end

function (c::CrossingEdgeSplit)(::GI.MultiPolygonTrait, mp)
    polys = collect(GI.getgeom(mp))
    fixed = [c(GI.PolygonTrait(), p) for p in polys]
    all(f === p for (f, p) in zip(fixed, polys)) && return mp
    out = Any[]
    for f in fixed
        GI.trait(f) isa GI.MultiPolygonTrait ? append!(out, GI.getgeom(f)) : push!(out, f)
    end
    return GI.MultiPolygon(out)
end

# The open lon/lat vertex list of a ring (closing duplicate dropped).
function _ring_lonlat_open(ring)
    ll = [(Float64(GI.x(p)), Float64(GI.y(p))) for p in GI.getpoint(ring)]
    length(ll) > 1 && ll[end] == ll[1] && pop!(ll)
    return ll
end

#=
Split one ring at its proper great-circle crossings: `nothing` when there
are none (the caller keeps the input, no copy), otherwise the open lon/lat
loops of the split. Detection is the `prepare`-validation predicate
(`_edges_cross_properly` over non-adjacent edge pairs, arc-extent pruned);
the isolation checks reject the tangled class upfront. Splitting recurses:
cut at one crossing, re-detect in each resulting loop — with a laminar
(non-interleaved), edge-disjoint crossing family every remaining crossing
falls wholly inside one loop, and each cut consumes one crossing, so the
recursion performs exactly `k` cuts for `k` crossings. The budget guards
the pathological case where a constructed (rounded) crossing vertex creates
a crossing that was not in the input family.
=#
function _split_ring_at_crossings(ring)
    ll = _ring_lonlat_open(ring)
    kp = [_spherical_kernel_point(p) for p in ll]
    _dedup_ring_vertices!(ll, kp)
    crossings = _ring_proper_crossings(kp)
    isempty(crossings) && return nothing
    _check_isolated_crossings(crossings)
    out = Vector{Vector{Tuple{Float64, Float64}}}()
    budget = Ref(length(crossings))
    _split_loops!(out, ll, kp, budget)
    return out
end

# Consecutive duplicate kernel vertices (repeated input points) collapse in
# tandem in both coordinate lists, so edge indices stay aligned.
function _dedup_ring_vertices!(ll, kp)
    i = 1
    while i <= length(kp) && length(kp) > 1
        j = mod1(i + 1, length(kp))
        if kp[i] == kp[j]
            deleteat!(ll, j)
            deleteat!(kp, j)
        else
            i += 1
        end
    end
    return nothing
end

# All proper great-circle crossings between non-adjacent edges of the open
# vertex list `kp` (edge k = kp[k] → kp[mod1(k + 1, n)]), as sorted index
# pairs. Arc-extent pruned: a nested pair loop below the accelerator
# threshold, an `RTree` self-join (the `prepare`-validation pattern) above.
function _ring_proper_crossings(kp::Vector)
    n = length(kp)
    out = Tuple{Int, Int}[]
    n < 4 && return out
    exts = [spherical_arc_extent(kp[k], kp[mod1(k + 1, n)]) for k in 1:n]
    if n < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS
        for i in 1:n, j in (i + 1):n
            Extents.intersects(exts[i], exts[j]) || continue
            _push_ring_crossing!(out, kp, n, i, j)
        end
    else
        tree = RTree(Unsorted(), collect(1:n); extents = exts, nodecapacity = 16)
        SpatialTreeInterface.dual_depth_first_search(Extents.intersects, tree, tree) do ia, ib
            ia < ib || return nothing
            _push_ring_crossing!(out, kp, n, tree.data[ia], tree.data[ib])
            return nothing
        end
        sort!(out)
    end
    return out
end

function _push_ring_crossing!(out, kp, n, i, j)
    a0 = kp[i]; a1 = kp[mod1(i + 1, n)]
    b0 = kp[j]; b1 = kp[mod1(j + 1, n)]
    #-- shared endpoints (adjacency, vertex touches) are not crossings
    (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1) && return nothing
    _edges_cross_properly(Spherical(), a0, a1, b0, b1; exact = True()) &&
        push!(out, (i, j))
    return nothing
end

@noinline _throw_tangled_crossings(why) = throw(ArgumentError(
    "CrossingEdgeSplit: the ring's proper crossings are tangled ($why); " *
    "this correction repairs isolated pairwise crossings (the needle/bowtie " *
    "class) only, and refuses rather than emit wrong geometry"))

# Repairability: each edge in at most one crossing, and no two crossings
# interleaved around the ring (as chords of the vertex circle they must
# nest or be disjoint — a laminar family — for the cut loops to be simple).
function _check_isolated_crossings(crossings)
    for a in eachindex(crossings)
        (i1, j1) = crossings[a]
        for b in (a + 1):lastindex(crossings)
            (i2, j2) = crossings[b]
            (i1 == i2 || i1 == j2 || j1 == i2 || j1 == j2) &&
                _throw_tangled_crossings("an edge crosses two others")
            (i1 < i2 < j1) == (i1 < j2 < j1) ||
                _throw_tangled_crossings("two crossings interleave")
        end
    end
    return nothing
end

function _split_loops!(out, ll::Vector, kp::Vector, budget::Base.RefValue{Int})
    crossings = _ring_proper_crossings(kp)
    if isempty(crossings)
        push!(out, ll)
        return nothing
    end
    _check_isolated_crossings(crossings)
    budget[] -= 1
    budget[] < 0 && _throw_tangled_crossings(
        "splitting produced more crossings than the input family")
    n = length(kp)
    (i, j) = crossings[1]
    x_ll = _ring_crossing_point(kp[i], kp[mod1(i + 1, n)], kp[j], kp[mod1(j + 1, n)])
    #-- the loop vertices stay lon/lat; their kernel points are re-derived
    #-- from the STORED coordinates, so the recursion (and any later
    #-- `prepare` of the output) sees exactly the emitted geometry
    x_kp = _spherical_kernel_point(x_ll)
    ja = (i + 1):j                                   # loop A: X, v[i+1] … v[j]
    jb = j == n ? (1:i) : [(j + 1):n; 1:i]           # loop B: X, v[j+1] … v[i]
    _split_loops!(out, [[x_ll]; ll[ja]], [[x_kp]; kp[ja]], budget)
    _split_loops!(out, [[x_ll]; ll[jb]], [[x_kp]; kp[jb]], budget)
    return nothing
end

# The crossing point of two properly crossing arcs, as Float64 lon/lat: the
# exact crossing direction (`_sph_crossing_dir` over the canonical crossing
# node — the candidate strictly interior to both arcs), rounded once.
function _ring_crossing_point(a0, a1, b0, b1)
    d = _sph_crossing_dir(True(), crossing_node(a0, a1, b0, b1))
    x = normalize(SVector{3, Float64}(Float64(d[1]), Float64(d[2]), Float64(d[3])))
    return GeographicFromUnitSphere()(x)
end

# The index of the shell loop whose enclosed region contains hole loop `h`:
# decided at the first hole vertex that is not on the candidate shell's
# boundary. `nothing` if every shell rejects it.
function _enclosing_shell(m::Spherical, shell_krs, h)
    for (s, kr) in enumerate(shell_krs)
        for v in h
            loc = rk_point_in_ring(m, _spherical_kernel_point(v), kr; exact = True())
            loc == LOC_BOUNDARY && continue
            loc == LOC_INTERIOR && return s
            break   # exterior of this shell: try the next one
        end
    end
    return nothing
end
