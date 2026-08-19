# NOTE: This functionality is experimental and may change at any time.

# # OverlayNG driver — the internal end-to-end overlay engine
#
# Phase 2b of the OverlayNG port (design doc §3, §4). Ties the phase-1 noding
# substrate, the phase-2a graph, and the phase-2b labeller/builders into the
# internal driver `_overlay_ng`. Ports the engine core of `OverlayNG.java`
# (`getResult` dispatch, `computeEdgeOverlay` phase order, `extractResult`
# dimensional priority) and `OverlayUtil.java` (`resultDimension`,
# `createEmptyResult`, result assembly), plus the input model of
# `InputGeometry.java`. Skipped per the design: ElevationModel,
# FastOverlayFilter, strict mode, precision/PM, RingClipper/LineLimiter.
#
# This is NOT public: no exports, no `@ref` docstrings. The public opt-in
# `OverlayNG{M}` algorithm lives in api.jl; the existing
# `intersection`/`union`/`difference` defaults are untouched.
#
# Accepted inputs: point (`Point`/`MultiPoint`), line
# (`LineString`/`LinearRing`/`MultiLineString`) and area
# (`Polygon`/`MultiPolygon`) geometries, in any A×B combination — point inputs
# are routed to `_overlay_points` / `_overlay_mixed_points` (phase 3). Geometry
# collections are rejected with a clear error.

# ## Input model (port of `InputGeometry`)

# The two overlay operands with their dimensions, lazily-built area locators, and
# empty flags. `a`/`b` are boxed (`Any`) — only the cold locator-build path reads
# them; the hot pipeline runs on the type-erased graph.
mutable struct _OverlayInput{M <: Manifold, E}
    m::M
    a::Any
    b::Any
    dim_a::Int
    dim_b::Int
    exact::E
    empty_a::Bool
    empty_b::Bool
    loc_a::Any   # Union{Nothing, IndexedPointInAreaLocator}, lazy
    loc_b::Any
end

@inline _input_dim(input::_OverlayInput, gi::Integer) = gi == 0 ? input.dim_a : input.dim_b
@inline _input_is_area(input::_OverlayInput, gi::Integer) = _input_dim(input, gi) == 2
@inline _input_is_line(input::_OverlayInput, gi::Integer) = _input_dim(input, gi) == 1
@inline _input_has_edges(input::_OverlayInput, gi::Integer) = _input_dim(input, gi) > 0

# Port of `InputGeometry.getAreaIndex`: the index of an area input, or `-1`.
@inline _input_area_index(input::_OverlayInput) =
    input.dim_a == 2 ? 0 : (input.dim_b == 2 ? 1 : -1)

# Port of `InputGeometry.locatePointInArea`: locate an emitted point against the
# ORIGINAL input area (design §3 amendment 7). Empty geometries locate EXTERIOR;
# the indexed locator is built once per side on first use.
function _input_locate_in_area(input::_OverlayInput, gi::Integer, pt)
    if gi == 0
        input.empty_a && return LOC_EXTERIOR
        input.loc_a === nothing &&
            (input.loc_a = IndexedPointInAreaLocator(input.m, input.a; exact = input.exact))
        return locate(input.loc_a, pt)
    else
        input.empty_b && return LOC_EXTERIOR
        input.loc_b === nothing &&
            (input.loc_b = IndexedPointInAreaLocator(input.m, input.b; exact = input.exact))
        return locate(input.loc_b, pt)
    end
end

# Dimension of an overlay operand: 2 (area), 1 (line), 0 (point). Geometry
# collections and other traits are unsupported here (phase 3).
function _overlay_dimension(geom)
    t = GI.trait(geom)
    if t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait
        return 2
    elseif t isa GI.LineStringTrait || t isa GI.LinearRingTrait || t isa GI.MultiLineStringTrait
        return 1
    elseif t isa GI.PointTrait || t isa GI.MultiPointTrait
        return 0
    end
    throw(ArgumentError("_overlay_ng: unsupported input geometry trait $(typeof(t))"))
end

# Whether an overlay operand is empty. `GI.isempty` falls back to `false` for
# every geometry type that does not implement it — including GeoInterface's own
# wrappers, which is exactly what this engine emits — so emptiness is also
# tested structurally. `PointTrait` has no `npoint`, and a point that does not
# report itself empty always carries a coordinate.
function _ov_isempty(geom)
    GI.isempty(geom) && return true
    GI.trait(geom) isa GI.PointTrait && return false
    return GI.npoint(geom) == 0
end

# ## The driver (port of `getResult` / `computeEdgeOverlay`)

"""
    _overlay_ng(m, op::_OverlayOpCode, a, b; exact=True(), tree_a=nothing, tree_b=nothing, target=nothing, point_type=_kernel_point_type(m))

Compute the overlay of `a` and `b` under `op` on manifold `m`, returning a
GeoInterface geometry. Internal engine entry point for OverlayNG — point, line
and area inputs are supported in any A×B combination. `tree_a`/`tree_b` accept
caller-prebuilt segment indices (threaded to the noding substrate). `target`
narrows the result to one dimension — see "Result targeting" below. `point_type`
is the output coordinate type (`OverlayNG`'s keyword of the same name).

`point_type` becomes a POSITIONAL `::Type{T}` one call in, for the reason
`NodedArrangement` does the same: a keyword-bound static parameter does not
specialize, and every result type in this file is a function of `T`.
"""
_overlay_ng(m::Manifold, op::_OverlayOpCode, a, b;
        exact = True(), tree_a = nothing, tree_b = nothing, target = nothing,
        point_type = _kernel_point_type(m)) =
    _overlay_ng(m, point_type, op, a, b, target, exact, tree_a, tree_b)

function _overlay_ng(m::Manifold, ::Type{T}, op::_OverlayOpCode, a, b,
        target, exact, tree_a, tree_b) where {T}
    tgt = _ov_target(target)
    dim_a = _overlay_dimension(a)
    dim_b = _overlay_dimension(b)

    input = _OverlayInput(m, a, b, dim_a, dim_b, exact, _ov_isempty(a), _ov_isempty(b),
                          nothing, nothing)

    #-- a target above the result dimension can match nothing whatever the inputs
    #-- contain, so it is answered without noding at all
    _target_above_result(tgt, op, dim_a, dim_b) && return _empty_result(T, op, input, tgt)

    #-- the planar input envelopes this op can use, computed once: the disjoint
    #-- short circuit and the clip pruning below are the two readers
    ea, eb = _overlay_envelopes(m, op, input)

    #-- empty-input / disjoint-envelope short circuits (port of isEmptyResult)
    _empty_result_short_circuit(m, op, input, ea, eb) && return _empty_result(T, op, input, tgt)

    #-- point paths (port of the `getResult` dispatch): points cannot node
    #-- anything, so they never reach the arrangement
    if dim_a == 0 && dim_b == 0
        return _overlay_points(m, T, op, a, b, tgt)                  # isAllPoints
    elseif dim_a == 0 || dim_b == 0
        return _overlay_mixed_points(m, T, op, a, b, dim_a, dim_b, tgt; exact)  # hasPoints
    end

    g = _overlay_marked_graph(m, T, op, a, b, input, exact, tree_a, tree_b, ea, eb)

    return _extract_result(m, op, g, input, tgt; exact)
end

# Node the inputs, build the graph, label it, and mark its result-area edges: the
# whole pipeline up to the point where a result is extracted from it. Shared with
# `intersection_area`, which extracts an area instead of a geometry.
function _overlay_marked_graph(m::Manifold, ::Type{T}, op::_OverlayOpCode, a, b, input,
        exact, tree_a, tree_b, ea, eb) where {T}
    clip_a, clip_b = _overlay_clip_envelopes(op, ea, eb)
    #-- the positional form: a `point_type` keyword would be re-boxed into the
    #-- kwargs `NamedTuple` as a bare `DataType` and `T` would stop propagating
    arr = _noded_arrangement(m, T, a, b, exact, tree_a, tree_b, clip_a, clip_b)
    g = OverlayGraph(m, arr; exact)

    _compute_labelling!(g, input)
    _mark_result_area_edges!(g, op)
    _unmark_duplicate_edges_from_result_area!(g)
    return g
end

# ## Result extraction (port of `extractResult`)

#=
The three builds are a *dependency chain*, not three independent passes, and that
is what decides which of them a `target` can skip:

  * polygons depend on nothing downstream;
  * lines depend on `has_result_area` — `isResultLine` subsumes a line that lies
    inside the result area — which is only known once the polygons are built;
  * points depend on the in-result marks of *both* earlier builds
    (`isResultPoint` rejects a node incident on any in-result edge).

So the elision is one-directional: an areal target drops both later builds, a
line target drops the point build, and a point target drops nothing. Building the
polygons anyway for a line or point target is not an oversight — `has_result_area`
is `!isempty(polys)`, and no cheaper test is equivalent (result-area edges can be
marked and still ring up to no polygon).
=#
function _extract_result(m::Manifold, op::_OverlayOpCode, g::OverlayGraph{P, T}, input,
        target = nothing; exact) where {P, T}
    result_area_edges = graph_result_area_edges(g)
    polys = _build_polygons(m, g, result_area_edges; exact)
    has_result_area = !isempty(polys)

    #-- only Intersection produces points from non-point inputs
    build_points = op == OVERLAY_INTERSECTION && _target_needs_dim(target, 0)

    #-- non-strict semantics always allow result lines. The line pass also runs
    #-- for a point target, whose lines are then discarded: `isResultPoint` reads
    #-- the in-result-LINE marks this pass writes, and without them every node on
    #-- a result line reports as an isolated intersection point.
    lines = (_target_needs_dim(target, 1) || build_points) ?
        _build_lines(m, g, input, has_result_area, op; exact) : _result_line_type(T)[]

    points = build_points ? _build_points(g) : T[]

    if _target_is_empty(target, polys, lines, points)
        return _resolve_empty_result(m, T, op, input, target)
    end
    return _target_result(target, polys, lines, points)
end

# ## Result component types
#
# The concrete wrappers every path in the engine emits, as functions of the
# output point type `T`: `_ring_to_polygon` and `_edge_line` build them from the
# graph, and `_output_component` (the mixed-point path) reconstructs the same
# ones from the input. Named because two guarantees rest on them — an empty
# result matches the type of a non-empty one at the same dimension
# (`_empty_geom`, `_target_result`), and a targeted result has one concrete
# return type.
#
# The `Is3D` flag is `_output_is3d(T)`, not `false`: a `UnitSphericalPoint`
# carries x/y/z, so a spherical result built on one IS a 3D geometry and every
# GeoInterface consumer has to be told so. Each of these is a literal `Bool`
# away from a concrete type, so inference folds them at every call site.
_result_ring_type(::Type{T}) where {T} =
    GI.LinearRing{_output_is3d(T), false, Vector{T}, Nothing, Nothing}
_result_poly_type(::Type{T}) where {T} =
    GI.Polygon{_output_is3d(T), false, Vector{_result_ring_type(T)}, Nothing, Nothing}
_result_line_type(::Type{T}) where {T} =
    GI.LineString{_output_is3d(T), false, Vector{T}, Nothing, Nothing}
_result_point_type(::Type{T}) where {T} = GI.Point{_output_is3d(T), false, T, Nothing}
_result_component_type(::Type{T}) where {T} =
    Union{_result_poly_type(T), _result_line_type(T), _result_point_type(T)}

# Port of `OverlayUtil.createResultGeometry` + `GeometryFactory.buildGeometry`:
# the most specific geometry over the A, L, P components.
function _create_result_geometry(::Type{T}, polys, lines, points) where {T}
    has_p = !isempty(polys)
    has_l = !isempty(lines)
    has_pt = !isempty(points)
    if has_p && !has_l && !has_pt
        return length(polys) == 1 ? polys[1] : GI.MultiPolygon(polys)
    elseif !has_p && has_l && !has_pt
        return length(lines) == 1 ? lines[1] : GI.MultiLineString(lines)
    elseif !has_p && !has_l && has_pt
        return length(points) == 1 ? GI.Point(points[1]) : GI.MultiPoint(points)
    end
    #=
    Mixed dimensions: a geometry collection in A, L, P order. The element type is
    the small Union of what the builders emit, NOT `Any`, and that buys two
    things. A `Vector{Any}` parent leaves the collection's `hasz`/`hasm`
    parameters uninferrable too — the wrapper detects them by inspecting the
    elements — so the whole wrapper type stays abstract even though it is always
    `{false, false}` at runtime; with the Union it infers concretely. And
    `GI.getgeom` on the result then iterates a three-way Union rather than `Any`,
    so a caller walking the components union-splits instead of dispatching
    dynamically on each one.

    The cost is that this now *rejects* a component of any other type rather than
    silently widening to `Any`. That is deliberate: these three types are the
    engine's contract (see above), and the targeted path already enforces it.
    =#
    comps = _result_component_type(T)[]
    append!(comps, polys)
    append!(comps, lines)
    for p in points
        push!(comps, GI.Point(p))
    end
    return GI.GeometryCollection(comps)
end

# ## Result targeting (the `target` keyword)

#=
`target` fixes the dimension of the result *up front*, instead of letting the
arrangement decide it. Untargeted, an overlay returns the most specific geometry
that fits, which may be a `GeometryCollection` when components of several
dimensions survive — two polygons that overlap *and* share a boundary segment
intersect in a polygon and a line. A caller who only ever wants the areal part
then has to destructure a collection whose shape depends on the data.

With a target it does not: the return type is a function of the target alone.
The accepted targets and their returns mirror the Foster–Hormann entry points —
a singular trait gives the `Vector` of atomic components, a `Multi` trait gives
the one multi-geometry:

    dim  singular target          -> returns          multi target                 -> returns
    2    `GI.PolygonTrait()`      -> Vector{<:Polygon}     `GI.MultiPolygonTrait()`    -> MultiPolygon
    1    `GI.LineStringTrait()`   -> Vector{<:LineString}  `GI.MultiLineStringTrait()` -> MultiLineString
    0    `GI.PointTrait()`        -> Vector{<:Point}       `GI.MultiPointTrait()`      -> MultiPoint

An empty result is the empty `Vector` or the empty multi-geometry of that same
concrete type — never a geometry of another dimension, which is what the
untargeted `_empty_geom` would give. That promise rests on the `_Result*`
component types above, which every path in the engine agrees on.

The consequence is that a targeted call is type STABLE: `target` is a singleton,
so `_target_result` resolves to one method, and both its branches (empty and
not) return the same concrete type. Untargeted, the return type is a Union of
all seven geometries the engine may emit, since inference cannot know which
dimensions survive.

Targeting also *removes work*, in two places: `_target_above_result` answers
without noding whenever the target is above the OGC result dimension (an areal
target on a line × area intersection is empty for every possible input), and
`_extract_result` skips the builds the target cannot want.
=#

const _OverlayTarget = Union{GI.PointTrait, GI.MultiPointTrait,
                             GI.LineStringTrait, GI.MultiLineStringTrait,
                             GI.PolygonTrait, GI.MultiPolygonTrait}

# Normalize whatever the caller passed to `nothing` or a bare trait instance.
# Trait types and `TraitTarget`s are accepted for parity with the Foster–Hormann
# entry points, which document `target::Type` and wrap in `TraitTarget`.
_ov_target(::Nothing) = nothing
_ov_target(t::_OverlayTarget) = t
_ov_target(::TraitTarget{T}) where {T} = _ov_target(T)
function _ov_target(::Type{T}) where {T <: GI.AbstractTrait}
    #-- a `TraitTarget` over a Union of traits has no single result dimension
    isconcretetype(T) || throw(ArgumentError(_bad_target_msg(T)))
    return _ov_target(T())
end
_ov_target(x) = throw(ArgumentError(_bad_target_msg(x)))

_bad_target_msg(x) =
    "OverlayNG: unsupported `target` $(x). An overlay result is targeted by " *
    "dimension, so `target` must be `nothing` (the default — the most specific " *
    "geometry over every dimension present) or exactly one of `GI.PolygonTrait()`, " *
    "`GI.MultiPolygonTrait()`, `GI.LineStringTrait()`, `GI.MultiLineStringTrait()`, " *
    "`GI.PointTrait()` or `GI.MultiPointTrait()`."

_target_dim(::Union{GI.PolygonTrait, GI.MultiPolygonTrait}) = 2
_target_dim(::Union{GI.LineStringTrait, GI.MultiLineStringTrait}) = 1
_target_dim(::Union{GI.PointTrait, GI.MultiPointTrait}) = 0

# Whether the target is above the OGC dimension of the result, i.e. provably
# unsatisfiable from the inputs' dimensions alone.
_target_above_result(::Nothing, op::_OverlayOpCode, dim_a, dim_b) = false
_target_above_result(t, op::_OverlayOpCode, dim_a, dim_b) =
    _target_dim(t) > _result_dimension(op, dim_a, dim_b)

# Whether components of dimension `d` can appear in the targeted result — the
# builder-elision test.
_target_needs_dim(::Nothing, d::Integer) = true
_target_needs_dim(t, d::Integer) = _target_dim(t) == d

_target_admits_area(::Nothing) = true
_target_admits_area(t) = _target_dim(t) == 2

_target_is_empty(::Nothing, polys, lines, points) =
    isempty(polys) && isempty(lines) && isempty(points)
_target_is_empty(::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, polys, lines, points) =
    isempty(polys)
_target_is_empty(::Union{GI.LineStringTrait, GI.MultiLineStringTrait}, polys, lines, points) =
    isempty(lines)
_target_is_empty(::Union{GI.PointTrait, GI.MultiPointTrait}, polys, lines, points) =
    isempty(points)

# The component lists are already concrete `Vector`s on every path but one: the
# mixed-point path substitutes `_NO_COMPONENTS` (`Any[]`) for the dimension its
# non-point input does not have, and that is always empty.
@inline _typed(::Type{T}, v::Vector{T}) where {T} = v
@inline _typed(::Type{T}, v) where {T} = T[c for c in v]

# Assemble the targeted result. `nothing` is the untargeted assembler and, like
# it, requires at least one component; the six trait methods are total.
#
# The output point type is read off `points`, which is the one argument every
# path supplies as a concretely-typed `Vector{T}` (the other two may be the
# shared untyped `_NO_COMPONENTS`).
_target_result(::Nothing, polys, lines, points::Vector{T}) where {T} =
    _create_result_geometry(T, polys, lines, points)

_target_result(::GI.PolygonTrait, polys, lines, points::Vector{T}) where {T} =
    _typed(_result_poly_type(T), polys)
_target_result(::GI.LineStringTrait, polys, lines, points::Vector{T}) where {T} =
    _typed(_result_line_type(T), lines)
_target_result(::GI.PointTrait, polys, lines, points::Vector{T}) where {T} =
    _result_point_type(T)[GI.Point(p) for p in points]

#-- the wrapper constructors read `first(geom)` to detect their element type, so
#-- an empty multi-geometry has to come from `_empty_geom`'s raw typed constructor
_target_result(::GI.MultiPolygonTrait, polys, lines, points::Vector{T}) where {T} =
    (v = _typed(_result_poly_type(T), polys); isempty(v) ? _empty_geom(T, 2) : GI.MultiPolygon(v))
_target_result(::GI.MultiLineStringTrait, polys, lines, points::Vector{T}) where {T} =
    (v = _typed(_result_line_type(T), lines); isempty(v) ? _empty_geom(T, 1) : GI.MultiLineString(v))
_target_result(::GI.MultiPointTrait, polys, lines, points::Vector{T}) where {T} =
    isempty(points) ? _empty_geom(T, 0) : GI.MultiPoint(points)

_empty_target_result(::Type{T}, t) where {T} =
    _target_result(t, _result_poly_type(T)[], _result_line_type(T)[], T[])

# Empty component lists, for the dimensions a path does not produce at all.
# Shared and never mutated — `_target_result` wraps but does not take ownership.
const _NO_COMPONENTS = Any[]

_empty_dim_result(::Type{T}, target, dim::Integer) where {T} =
    target === nothing ? _empty_geom(T, dim) : _empty_target_result(T, target)

# The shared result tail of the three overlay paths: assemble the targeted
# result, falling back to the empty geometry of the OGC result dimension. Only
# `_extract_result` needs more than this — it interposes the full-sphere check.
function _dimensional_result(target, dim::Integer, polys, lines, points::Vector{T}) where {T}
    _target_is_empty(target, polys, lines, points) && return _empty_dim_result(T, target, dim)
    return _target_result(target, polys, lines, points)
end

# ## Empty / full-sphere handling

# Port of `OverlayUtil.isEmptyResult` (the input-driven short circuit). Whether
# the inputs alone force an empty result, so the pipeline need not run. `ea`/`eb`
# are the input envelopes `_overlay_envelopes` produced for this op, or `nothing`.
function _empty_result_short_circuit(m::Manifold, op::_OverlayOpCode, input::_OverlayInput,
        ea, eb)
    if op == OVERLAY_INTERSECTION
        (input.empty_a || input.empty_b) && return true
        #-- disjoint-envelope reject (planar only; the spherical box is unreliable)
        _env_disjoint(ea, eb) && return true
    elseif op == OVERLAY_DIFFERENCE
        input.empty_a && return true
    else # UNION / SYMDIFFERENCE
        input.empty_a && input.empty_b && return true
    end
    return false
end

@inline _env_disjoint(ea, eb) =
    !(ea === nothing || eb === nothing) && !Extents.intersects(ea, eb)

# ## Clip pruning (the construct-free `RingClipper`)

#=
GEOS clips both inputs to the intersection envelope before noding
(`RingClipper`), which is most of why its `intersection` is several times cheaper
than its own `union`. Closing a clipped ring along the box constructs coordinates
that then feed the noder, so that form cannot be adopted here (design §0). The
construct-free equivalent is to prune rather than clip: a parent segment whose
bbox misses the box simply produces no `NodedEdge`, and the surviving chains stay
OPEN — see the long note in noding/split.jl for why the arrangement tolerates
that, and `_compute_labelling!` for the one pass that has to know.

What this function owns is the remaining premise: **every edge that can be in the
result lies inside the box it prunes against**.

  * INTERSECTION — the result is contained in `env(A) ∩ env(B)`, and so is every
    piece of either boundary that bounds it. Prune both sides to that box.
  * DIFFERENCE (A − B) — the result is contained in A. Its boundary is made of
    pieces of ∂A (so A is never pruned) and pieces of ∂B lying inside A, hence
    inside `env(A)`. Prune only B, to `env(A)`.
  * UNION / SYMDIFFERENCE — the result reaches everywhere either input does, so
    there is no box to prune against. GEOS does not clip these either.

Spherical is excluded: the lon/lat box is unreliable across the antimeridian and
the poles, which is the same reason the disjoint-envelope short circuit above is
planar-only. A non-`_OverlayOpCode` `op` (the labeller accepts any
`(loc0, loc1) -> Bool`) is excluded too — an arbitrary predicate has no result
containment to appeal to.
=#

# The input envelopes this op can make use of, or `nothing`s. Computing them here
# rather than inside each consumer keeps `GI.extent` to one traversal per side,
# and keeps it off the ops that have no use for it. An EMPTY operand is excluded
# too: it can have no envelope (LibGEOS raises on one), and every op that would
# read it short circuits on the emptiness first.
@inline function _overlay_envelopes(m::Manifold, op::_OverlayOpCode, input::_OverlayInput)
    m isa Planar || return (nothing, nothing)
    if op == OVERLAY_INTERSECTION
        (input.empty_a || input.empty_b) && return (nothing, nothing)
        return (GI.extent(input.a), GI.extent(input.b))
    elseif op == OVERLAY_DIFFERENCE
        input.empty_a && return (nothing, nothing)
        return (GI.extent(input.a), nothing)
    end
    return (nothing, nothing)
end
@inline _overlay_envelopes(m::Manifold, op, input::_OverlayInput) = (nothing, nothing)

# The per-side clip boxes for `op`, given the envelopes above.
@inline function _overlay_clip_envelopes(op::_OverlayOpCode, ea, eb)
    if op == OVERLAY_INTERSECTION
        e = _clip_box(ea, eb)
        return (e, e)
    elseif op == OVERLAY_DIFFERENCE
        return (nothing, _clip_box(ea))
    end
    return (nothing, nothing)
end
@inline _overlay_clip_envelopes(op, ea, eb) = (nothing, nothing)

# An envelope as the concrete XY box `_seg_in_clip` tests against.
@inline _clip_box(::Nothing) = nothing
@inline _clip_box(e) = Extents.Extent((X = (Float64(e.X[1]), Float64(e.X[2])),
                                       Y = (Float64(e.Y[1]), Float64(e.Y[2]))))

@inline _clip_box(::Nothing, eb) = nothing
@inline _clip_box(ea, ::Nothing) = nothing
@inline _clip_box(::Nothing, ::Nothing) = nothing
@inline _clip_box(ea, eb) = Extents.Extent((
    X = (max(Float64(ea.X[1]), Float64(eb.X[1])), min(Float64(ea.X[2]), Float64(eb.X[2]))),
    Y = (max(Float64(ea.Y[1]), Float64(eb.Y[1])), min(Float64(ea.Y[2]), Float64(eb.Y[2])))))

# Port of `OverlayUtil.resultDimension`.
function _result_dimension(op::_OverlayOpCode, d0::Integer, d1::Integer)
    op == OVERLAY_INTERSECTION && return min(d0, d1)
    op == OVERLAY_UNION && return max(d0, d1)
    op == OVERLAY_DIFFERENCE && return d0
    return max(d0, d1) # SYMDIFFERENCE
end

#=
Whether `op` can return the full sphere at all, from the op alone.

`INTERSECTION` and `DIFFERENCE` both return a SUBSET of A, and A is a
ring-bounded polygon — it has a nonempty boundary, so it is not the full sphere
and neither is any subset of it. Only the two ops that can return a superset of
an input (`UNION`, `SYMDIFFERENCE`) can cover everything.

This is not an optimization. `_covers_everything` probes by locating a single
input VERTEX, and a vertex lies on its own input's boundary by construction —
the one place where the op's value is not representative of the neighbourhood.
Two lon/lat cells sharing an edge put that vertex on BOTH boundaries, so the
probe reported their (empty) intersection as covering the sphere. Restricting
the question to the ops that can actually answer "yes" removes the false
positive at its source rather than hardening the probe.
=#
@inline _op_can_cover_everything(op::_OverlayOpCode) =
    op == OVERLAY_UNION || op == OVERLAY_SYMDIFFERENCE

# Resolve a pipeline that produced no components. On the plane an empty result is
# always the empty geometry. On the sphere (design §3 amendment 6) a boundaryless
# area result is ambiguous between empty and the whole sphere; disambiguate by
# locating one input vertex under the op semantics, and reject a full-sphere
# result — it has no representation here (see `_FULL_SPHERE_MSG`).
#
# The rejection is areal, so it is skipped for a target that excludes areas: the
# lines or points of such an overlay are perfectly representable, and the caller
# who asked only for those has no stake in the region being unnameable.
function _resolve_empty_result(m::Manifold, ::Type{T}, op::_OverlayOpCode,
        input::_OverlayInput, target = nothing) where {T}
    if m isa Spherical && _target_admits_area(target) &&
       _result_dimension(op, input.dim_a, input.dim_b) == 2 &&
       _op_can_cover_everything(op) && _covers_everything(m, op, input)
        throw(ArgumentError(_FULL_SPHERE_MSG))
    end
    return _empty_result(T, op, input, target)
end

#=
The full sphere is the one areal region on `Spherical()` that GeometryOps cannot
denote, and that is a property of the geometry model, not of this engine:

  * a polygon denotes the region bounded by its rings. Under `Spherical(;
    oriented = false)` (the default, and the ecosystem's — S2, R/s2, Python)
    that is the *enclosed* region of the ring; under `oriented = true` it is the
    region to the ring's left. Either way a nonempty boundary is required, and
    the full sphere has none.
  * a ring-free polygon is exactly how an *empty* geometry is spelled
    (`GI.isempty` is true for it), so it cannot also mean the full sphere.

S2 solves this with a sentinel (`S2Loop::kFull()`, a one-vertex loop that every
predicate special-cases). Adding such a sentinel is a change to the meaning of
`Spherical(; oriented)` across `area`, the predicates and the locators — far
outside overlay — so this engine raises instead of guessing. Reaching this error
requires an overlay whose result has *no* boundary at all and covers everything:
in practice a union/symdifference of complementary hemispheres, or a union whose
operands' boundaries cancel exactly. If you hit it, subtract the region you care
about instead (`difference(OverlayNG(Spherical()), whole, part)` is
representable whenever `whole` is).
=#
const _FULL_SPHERE_MSG =
    "OverlayNG: the overlay result covers the whole sphere. GeometryOps has no " *
    "representation for the full sphere — a polygon denotes the region bounded by " *
    "its rings, and a ring-free polygon already means the empty geometry — so the " *
    "result cannot be returned. Reformulate the operation (e.g. as a `difference` " *
    "from the covering region) or work with the complement."

# Whether a boundaryless result covers the entire manifold: since there is no
# result boundary the result is uniform, so evaluating the op at any single point
# decides it. Uses a vertex of an area input (boundary counts as interior).
function _covers_everything(m::Manifold, op::_OverlayOpCode, input::_OverlayInput)
    p = _first_area_vertex(m, input)
    loc0 = _input_is_area(input, 0) ? _input_locate_in_area(input, 0, p) : LOC_EXTERIOR
    loc1 = _input_is_area(input, 1) ? _input_locate_in_area(input, 1, p) : LOC_EXTERIOR
    return _is_result_of_op(op, loc0, loc1)
end

# The vertex goes straight to the point-in-area locators, so it has to be in the
# manifold's KERNEL coordinates, not a bare `(x, y)` pair off the input. On the
# sphere an input vertex may already be xyz (`UnitSphericalPoint`); keeping only
# its first two components drops z, and the locator then reads the pair back as
# lon/lat — a different point on a different part of the sphere. Two lat/lon grid
# cells meeting at the equator located that ghost point as BOUNDARY of both and
# reported their (empty) intersection as covering the whole sphere.
function _first_area_vertex(m::Manifold, input::_OverlayInput)
    geom = _input_is_area(input, 0) ? input.a : input.b
    return _to_kernel_point(m, first(GI.getpoint(geom)))
end

_empty_result(::Type{T}, op::_OverlayOpCode, input::_OverlayInput, target = nothing) where {T} =
    _empty_dim_result(T, target, _result_dimension(op, input.dim_a, input.dim_b))

# Port of `OverlayUtil.createEmptyResult`: an empty geometry of the given
# dimension (2 → area, 1 → line, 0 → point). GeoInterface's auto-detecting
# wrapper constructors inspect `first(geom)`, so an empty geometry must be built
# through the raw typed (`{Z,M,T,E,C}`) constructor.
#
# Deviation from Java: JTS returns *atomic* empties (`POLYGON EMPTY`,
# `LINESTRING EMPTY`, `POINT EMPTY`); this engine returns the multi form at each
# dimension, matching the non-empty results it builds (which are Multi whenever
# the component count is not exactly one).
#
# The element types are `_Result*`, so an empty multi-geometry has the same
# concrete type as the non-empty one of that dimension.
function _empty_geom(::Type{T}, dim::Integer) where {T}
    Z = _output_is3d(T)
    if dim == 2
        return GI.MultiPolygon{Z, false, Vector{_result_poly_type(T)}, Nothing, Nothing}(
            _result_poly_type(T)[], nothing, nothing)
    elseif dim == 1
        return GI.MultiLineString{Z, false, Vector{_result_line_type(T)}, Nothing, Nothing}(
            _result_line_type(T)[], nothing, nothing)
    end
    #-- dim == 0 (point)
    return GI.MultiPoint{Z, false, Vector{T}, Nothing, Nothing}(T[], nothing, nothing)
end
