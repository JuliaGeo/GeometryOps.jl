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
    _overlay_ng(m, op::_OverlayOpCode, a, b; exact=True(), tree_a=nothing, tree_b=nothing)

Compute the overlay of `a` and `b` under `op` on manifold `m`, returning a
GeoInterface geometry. Internal engine entry point for OverlayNG — point, line
and area inputs are supported in any A×B combination. `tree_a`/`tree_b` accept
caller-prebuilt segment indices (threaded to the noding substrate).
"""
function _overlay_ng(m::Manifold, op::_OverlayOpCode, a, b;
        exact = True(), tree_a = nothing, tree_b = nothing)
    dim_a = _overlay_dimension(a)
    dim_b = _overlay_dimension(b)

    input = _OverlayInput(m, a, b, dim_a, dim_b, exact, _ov_isempty(a), _ov_isempty(b),
                          nothing, nothing)

    #-- empty-input / disjoint-envelope short circuits (port of isEmptyResult)
    er = _empty_result_short_circuit(m, op, input)
    er === nothing || return er

    #-- point paths (port of the `getResult` dispatch): points cannot node
    #-- anything, so they never reach the arrangement
    if dim_a == 0 && dim_b == 0
        return _overlay_points(m, op, a, b)                  # isAllPoints
    elseif dim_a == 0 || dim_b == 0
        return _overlay_mixed_points(m, op, a, b, dim_a, dim_b; exact)  # hasPoints
    end

    arr = NodedArrangement(m, a, b; exact, tree_a, tree_b)
    g = OverlayGraph(m, arr; exact)

    _compute_labelling!(g, input)
    _mark_result_area_edges!(g, op)
    _unmark_duplicate_edges_from_result_area!(g)

    return _extract_result(m, op, g, input; exact)
end

# ## Result extraction (port of `extractResult`)

function _extract_result(m::Manifold, op::_OverlayOpCode, g::OverlayGraph, input; exact)
    result_area_edges = graph_result_area_edges(g)
    polys = _build_polygons(m, g, result_area_edges; exact)
    has_result_area = !isempty(polys)

    #-- non-strict semantics always allow result lines
    lines = _build_lines(m, g, input, has_result_area, op; exact)

    #-- only Intersection produces points from non-point inputs
    points = op == OVERLAY_INTERSECTION ? _build_points(g) : Tuple{Float64, Float64}[]

    if isempty(polys) && isempty(lines) && isempty(points)
        return _resolve_empty_result(m, op, input)
    end
    return _create_result_geometry(polys, lines, points)
end

# Port of `OverlayUtil.createResultGeometry` + `GeometryFactory.buildGeometry`:
# the most specific geometry over the A, L, P components.
function _create_result_geometry(polys, lines, points)
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
    #-- mixed dimensions: a geometry collection in A, L, P order
    comps = Any[]
    append!(comps, polys)
    append!(comps, lines)
    for p in points
        push!(comps, GI.Point(p))
    end
    return GI.GeometryCollection(comps)
end

# ## Empty / full-sphere handling

# Port of `OverlayUtil.isEmptyResult` (the input-driven short circuit). Returns
# an empty result geometry, or `nothing` if the pipeline must run.
function _empty_result_short_circuit(m::Manifold, op::_OverlayOpCode, input::_OverlayInput)
    if op == OVERLAY_INTERSECTION
        (input.empty_a || input.empty_b) && return _empty_result(op, input)
        #-- disjoint-envelope reject (planar only; the spherical box is unreliable)
        m isa Planar && _env_disjoint(input.a, input.b) && return _empty_result(op, input)
    elseif op == OVERLAY_DIFFERENCE
        input.empty_a && return _empty_result(op, input)
    else # UNION / SYMDIFFERENCE
        input.empty_a && input.empty_b && return _empty_result(op, input)
    end
    return nothing
end

@inline function _env_disjoint(a, b)
    ea = GI.extent(a); eb = GI.extent(b)
    (ea === nothing || eb === nothing) && return false
    return !Extents.intersects(ea, eb)
end

# Port of `OverlayUtil.resultDimension`.
function _result_dimension(op::_OverlayOpCode, d0::Integer, d1::Integer)
    op == OVERLAY_INTERSECTION && return min(d0, d1)
    op == OVERLAY_UNION && return max(d0, d1)
    op == OVERLAY_DIFFERENCE && return d0
    return max(d0, d1) # SYMDIFFERENCE
end

# Resolve a pipeline that produced no components. On the plane an empty result is
# always the empty geometry. On the sphere (design §3 amendment 6) a boundaryless
# area result is ambiguous between empty and the whole sphere; disambiguate by
# locating one input vertex under the op semantics, and reject a full-sphere
# result — it has no representation here (see `_FULL_SPHERE_MSG`).
function _resolve_empty_result(m::Manifold, op::_OverlayOpCode, input::_OverlayInput)
    if m isa Spherical &&
       _result_dimension(op, input.dim_a, input.dim_b) == 2 &&
       _covers_everything(m, op, input)
        throw(ArgumentError(_FULL_SPHERE_MSG))
    end
    return _empty_result(op, input)
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
    p = _first_area_vertex(input)
    loc0 = _input_is_area(input, 0) ? _input_locate_in_area(input, 0, p) : LOC_EXTERIOR
    loc1 = _input_is_area(input, 1) ? _input_locate_in_area(input, 1, p) : LOC_EXTERIOR
    return _is_result_of_op(op, loc0, loc1)
end

function _first_area_vertex(input::_OverlayInput)
    geom = _input_is_area(input, 0) ? input.a : input.b
    p = first(GI.getpoint(geom))
    return (Float64(GI.x(p)), Float64(GI.y(p)))
end

function _empty_result(op::_OverlayOpCode, input::_OverlayInput)
    dim = _result_dimension(op, input.dim_a, input.dim_b)
    return _empty_geom(dim)
end

# Port of `OverlayUtil.createEmptyResult`: an empty geometry of the given
# dimension (2 → area, 1 → line, 0 → point). GeoInterface's auto-detecting
# wrapper constructors inspect `first(geom)`, so an empty geometry must be built
# through the raw typed (`{Z,M,T,E,C}`) constructor.
#
# Deviation from Java: JTS returns *atomic* empties (`POLYGON EMPTY`,
# `LINESTRING EMPTY`, `POINT EMPTY`); this engine returns the multi form at each
# dimension, matching the non-empty results it builds (which are Multi whenever
# the component count is not exactly one).
const _EmptyRing = Vector{Tuple{Float64, Float64}}
const _EmptyPoly = GI.Polygon{false, false, Vector{_EmptyRing}, Nothing, Nothing}
const _EmptyLine = GI.LineString{false, false, Vector{Tuple{Float64, Float64}}, Nothing, Nothing}

function _empty_geom(dim::Integer)
    if dim == 2
        return GI.MultiPolygon{false, false, Vector{_EmptyPoly}, Nothing, Nothing}(
            _EmptyPoly[], nothing, nothing)
    elseif dim == 1
        return GI.MultiLineString{false, false, Vector{_EmptyLine}, Nothing, Nothing}(
            _EmptyLine[], nothing, nothing)
    end
    #-- dim == 0 (point)
    return GI.MultiPoint{false, false, Vector{Tuple{Float64, Float64}}, Nothing, Nothing}(
        Tuple{Float64, Float64}[], nothing, nothing)
end
