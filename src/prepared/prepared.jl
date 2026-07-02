#=
# Prepared geometry: builder and v1 preparations

`prepare(geom; preps, manifold)` rebuilds the GeoInterface tree once, bottom-up, wrapping
each node in `GeometryOpsCore.Prepared` with a manifold-aware cached extent
(`rk_interaction_bounds`), and attaching each spec's built preparation at the *highest*
tree level whose trait matches `appliesto(spec)` (topmost-wins). Modeled on relateng's
`_relate_cache_extents` (`relate_geometry.jl`), which it subsumes for `Prepared` inputs.

Design: `docs/plans/2026-07-01-prepared-geometry-design.md`.
=#

export RingEdgeIndex, ChildTree, PointInArea
export SpatialEdgeIndex, SpatialIndex, PointInAreaIndex

"""
    prepare(geom; preps::Tuple = (), manifold::Manifold = Planar())

Wrap `geom` (and, recursively, its parts) in [`Prepared`](@ref), eagerly building the
preparations requested in `preps` (a tuple of `PreparationSpec`s such as
[`RingEdgeIndex`](@ref), [`ChildTree`](@ref), [`PointInArea`](@ref)).

Every wrapped node caches its manifold-aware extent. A spec is consumed at the highest
matching level of the tree. Preparations do not survive coordinate transformations —
re-`prepare` after `apply`-style rebuilds.
"""
prepare(geom; preps::Tuple = (), manifold::Manifold = Planar()) =
    _prepare(manifold, preps, GI.trait(geom), geom)

#-- split the spec tuple: consumed at this node vs offered to children (topmost-wins)
_specs_here(specs::Tuple, trait) = filter(s -> trait isa appliesto(s), specs)
_specs_below(specs::Tuple, trait) = filter(s -> !(trait isa appliesto(s)), specs)

#-- points and multipoints pass through: their extent is themselves
_prepare(::Manifold, specs, ::Union{GI.AbstractPointTrait, GI.AbstractMultiPointTrait}, geom) = geom

#-- curve leaves: the only level where coordinates are read
function _prepare(m::Manifold, specs, trait::GI.AbstractCurveTrait, line)
    GI.isempty(line) && return line
    return _wrap(m, _specs_here(specs, trait), line, rk_interaction_bounds(m, line))
end

function _prepare(m::Manifold, specs, trait::GI.AbstractPolygonTrait, poly)
    GI.isempty(poly) && return poly
    below = _specs_below(specs, trait)
    rings = [_prepare(m, below, GI.trait(r), r) for r in GI.getring(poly)]
    rebuilt = GI.geointerface_geomtype(trait)(rings; crs = GI.crs(poly))
    ext = _union_prepared_extents(m, rings)
    ext === nothing && return rebuilt
    return _wrap(m, _specs_here(specs, trait), rebuilt, ext)
end

#-- collections (covers Multi* types too)
function _prepare(m::Manifold, specs, trait::GI.AbstractGeometryCollectionTrait, geom)
    GI.isempty(geom) && return geom
    below = _specs_below(specs, trait)
    children = [_prepare(m, below, GI.trait(g), g) for g in GI.getgeom(geom)]
    rebuilt = GI.geointerface_geomtype(trait)(children; crs = GI.crs(geom))
    ext = _union_prepared_extents(m, children)
    ext === nothing && return rebuilt
    return _wrap(m, _specs_here(specs, trait), rebuilt, ext)
end

#-- unknown traits pass through untouched
_prepare(::Manifold, specs, ::GI.AbstractTrait, geom) = geom

function _wrap(m::Manifold, specs::Tuple, geom, ext)
    return Prepared(geom, m, _build_preps(specs, m, geom), ext)
end

#-- build each consumed spec; a spec may decline (empty geometry) by returning `nothing`
_build_preps(::Tuple{}, m, geom) = ()
function _build_preps(specs::Tuple, m, geom)
    p = buildprep(first(specs), m, geom)
    rest = _build_preps(Base.tail(specs), m, geom)
    return p === nothing ? rest : (p, rest...)
end

function _union_prepared_extents(m::Manifold, children)
    ext = nothing
    for c in children
        ce = if c isa Prepared
            c.extent
        elseif GI.isempty(c)
            nothing
        else
            rk_interaction_bounds(m, c)
        end
        ce === nothing && continue
        ext = ext === nothing ? ce : Extents.union(ext, ce)
    end
    return ext
end

# ## `RingEdgeIndex` — per-ring segment extent tree

"""
    SpatialEdgeIndex(tree)

Built preparation: an extent tree (anything implementing SpatialTreeInterface) over the
segments of the wrapped curve, in traversal order. Capability: `SpatialEdgeIndexLike`.
"""
struct SpatialEdgeIndex{T} <: AbstractPreparation
    tree::T
end
preptrait(::SpatialEdgeIndex) = SpatialEdgeIndexLike()

"""
    RingEdgeIndex(; nodecapacity = 16)

Spec: build a `NaturalIndex` over each linear ring's segment extents (manifold-aware:
3D arc extents on `Spherical`). Replaces the `NaturallyIndexedRing` experiment.
"""
struct RingEdgeIndex <: PreparationSpec
    nodecapacity::Int
end
RingEdgeIndex(; nodecapacity::Integer = 16) = RingEdgeIndex(Int(nodecapacity))
appliesto(::RingEdgeIndex) = GI.LinearRingTrait

function buildprep(spec::RingEdgeIndex, m::Manifold, ring)
    exts = _curve_segment_extents(m, ring)
    isempty(exts) && return nothing
    return SpatialEdgeIndex(NaturalIndex(exts; nodecapacity = spec.nodecapacity))
end

#-- per-segment manifold-aware extents of a curve, in point order (planar 2D boxes,
#-- 3D great-circle arc extents on the sphere). Mirrors the kernel's segment-extent path
#-- (`_segment_extent_table`) but per-ring and without the relateng owners table.
function _curve_segment_extents(m::Manifold, curve)
    exts = _segment_extent_type(m)[]
    n = GI.npoint(curve)
    n < 2 && return exts
    sizehint!(exts, n - 1)
    prev = _to_kernel_point(m, GI.getpoint(curve, 1))
    for i in 2:n
        cur = _to_kernel_point(m, GI.getpoint(curve, i))
        push!(exts, _segment_extent(m, prev, cur))
        prev = cur
    end
    return exts
end

# ## `ChildTree` — extent tree over child elements

"""
    SpatialIndex(tree)

Built preparation: an extent tree over the wrapped geometry's *child elements* (rings of a
polygon, polygons of a multipolygon, members of a collection). Empty children are skipped,
so leaf `i` is the `i`-th *non-empty* child in `GI.getgeom` order — not necessarily original
child `i`. Callers whose geometries may contain empty children must account for that shift.
Capability: `SpatialIndexLike`.
"""
struct SpatialIndex{T} <: AbstractPreparation
    tree::T
end
preptrait(::SpatialIndex) = SpatialIndexLike()

"""
    ChildTree(; nodecapacity = 32)

Spec: build a `NaturalIndex` over child-element extents. Attaches (topmost-wins) to
polygons, multi-geometries, and geometry collections. Empty children are skipped, so leaf
`i` in the built [`SpatialIndex`](@ref) is the `i`-th *non-empty* child in `GI.getgeom`
order — callers with possibly-empty children must account for the shift.
"""
struct ChildTree <: PreparationSpec
    nodecapacity::Int
end
ChildTree(; nodecapacity::Integer = 32) = ChildTree(Int(nodecapacity))
appliesto(::ChildTree) = Union{GI.PolygonTrait, GI.AbstractMultiPolygonTrait,
    GI.AbstractMultiCurveTrait, GI.GeometryCollectionTrait}

#-- `_segment_extent_type(m)` doubles as "the manifold's extent type" — segment and
#-- interaction bounds share it. Empty children are skipped, so if empties make child
#-- indices ambiguous for a caller, that caller filters; v1 keeps build simple.
function buildprep(spec::ChildTree, m::Manifold, geom)
    exts = _segment_extent_type(m)[]
    for c in GI.getgeom(geom)
        GI.isempty(c) && continue
        push!(exts, c isa Prepared ? c.extent : rk_interaction_bounds(m, c))
    end
    isempty(exts) && return nothing
    return SpatialIndex(NaturalIndex(exts; nodecapacity = spec.nodecapacity))
end

# ## `PointInArea` — planar point-in-polygonal-area locator

"""
    PointInAreaIndex(locator)

Built preparation: an `IndexedPointInAreaLocator` (sorted leaves) over the wrapped
polygonal geometry. Capability: `PointInAreaLike`. `Planar` only.
"""
struct PointInAreaIndex{L} <: AbstractPreparation
    locator::L
end
preptrait(::PointInAreaIndex) = PointInAreaLike()

"""
    PointInArea(; exact = True())

Spec: build an indexed point-in-area locator for a polygon or multipolygon. `Planar` only —
the spherical point-in-area index is deliberately unimplemented (relateng falls back to an
O(n) `rk_point_in_ring` walk on `Spherical`).

Reuse note: relateng's point locator reuses this prepared locator only when the algorithm's
`exact` type matches the prep's (both default `True()`); otherwise it builds its own.
"""
struct PointInArea{E} <: PreparationSpec
    exact::E
end
PointInArea(; exact = True()) = PointInArea(exact)
appliesto(::PointInArea) = Union{GI.PolygonTrait, GI.MultiPolygonTrait}

#-- Planar only: the rebuilt polygon/multipolygon whose rings are themselves `Prepared` is
#-- fed straight to the locator (GI forwarding makes the prepared children transparent).
buildprep(spec::PointInArea, m::Planar, geom) =
    PointInAreaIndex(IndexedPointInAreaLocator(m, geom; exact = spec.exact, sort_leaves = true))
buildprep(::PointInArea, m::Manifold, _) =
    throw(ArgumentError("`PointInArea` requires the `Planar` manifold (as of now); got \
        `$(typeof(m))`. A point-in-area index is implemented for `Planar` only — for any \
        other manifold omit this spec (on `Spherical`, relateng falls back to an O(n) ring \
        walk instead)."))
