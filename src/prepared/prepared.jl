#=
# Prepared geometry: builder and v1 preparations

`prepare(geom; preps, manifold)` rebuilds the GeoInterface tree once, bottom-up, wrapping
each node in `GeometryOpsCore.Prepared` with a manifold-aware cached extent
(`rk_interaction_bounds`), and attaching each spec's built preparation at the *highest*
tree level whose trait matches `appliesto(spec)` (topmost-wins). Modeled on relateng's
`_relate_cache_extents` (`relate_geometry.jl`), which it subsumes for `Prepared` inputs.

Design: `docs/plans/2026-07-01-prepared-geometry-design.md`.
=#

export RingEdgeIndex
export SpatialEdgeIndex
# NOTE: the specs `ChildTree`/`PointInArea` and their built-prep types (`SpatialIndex`,
# `PointInAreaIndex`) arrive in Tasks 6–7; their `export` lines are added alongside them.

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
