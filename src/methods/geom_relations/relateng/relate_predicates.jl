# # Named DE-9IM relate predicates
#=
Port of JTS `RelatePredicate.java` (the `IMPredicate` kinds — the
`intersects`/`disjoint` `BasicPredicate` kinds live in
`topology_predicate.jl`), `IMPatternMatcher.java`,
`IntersectionMatrixPattern.java` and `RelateMatrixPredicate.java`.

Each JTS anonymous inner class becomes a kind singleton parameterizing
`IMPredicate{K}`, with its `requireX` overrides, `init` overrides,
`isDetermined` and `valueIM` ported method-for-method, in the same
order as `RelatePredicate.java`. The `value_im` matrix queries
(`is_contains` etc.) are the ports of `geom/IntersectionMatrix.java`
named-relationship methods in `de9im.jl`.
=#

# Which input geometry a source flag refers to (JTS `RelateGeometry.GEOM_A/GEOM_B`).
const GEOM_A = true
const GEOM_B = false

# Envelope helpers for `init_bounds!` (JTS `Envelope.isNull` / `Envelope.equals`).
# A null (empty-geometry) extent is represented as `nothing` or as an
# extent with an inverted X interval, matching the JTS null-envelope convention.
ext_isnull(::Nothing) = true
ext_isnull(ext) = ext.X[2] < ext.X[1]
function ext_equals(extA, extB)
    (ext_isnull(extA) || ext_isnull(extB)) && return ext_isnull(extA) && ext_isnull(extB)
    return extA.X == extB.X && extA.Y == extB.Y
end

#=
## `contains` (RelatePredicate.java `contains()`)
=#
struct ContainsPred end
pred_contains() = IMPredicate(ContainsPred())
predicate_name(::IMPredicate{ContainsPred}) = "contains"
require_covers(::Type{IMPredicate{ContainsPred}}, is_source_a::Bool) = is_source_a == GEOM_A
# only need to check B against Exterior of A
require_exterior_check(::Type{IMPredicate{ContainsPred}}, is_source_a::Bool) = is_source_a == GEOM_B
init_dims_kind!(p::IMPredicate{ContainsPred}) =
    require!(p, is_dims_compatible_with_covers(p.dimA, p.dimB))
init_bounds!(p::IMPredicate{ContainsPred}, extA, extB) = require_covers!(p, extA, extB)
is_determined(p::IMPredicate{ContainsPred}) = intersects_exterior_of(p, GEOM_A)
value_im(p::IMPredicate{ContainsPred}) = is_contains(p.im)

#=
## `within` (RelatePredicate.java `within()`)
=#
struct WithinPred end
pred_within() = IMPredicate(WithinPred())
predicate_name(::IMPredicate{WithinPred}) = "within"
require_covers(::Type{IMPredicate{WithinPred}}, is_source_a::Bool) = is_source_a == GEOM_B
# only need to check A against Exterior of B
require_exterior_check(::Type{IMPredicate{WithinPred}}, is_source_a::Bool) = is_source_a == GEOM_A
init_dims_kind!(p::IMPredicate{WithinPred}) =
    require!(p, is_dims_compatible_with_covers(p.dimB, p.dimA))
init_bounds!(p::IMPredicate{WithinPred}, extA, extB) = require_covers!(p, extB, extA)
is_determined(p::IMPredicate{WithinPred}) = intersects_exterior_of(p, GEOM_B)
value_im(p::IMPredicate{WithinPred}) = is_within(p.im)

#=
## `covers` (RelatePredicate.java `covers()`)
=#
struct CoversPred end
pred_covers() = IMPredicate(CoversPred())
predicate_name(::IMPredicate{CoversPred}) = "covers"
require_covers(::Type{IMPredicate{CoversPred}}, is_source_a::Bool) = is_source_a == GEOM_A
# only need to check B against Exterior of A
require_exterior_check(::Type{IMPredicate{CoversPred}}, is_source_a::Bool) = is_source_a == GEOM_B
init_dims_kind!(p::IMPredicate{CoversPred}) =
    require!(p, is_dims_compatible_with_covers(p.dimA, p.dimB))
init_bounds!(p::IMPredicate{CoversPred}, extA, extB) = require_covers!(p, extA, extB)
is_determined(p::IMPredicate{CoversPred}) = intersects_exterior_of(p, GEOM_A)
value_im(p::IMPredicate{CoversPred}) = is_covers(p.im)

#=
## `coveredBy` (RelatePredicate.java `coveredBy()`)
=#
struct CoveredByPred end
pred_coveredby() = IMPredicate(CoveredByPred())
predicate_name(::IMPredicate{CoveredByPred}) = "coveredBy"
require_covers(::Type{IMPredicate{CoveredByPred}}, is_source_a::Bool) = is_source_a == GEOM_B
# only need to check A against Exterior of B
require_exterior_check(::Type{IMPredicate{CoveredByPred}}, is_source_a::Bool) = is_source_a == GEOM_A
init_dims_kind!(p::IMPredicate{CoveredByPred}) =
    require!(p, is_dims_compatible_with_covers(p.dimB, p.dimA))
init_bounds!(p::IMPredicate{CoveredByPred}, extA, extB) = require_covers!(p, extB, extA)
is_determined(p::IMPredicate{CoveredByPred}) = intersects_exterior_of(p, GEOM_B)
value_im(p::IMPredicate{CoveredByPred}) = is_coveredby(p.im)

#=
## `crosses` (RelatePredicate.java `crosses()`)
=#
struct CrossesPred end
pred_crosses() = IMPredicate(CrossesPred())
predicate_name(::IMPredicate{CrossesPred}) = "crosses"
function init_dims_kind!(p::IMPredicate{CrossesPred})
    is_both_points_or_areas = (p.dimA == DIM_P && p.dimB == DIM_P) ||
        (p.dimA == DIM_A && p.dimB == DIM_A)
    require!(p, !is_both_points_or_areas)
end
function is_determined(p::IMPredicate{CrossesPred})
    if p.dimA == DIM_L && p.dimB == DIM_L
        # L/L interaction can only be dim = P
        get_dimension(p, LOC_INTERIOR, LOC_INTERIOR) > DIM_P && return true
    elseif p.dimA < p.dimB
        if is_intersects_entry(p, LOC_INTERIOR, LOC_INTERIOR) &&
           is_intersects_entry(p, LOC_INTERIOR, LOC_EXTERIOR)
            return true
        end
    elseif p.dimA > p.dimB
        if is_intersects_entry(p, LOC_INTERIOR, LOC_INTERIOR) &&
           is_intersects_entry(p, LOC_EXTERIOR, LOC_INTERIOR)
            return true
        end
    end
    return false
end
value_im(p::IMPredicate{CrossesPred}) = is_crosses(p.im, p.dimA, p.dimB)

#=
## `equalsTopo` (RelatePredicate.java `equalsTopo()`)
=#
struct EqualsTopoPred end
pred_equalstopo() = IMPredicate(EqualsTopoPred())
predicate_name(::IMPredicate{EqualsTopoPred}) = "equals"
# don't require equal dims, because EMPTY = EMPTY for all dims
init_dims_kind!(p::IMPredicate{EqualsTopoPred}) = nothing
# allow EMPTY = EMPTY
require_interaction(::Type{IMPredicate{EqualsTopoPred}}) = false
function init_bounds!(p::IMPredicate{EqualsTopoPred}, extA, extB)
    # handle EMPTY = EMPTY cases
    set_value_if!(p, true, ext_isnull(extA) && ext_isnull(extB))
    require!(p, ext_equals(extA, extB))
    return nothing
end
function is_determined(p::IMPredicate{EqualsTopoPred})
    is_either_exterior_intersects =
        is_intersects_entry(p, LOC_INTERIOR, LOC_EXTERIOR) ||
        is_intersects_entry(p, LOC_BOUNDARY, LOC_EXTERIOR) ||
        is_intersects_entry(p, LOC_EXTERIOR, LOC_INTERIOR) ||
        is_intersects_entry(p, LOC_EXTERIOR, LOC_BOUNDARY)
    return is_either_exterior_intersects
end
value_im(p::IMPredicate{EqualsTopoPred}) = is_equals(p.im, p.dimA, p.dimB)

#=
## `overlaps` (RelatePredicate.java `overlaps()`)
=#
struct OverlapsPred end
pred_overlaps() = IMPredicate(OverlapsPred())
predicate_name(::IMPredicate{OverlapsPred}) = "overlaps"
init_dims_kind!(p::IMPredicate{OverlapsPred}) = require!(p, p.dimA == p.dimB)
function is_determined(p::IMPredicate{OverlapsPred})
    if p.dimA == DIM_A || p.dimA == DIM_P
        if is_intersects_entry(p, LOC_INTERIOR, LOC_INTERIOR) &&
           is_intersects_entry(p, LOC_INTERIOR, LOC_EXTERIOR) &&
           is_intersects_entry(p, LOC_EXTERIOR, LOC_INTERIOR)
            return true
        end
    end
    if p.dimA == DIM_L
        if is_dimension_entry(p, LOC_INTERIOR, LOC_INTERIOR, DIM_L) &&
           is_intersects_entry(p, LOC_INTERIOR, LOC_EXTERIOR) &&
           is_intersects_entry(p, LOC_EXTERIOR, LOC_INTERIOR)
            return true
        end
    end
    return false
end
value_im(p::IMPredicate{OverlapsPred}) = is_overlaps(p.im, p.dimA, p.dimB)

#=
## `touches` (RelatePredicate.java `touches()`)
=#
struct TouchesPred end
pred_touches() = IMPredicate(TouchesPred())
predicate_name(::IMPredicate{TouchesPred}) = "touches"
function init_dims_kind!(p::IMPredicate{TouchesPred})
    # Points have only interiors, so cannot touch
    is_both_points = p.dimA == DIM_P && p.dimB == DIM_P
    require!(p, !is_both_points)
end
function is_determined(p::IMPredicate{TouchesPred})
    # for touches interiors cannot intersect
    is_interiors_intersects = is_intersects_entry(p, LOC_INTERIOR, LOC_INTERIOR)
    return is_interiors_intersects
end
value_im(p::IMPredicate{TouchesPred}) = is_touches(p.im, p.dimA, p.dimB)

#=
## `IMPatternMatcher` (IMPatternMatcher.java)

A predicate that matches a DE-9IM pattern. Unlike the named kinds above
this is a standalone mutable struct (not an `IMPredicate` kind), because
its `require_interaction` flag depends on runtime data (the pattern
matrix), via the instance-level requirement-flag methods; the small
`IMPredicate` state-machine methods are mirrored below (Java gets them
by inheritance).
=#
mutable struct IMPatternMatcher <: TopologyPredicate
    const im_pattern::String
    const pattern_matrix::DE9IM
    dimA::Int8
    dimB::Int8
    im::DE9IM
    value::Int8
end
IMPatternMatcher(im_pattern::AbstractString) =
    IMPatternMatcher(String(im_pattern), DE9IM(im_pattern), DIM_UNKNOWN, DIM_UNKNOWN,
        # E/E is always dim = 2 (IMPredicate constructor)
        with_entry(DE9IM(), LOC_EXTERIOR, LOC_EXTERIOR, DIM_A), TRI_UNKNOWN)

predicate_name(::IMPatternMatcher) = "IMPattern"

# RelatePredicate.java `matches(String)` factory.
pred_matches(im_pattern::AbstractString) = IMPatternMatcher(im_pattern)

function init_dims!(p::IMPatternMatcher, dimA::Integer, dimB::Integer)
    p.dimA = dimA
    p.dimB = dimB
    return nothing
end

# if pattern specifies any non-E/non-E interaction, envelopes must not be disjoint
# (the Java method also starts with `super.init(dimA, dimB)`, which only
# re-assigns the already-set dims — a no-op not reproduced here)
function init_bounds!(p::IMPatternMatcher, extA, extB)
    requires_interaction = im_requires_interaction(p.pattern_matrix)
    is_disjoint = !ext_intersects(extA, extB)
    set_value_if!(p, false, requires_interaction && is_disjoint)
    return nothing
end

require_interaction(p::IMPatternMatcher) = im_requires_interaction(p.pattern_matrix)

# IMPatternMatcher.java static `requireInteraction(IntersectionMatrix)`
function im_requires_interaction(im::DE9IM)
    requires_interaction =
        _is_interaction(im[LOC_INTERIOR, LOC_INTERIOR]) ||
        _is_interaction(im[LOC_INTERIOR, LOC_BOUNDARY]) ||
        _is_interaction(im[LOC_BOUNDARY, LOC_INTERIOR]) ||
        _is_interaction(im[LOC_BOUNDARY, LOC_BOUNDARY])
    return requires_interaction
end
_is_interaction(im_dim::Integer) = im_dim == DIM_TRUE || im_dim >= DIM_P

# Mirrors of the inherited IMPredicate state-machine methods.
function update_dim!(p::IMPatternMatcher, locA, locB, dim)
    # only record an increased dimension value
    if dim > p.im[locA, locB]
        p.im = with_entry(p.im, locA, locB, dim)
        # set value if predicate value can be known
        if is_determined(p)
            set_value!(p, value_im(p))
        end
    end
    return nothing
end
get_dimension(p::IMPatternMatcher, locA, locB) = p.im[locA, locB]
finish!(p::IMPatternMatcher) = set_value!(p, value_im(p))

function is_determined(p::IMPatternMatcher)
    #=
    Matrix entries only increase in dimension as topology is computed.
    The predicate can be short-circuited (as false) if
    any computed entry is greater than the mask value.
    =#
    for i in 0:2, j in 0:2
        pattern_entry = p.pattern_matrix[i, j]
        pattern_entry == DIM_DONTCARE && continue
        matrix_val = get_dimension(p, i, j)
        if pattern_entry == DIM_TRUE
            # mask entry TRUE requires a known matrix entry
            matrix_val < 0 && return false
        elseif matrix_val > pattern_entry
            # result is known (false) if matrix entry has exceeded mask
            return true
        end
    end
    return false
end

value_im(p::IMPatternMatcher) = im_matches(p.im, p.im_pattern)

Base.show(io::IO, p::IMPatternMatcher) =
    print(io, predicate_name(p), "(", p.im_pattern, ")")

#=
## DE-9IM matrix pattern constants (IntersectionMatrixPattern.java)
=#

# Detects whether two polygonal geometries are adjacent along an edge,
# but do not overlap.
const IM_PATTERN_ADJACENT = "F***1****"
# Detects a geometry which properly contains another geometry (i.e. which
# lies entirely in the interior of the first geometry).
const IM_PATTERN_CONTAINS_PROPERLY = "T**FF*FF*"
# Detects if two geometries intersect in their interiors.
const IM_PATTERN_INTERIOR_INTERSECTS = "T********"

#=
## `RelateMatrixPredicate` (RelateMatrixPredicate.java)

Evaluates the full relate intersection matrix: it is never determined
early, so the entire matrix is computed. `result_im` returns the
accumulated matrix (which may be only partially complete before
`finish!` has been called).
=#
struct RelateMatrixPred end
RelateMatrixPredicate() = IMPredicate(RelateMatrixPred())
predicate_name(::IMPredicate{RelateMatrixPred}) = "relateMatrix"
# ensure entire matrix is computed
require_interaction(::Type{IMPredicate{RelateMatrixPred}}) = false
# ensure entire matrix is computed
is_determined(::IMPredicate{RelateMatrixPred}) = false
# indicates full matrix is being evaluated
value_im(::IMPredicate{RelateMatrixPred}) = false
# Gets the current state of the IM matrix (JTS `getIM`).
result_im(p::IMPredicate{RelateMatrixPred}) = p.im
