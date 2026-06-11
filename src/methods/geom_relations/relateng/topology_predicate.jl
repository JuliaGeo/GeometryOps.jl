# # Topology predicate framework
#=
Port of JTS `TopologyPredicate`, `BasicPredicate`, `IMPredicate`
(operation/relateng), plus the `intersects` and `disjoint`
`BasicPredicate` kinds from `RelatePredicate.java` (the other named
predicates are `IMPredicate` kinds and live in `relate_predicates.jl`).

Julia has no field inheritance, so JTS's class triangle becomes two
kind-parameterized mutable structs (`BasicPredicate{K}`, `IMPredicate{K}`).
Per-kind behavior (`is_determined`, `value_im`, requirement flags, init
overrides) dispatches on the kind singleton; requirement flags are pure
functions of the *type*, so evaluation specializes per predicate.
=#

#=
## `TopologyPredicate` API (TopologyPredicate.java)

The abstract supertype for strategy types implementing spatial predicates
based on the DE-9IM topology model. Concrete predicates implement:
`predicate_name(p)`, `update_dim!(p, locA, locB, dim)`, `finish!(p)`,
`is_known(p)`, `predicate_value(p)`, and may override the requirement
flags and `init_dims!`/`init_bounds!` hooks below.
=#
abstract type TopologyPredicate end

# Whether the predicate requires self-noding for geometries with
# crossing edges. JTS default: true.
require_self_noding(::Type{<:TopologyPredicate}) = true
# Whether the predicate requires interaction between the input geometries
# (i.e. some entry of IM[I/B, I/B] >= 0). JTS default: true.
require_interaction(::Type{<:TopologyPredicate}) = true
# Whether the predicate requires the source to cover the target.
# JTS default: false.
require_covers(::Type{<:TopologyPredicate}, is_source_a::Bool) = false
# Whether the predicate requires checking if the source input intersects
# the exterior of the target input. JTS default: true.
require_exterior_check(::Type{<:TopologyPredicate}, is_source_a::Bool) = true

# Instance-level forwarding of the requirement flags. Most predicates'
# flags are pure functions of the type (so evaluation can specialize on
# them), but runtime-data-dependent predicates (e.g. `IMPatternMatcher`,
# whose `require_interaction` is computed from its pattern matrix) can
# override these instance methods.
require_self_noding(p::TopologyPredicate) = require_self_noding(typeof(p))
require_interaction(p::TopologyPredicate) = require_interaction(typeof(p))
require_covers(p::TopologyPredicate, is_source_a::Bool) = require_covers(typeof(p), is_source_a)
require_exterior_check(p::TopologyPredicate, is_source_a::Bool) = require_exterior_check(typeof(p), is_source_a)
# Initializes the predicate for a specific geometric case from the input
# dimensions. Default: dimensions provide no information.
init_dims!(p::TopologyPredicate, dimA::Integer, dimB::Integer) = nothing
# Initializes the predicate from the input bounds (JTS `init(Envelope, Envelope)`).
# Default: bounds provide no information.
init_bounds!(p::TopologyPredicate, extA, extB) = nothing

#=
## Tri-state value and `BasicPredicate` (BasicPredicate.java)

The base for relate predicates with a boolean value, with tri-state logic
to detect when the final value has been determined.
=#

const TRI_UNKNOWN = Int8(-1)
const TRI_FALSE   = Int8(0)
const TRI_TRUE    = Int8(1)

# Tests if two geometries intersect based on an interaction at given locations.
is_intersection(locA::Integer, locB::Integer) =
    locA != LOC_EXTERIOR && locB != LOC_EXTERIOR

ext_intersects(extA, extB) = Extents.intersects(extA, extB)
ext_covers(extA, extB) = Extents.covers(extA, extB)

mutable struct BasicPredicate{K} <: TopologyPredicate
    const kind::K
    value::Int8
end
BasicPredicate(kind) = BasicPredicate(kind, TRI_UNKNOWN)

is_known(p::TopologyPredicate) = p.value != TRI_UNKNOWN
predicate_value(p::TopologyPredicate) = p.value == TRI_TRUE

# Updates the predicate value to the given state if it is currently unknown.
# (JTS `setValue`: doesn't change an already-known value.)
function set_value!(p::TopologyPredicate, val::Bool)
    is_known(p) && return nothing
    p.value = val ? TRI_TRUE : TRI_FALSE
    return nothing
end
set_value_if!(p::TopologyPredicate, val::Bool, cond::Bool) =
    cond ? set_value!(p, val) : nothing
require!(p::TopologyPredicate, cond::Bool) =
    cond ? nothing : set_value!(p, false)
require_covers!(p::TopologyPredicate, extA, extB) =
    require!(p, ext_covers(extA, extB))

#=
## `intersects` and `disjoint` kinds (RelatePredicate.java)

These are the only named predicates which are plain `BasicPredicate`s
in JTS (everything else tracks an intersection matrix).
=#

struct IntersectsPred end
struct DisjointPred end

# intersects (RelatePredicate.java `intersects()`)
pred_intersects() = BasicPredicate(IntersectsPred())
predicate_name(::BasicPredicate{IntersectsPred}) = "intersects"
# self-noding is not required to check for a simple interaction
require_self_noding(::Type{BasicPredicate{IntersectsPred}}) = false
# intersects only requires testing interaction
require_exterior_check(::Type{BasicPredicate{IntersectsPred}}, is_source_a::Bool) = false
init_bounds!(p::BasicPredicate{IntersectsPred}, extA, extB) =
    require!(p, ext_intersects(extA, extB))
update_dim!(p::BasicPredicate{IntersectsPred}, locA, locB, dim) =
    set_value_if!(p, true, is_intersection(locA, locB))
# if no intersecting locations were found
finish!(p::BasicPredicate{IntersectsPred}) = set_value!(p, false)

# disjoint (RelatePredicate.java `disjoint()`)
pred_disjoint() = BasicPredicate(DisjointPred())
predicate_name(::BasicPredicate{DisjointPred}) = "disjoint"
# self-noding is not required to check for a simple interaction
require_self_noding(::Type{BasicPredicate{DisjointPred}}) = false
require_interaction(::Type{BasicPredicate{DisjointPred}}) = false
# disjoint only requires testing interaction
require_exterior_check(::Type{BasicPredicate{DisjointPred}}, is_source_a::Bool) = false
init_bounds!(p::BasicPredicate{DisjointPred}, extA, extB) =
    set_value_if!(p, true, !ext_intersects(extA, extB))
update_dim!(p::BasicPredicate{DisjointPred}, locA, locB, dim) =
    set_value_if!(p, false, is_intersection(locA, locB))
# if no intersecting locations were found
finish!(p::BasicPredicate{DisjointPred}) = set_value!(p, true)

#=
## `IMPredicate` core (IMPredicate.java)

The base for predicates which are determined using entries in an
intersection matrix. Each kind must implement `is_determined(p)` and
`value_im(p)`; the kinds themselves are ported in `relate_predicates.jl`.
=#

# allow Points coveredBy zero-length Lines
is_dims_compatible_with_covers(dim0::Integer, dim1::Integer) =
    (dim0 == DIM_P && dim1 == DIM_L) ? true : dim0 >= dim1

const DIM_UNKNOWN = DIM_DONTCARE   # JTS IMPredicate.DIM_UNKNOWN = Dimension.DONTCARE

mutable struct IMPredicate{K} <: TopologyPredicate
    const kind::K
    dimA::Int8
    dimB::Int8
    im::DE9IM
    value::Int8
end
# JTS `IntersectionMatrix()` initializes all entries to `Dimension.FALSE`,
# then the IMPredicate constructor presets E/E, which is always dim = 2.
IMPredicate(kind) = IMPredicate(kind, DIM_UNKNOWN, DIM_UNKNOWN,
    with_entry(DE9IM(), LOC_EXTERIOR, LOC_EXTERIOR, DIM_A), TRI_UNKNOWN)

function init_dims!(p::IMPredicate, dimA::Integer, dimB::Integer)
    p.dimA = dimA
    p.dimB = dimB
    init_dims_kind!(p)   # per-kind hook (JTS subclasses override `init` and call super)
    return nothing
end
init_dims_kind!(p::IMPredicate) = nothing

function update_dim!(p::IMPredicate, locA, locB, dim)
    # only record an increased dimension value
    if is_dim_changed(p, locA, locB, dim)
        p.im = with_entry(p.im, locA, locB, dim)
        # set value if predicate value can be known
        if is_determined(p)
            set_value!(p, value_im(p))
        end
    end
    return nothing
end

is_dim_changed(p::IMPredicate, locA, locB, dim) = dim > p.im[locA, locB]

# Tests whether predicate evaluation can be short-circuited due to the
# current state of the matrix providing enough information to determine
# the predicate value. Implemented per kind.
function is_determined end
# Gets the value of the predicate according to the current intersection
# matrix state. Implemented per kind.
function value_im end

# Tests whether the exterior of the specified input geometry
# is intersected by any part of the other input.
intersects_exterior_of(p::IMPredicate, is_a::Bool) = is_a ?
    (is_intersects_entry(p, LOC_EXTERIOR, LOC_INTERIOR) || is_intersects_entry(p, LOC_EXTERIOR, LOC_BOUNDARY)) :
    (is_intersects_entry(p, LOC_INTERIOR, LOC_EXTERIOR) || is_intersects_entry(p, LOC_BOUNDARY, LOC_EXTERIOR))

is_intersects_entry(p::IMPredicate, locA, locB) = p.im[locA, locB] >= DIM_P
is_known_entry(p::IMPredicate, locA, locB) = p.im[locA, locB] != DIM_UNKNOWN
is_dimension_entry(p::IMPredicate, locA, locB, dim) = p.im[locA, locB] == dim
get_dimension(p::IMPredicate, locA, locB) = p.im[locA, locB]

# Sets the final value based on the state of the IM.
finish!(p::IMPredicate) = set_value!(p, value_im(p))

Base.show(io::IO, p::IMPredicate) =
    print(io, predicate_name(p), ": ", string(p.im))
