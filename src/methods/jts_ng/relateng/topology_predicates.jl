# # RelateNG topology predicate strategies

abstract type TopologyPredicate end

const RELATE_PREDICATE_UNKNOWN = nothing

mutable struct RelateInteractionPredicate <: TopologyPredicate
    name::Symbol
    known_value::Union{Nothing,Bool}
end

mutable struct RelateIMPredicate <: TopologyPredicate
    name::Symbol
    matrix::IntersectionMatrix
    dim_a::TopologicalDimension
    dim_b::TopologicalDimension
    known_value::Union{Nothing,Bool}
end

mutable struct RelatePatternPredicate <: TopologyPredicate
    pattern::String
    matrix::IntersectionMatrix
    dim_a::TopologicalDimension
    dim_b::TopologicalDimension
    known_value::Union{Nothing,Bool}
end

mutable struct RelateMatrixPredicate <: TopologyPredicate
    matrix::IntersectionMatrix
    dim_a::TopologicalDimension
    dim_b::TopologicalDimension
end

function _relate_im_matrix()
    matrix = IntersectionMatrix()
    matrix[loc_exterior, loc_exterior] = dim_area
    return matrix
end

relate_intersects_predicate() = RelateInteractionPredicate(:intersects, RELATE_PREDICATE_UNKNOWN)
relate_disjoint_predicate() = RelateInteractionPredicate(:disjoint, RELATE_PREDICATE_UNKNOWN)
relate_contains_predicate() = RelateIMPredicate(:contains, _relate_im_matrix(), dim_false, dim_false, nothing)
relate_within_predicate() = RelateIMPredicate(:within, _relate_im_matrix(), dim_false, dim_false, nothing)
relate_covers_predicate() = RelateIMPredicate(:covers, _relate_im_matrix(), dim_false, dim_false, nothing)
relate_coveredby_predicate() = RelateIMPredicate(:coveredby, _relate_im_matrix(), dim_false, dim_false, nothing)
relate_crosses_predicate() = RelateIMPredicate(:crosses, _relate_im_matrix(), dim_false, dim_false, nothing)
relate_equals_predicate() = RelateIMPredicate(:equals, _relate_im_matrix(), dim_false, dim_false, nothing)
relate_overlaps_predicate() = RelateIMPredicate(:overlaps, _relate_im_matrix(), dim_false, dim_false, nothing)
relate_touches_predicate() = RelateIMPredicate(:touches, _relate_im_matrix(), dim_false, dim_false, nothing)
relate_matches_predicate(pattern::AbstractString) =
    RelatePatternPredicate(String(pattern), _relate_im_matrix(), dim_false, dim_false, nothing)
relate_matrix_predicate() = RelateMatrixPredicate(_relate_im_matrix(), dim_false, dim_false)

predicate_name(predicate::RelateInteractionPredicate) = predicate.name
predicate_name(predicate::RelateIMPredicate) = predicate.name
predicate_name(::RelatePatternPredicate) = :matches
predicate_name(::RelateMatrixPredicate) = :matrix

predicate_is_known(predicate::Union{RelateInteractionPredicate,RelateIMPredicate,RelatePatternPredicate}) =
    !isnothing(predicate.known_value)
predicate_is_known(::RelateMatrixPredicate) = false

function predicate_value(predicate::Union{RelateInteractionPredicate,RelateIMPredicate,RelatePatternPredicate})
    predicate_is_known(predicate) ||
        throw(ArgumentError("Predicate value is not known yet."))
    return predicate.known_value
end

predicate_value(::RelateMatrixPredicate) = false
predicate_matrix(predicate::Union{RelateIMPredicate,RelatePatternPredicate,RelateMatrixPredicate}) =
    predicate.matrix

require_self_noding(predicate::RelateInteractionPredicate) = false
require_self_noding(predicate::TopologyPredicate) = true

require_interaction(predicate::RelateInteractionPredicate) = predicate.name == :intersects
require_interaction(predicate::RelateMatrixPredicate) = false
require_interaction(predicate::RelatePatternPredicate) = _relate_pattern_requires_interaction(predicate.pattern)
require_interaction(predicate::TopologyPredicate) = true

require_covers(predicate::TopologyPredicate, input_side::NGInputSide) = false
require_covers(predicate::RelateIMPredicate, input_side::NGInputSide) =
    predicate.name in (:contains, :covers) ? input_side == input_a :
    predicate.name in (:within, :coveredby) ? input_side == input_b :
    false

require_exterior_check(predicate::RelateInteractionPredicate, input_side::NGInputSide) = false
require_exterior_check(predicate::TopologyPredicate, input_side::NGInputSide) = true
require_exterior_check(predicate::RelateIMPredicate, input_side::NGInputSide) =
    predicate.name in (:contains, :covers) ? input_side == input_b :
    predicate.name in (:within, :coveredby) ? input_side == input_a :
    true

function relate_init_predicate!(
    predicate::TopologyPredicate,
    geom_a::RelateGeometry,
    geom_b::RelateGeometry,
)
    relate_init_dimensions!(predicate, relate_dimension_real(geom_a), relate_dimension_real(geom_b))
    relate_init_extents!(predicate, geom_a.extent, geom_b.extent)
    return predicate
end

relate_init_dimensions!(predicate::RelateInteractionPredicate, dim_a, dim_b) = predicate
relate_init_extents!(predicate::TopologyPredicate, extent_a, extent_b) = predicate

function relate_init_dimensions!(predicate::Union{RelateIMPredicate,RelatePatternPredicate,RelateMatrixPredicate}, dim_a, dim_b)
    predicate.dim_a = dim_a
    predicate.dim_b = dim_b
    if predicate isa RelateIMPredicate
        _relate_apply_dimension_requirements!(predicate)
    end
    return predicate
end

function relate_init_extents!(predicate::RelateInteractionPredicate, extent_a, extent_b)
    _relate_extents_known(extent_a, extent_b) || return predicate
    if predicate.name == :intersects
        _relate_set_value_if_unknown!(predicate, false, !Extents.intersects(extent_a, extent_b))
    elseif predicate.name == :disjoint
        _relate_set_value_if_unknown!(predicate, true, Extents.disjoint(extent_a, extent_b))
    end
    return predicate
end

function relate_init_extents!(predicate::RelateIMPredicate, extent_a, extent_b)
    _relate_extents_known(extent_a, extent_b) || return predicate
    if predicate.name in (:contains, :covers)
        _relate_require!(predicate, _relate_extent_covers(extent_a, extent_b))
    elseif predicate.name in (:within, :coveredby)
        _relate_require!(predicate, _relate_extent_covers(extent_b, extent_a))
    elseif predicate.name == :equals
        _relate_require!(predicate, extent_a == extent_b)
    end
    return predicate
end

function relate_init_extents!(predicate::RelatePatternPredicate, extent_a, extent_b)
    _relate_extents_known(extent_a, extent_b) || return predicate
    _relate_set_value_if_unknown!(
        predicate,
        false,
        _relate_pattern_requires_interaction(predicate.pattern) && Extents.disjoint(extent_a, extent_b),
    )
    return predicate
end

function _relate_extents_known(extent_a, extent_b)
    return !isnothing(extent_a) && !isnothing(extent_b)
end

function _relate_extent_covers(source, target)
    return source.X[1] <= target.X[1] &&
        source.X[2] >= target.X[2] &&
        source.Y[1] <= target.Y[1] &&
        source.Y[2] >= target.Y[2]
end

function relate_update_dimension!(
    predicate::RelateInteractionPredicate,
    loc_a::TopologicalLocation,
    loc_b::TopologicalLocation,
    dimension::TopologicalDimension,
)
    if loc_a != loc_exterior && loc_b != loc_exterior && is_true_dimension(dimension)
        if predicate.name == :intersects
            _relate_set_value_if_unknown!(predicate, true)
        elseif predicate.name == :disjoint
            _relate_set_value_if_unknown!(predicate, false)
        end
    end
    return predicate
end

function relate_update_dimension!(
    predicate::Union{RelateIMPredicate,RelatePatternPredicate,RelateMatrixPredicate},
    loc_a::TopologicalLocation,
    loc_b::TopologicalLocation,
    dimension::TopologicalDimension,
)
    previous = predicate.matrix[loc_a, loc_b]
    dimension_value(dimension) > dimension_value(previous) || return predicate

    predicate.matrix[loc_a, loc_b] = dimension
    if predicate isa Union{RelateIMPredicate,RelatePatternPredicate} &&
            _relate_predicate_determined(predicate)
        _relate_set_value_if_unknown!(predicate, _relate_value_from_matrix(predicate))
    end
    return predicate
end

function relate_finish!(predicate::RelateInteractionPredicate)
    if predicate.name == :intersects
        _relate_set_value_if_unknown!(predicate, false)
    elseif predicate.name == :disjoint
        _relate_set_value_if_unknown!(predicate, true)
    end
    return predicate
end

function relate_finish!(predicate::Union{RelateIMPredicate,RelatePatternPredicate})
    _relate_set_value_if_unknown!(predicate, _relate_value_from_matrix(predicate))
    return predicate
end

relate_finish!(predicate::RelateMatrixPredicate) = predicate

function _relate_set_value_if_unknown!(predicate, value::Bool, condition::Bool = true)
    condition || return predicate
    predicate_is_known(predicate) && return predicate
    predicate.known_value = value
    return predicate
end

_relate_require!(predicate, condition::Bool) =
    _relate_set_value_if_unknown!(predicate, false, !condition)

function _relate_apply_dimension_requirements!(predicate::RelateIMPredicate)
    if predicate.name in (:contains, :covers)
        _relate_require!(predicate, _relate_dims_compatible_with_covers(predicate.dim_a, predicate.dim_b))
    elseif predicate.name in (:within, :coveredby)
        _relate_require!(predicate, _relate_dims_compatible_with_covers(predicate.dim_b, predicate.dim_a))
    elseif predicate.name == :crosses
        both_points_or_areas =
            (predicate.dim_a == dim_point && predicate.dim_b == dim_point) ||
            (predicate.dim_a == dim_area && predicate.dim_b == dim_area)
        _relate_require!(predicate, !both_points_or_areas)
    elseif predicate.name == :overlaps
        _relate_require!(predicate, predicate.dim_a == predicate.dim_b)
    elseif predicate.name == :touches
        _relate_require!(predicate, !(predicate.dim_a == dim_point && predicate.dim_b == dim_point))
    end
    return predicate
end

function _relate_dims_compatible_with_covers(dim0::TopologicalDimension, dim1::TopologicalDimension)
    dim0 == dim_point && dim1 == dim_line && return true
    return dimension_value(dim0) >= dimension_value(dim1)
end

function _relate_predicate_determined(predicate::RelateIMPredicate)
    if predicate.name in (:contains, :covers)
        return _relate_intersects_exterior_of(predicate, input_a)
    elseif predicate.name in (:within, :coveredby)
        return _relate_intersects_exterior_of(predicate, input_b)
    elseif predicate.name == :crosses
        if predicate.dim_a == dim_line && predicate.dim_b == dim_line
            return dimension_value(predicate.matrix[loc_interior, loc_interior]) > dimension_value(dim_point)
        elseif dimension_value(predicate.dim_a) < dimension_value(predicate.dim_b)
            return _relate_matrix_intersects(predicate, loc_interior, loc_interior) &&
                _relate_matrix_intersects(predicate, loc_interior, loc_exterior)
        elseif dimension_value(predicate.dim_a) > dimension_value(predicate.dim_b)
            return _relate_matrix_intersects(predicate, loc_interior, loc_interior) &&
                _relate_matrix_intersects(predicate, loc_exterior, loc_interior)
        end
    elseif predicate.name == :equals
        return _relate_matrix_intersects(predicate, loc_interior, loc_exterior) ||
            _relate_matrix_intersects(predicate, loc_boundary, loc_exterior) ||
            _relate_matrix_intersects(predicate, loc_exterior, loc_interior) ||
            _relate_matrix_intersects(predicate, loc_exterior, loc_boundary)
    elseif predicate.name == :overlaps
        if predicate.dim_a in (dim_area, dim_point)
            return _relate_matrix_intersects(predicate, loc_interior, loc_interior) &&
                _relate_matrix_intersects(predicate, loc_interior, loc_exterior) &&
                _relate_matrix_intersects(predicate, loc_exterior, loc_interior)
        elseif predicate.dim_a == dim_line
            return predicate.matrix[loc_interior, loc_interior] == dim_line &&
                _relate_matrix_intersects(predicate, loc_interior, loc_exterior) &&
                _relate_matrix_intersects(predicate, loc_exterior, loc_interior)
        end
    elseif predicate.name == :touches
        return _relate_matrix_intersects(predicate, loc_interior, loc_interior)
    end
    return false
end

function _relate_predicate_determined(predicate::RelatePatternPredicate)
    for row in (loc_interior, loc_boundary, loc_exterior)
        for col in (loc_interior, loc_boundary, loc_exterior)
            pattern_char = predicate.pattern[_matrix_index(row, col)]
            pattern_char == '*' && continue

            matrix_dim = predicate.matrix[row, col]
            if uppercase(pattern_char) == 'T'
                is_true_dimension(matrix_dim) || continue
            elseif dimension_value(matrix_dim) > dimension_value(dimension_from_char(pattern_char))
                return true
            end
        end
    end
    return false
end

function _relate_intersects_exterior_of(predicate::RelateIMPredicate, input_side::NGInputSide)
    if input_side == input_a
        return _relate_matrix_intersects(predicate, loc_exterior, loc_interior) ||
            _relate_matrix_intersects(predicate, loc_exterior, loc_boundary)
    else
        return _relate_matrix_intersects(predicate, loc_interior, loc_exterior) ||
            _relate_matrix_intersects(predicate, loc_boundary, loc_exterior)
    end
end

_relate_matrix_intersects(predicate, row, col) = is_true_dimension(predicate.matrix[row, col])

function _relate_value_from_matrix(predicate::RelatePatternPredicate)
    return matches(predicate.matrix, predicate.pattern)
end

function _relate_value_from_matrix(predicate::RelateIMPredicate)
    matrix = predicate.matrix
    if predicate.name == :contains
        return matches(matrix, "T*****FF*")
    elseif predicate.name == :within
        return matches(matrix, "T*F**F***")
    elseif predicate.name == :covers
        return any(pattern -> matches(matrix, pattern), ("T*****FF*", "*T****FF*", "***T**FF*", "****T*FF*"))
    elseif predicate.name == :coveredby
        return any(pattern -> matches(matrix, pattern), ("T*F**F***", "*TF**F***", "**FT*F***", "**F*TF***"))
    elseif predicate.name == :crosses
        return _relate_matrix_crosses(matrix, predicate.dim_a, predicate.dim_b)
    elseif predicate.name == :equals
        return matches(matrix, "T*F**FFF*")
    elseif predicate.name == :overlaps
        return _relate_matrix_overlaps(matrix, predicate.dim_a)
    elseif predicate.name == :touches
        return !(predicate.dim_a == dim_point && predicate.dim_b == dim_point) &&
            any(pattern -> matches(matrix, pattern), ("FT*******", "F**T*****", "F***T****"))
    end
    return false
end

function _relate_matrix_crosses(matrix::IntersectionMatrix, dim_a, dim_b)
    if dim_a == dim_line && dim_b == dim_line
        return matches(matrix, "0********")
    elseif dimension_value(dim_a) < dimension_value(dim_b)
        return matches(matrix, "T*T******")
    elseif dimension_value(dim_a) > dimension_value(dim_b)
        return matches(matrix, "T*****T**")
    end
    return false
end

function _relate_matrix_overlaps(matrix::IntersectionMatrix, dim)
    dim in (dim_point, dim_area) && return matches(matrix, "T*T***T**")
    dim == dim_line && return matches(matrix, "1*T***T**")
    return false
end

function _relate_pattern_requires_interaction(pattern::AbstractString)
    length(pattern) == _DE9IM_SIZE ||
        throw(ArgumentError("A DE-9IM pattern must have 9 characters."))

    for idx in (_matrix_index(loc_interior, loc_interior),
        _matrix_index(loc_interior, loc_boundary),
        _matrix_index(loc_boundary, loc_interior),
        _matrix_index(loc_boundary, loc_boundary))
        c = uppercase(pattern[idx])
        c == 'T' || c in ('0', '1', '2') || continue
        return true
    end
    return false
end
