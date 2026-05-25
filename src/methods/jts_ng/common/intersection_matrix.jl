# # DE-9IM matrix substrate

const _DE9IM_SIZE = 9

"""
    IntersectionMatrix()

Mutable DE-9IM matrix storing the maximum known dimension for each
interior/boundary/exterior cell.
"""
mutable struct IntersectionMatrix
    cells::StaticArrays.MVector{_DE9IM_SIZE, TopologicalDimension}
end

function IntersectionMatrix()
    return IntersectionMatrix(
        StaticArrays.MVector{_DE9IM_SIZE, TopologicalDimension}(
            ntuple(_ -> dim_false, _DE9IM_SIZE),
        ),
    )
end

function IntersectionMatrix(pattern::AbstractString)
    length(pattern) == _DE9IM_SIZE ||
        throw(ArgumentError("A DE-9IM matrix string must have 9 characters."))

    matrix = IntersectionMatrix()
    for (i, c) in enumerate(pattern)
        matrix.cells[i] = dimension_from_char(c)
    end
    return matrix
end

function _matrix_index(row::TopologicalLocation, col::TopologicalLocation)
    return (location_index(row) - 1) * 3 + location_index(col)
end

Base.getindex(matrix::IntersectionMatrix, row::TopologicalLocation, col::TopologicalLocation) =
    matrix.cells[_matrix_index(row, col)]

function Base.setindex!(
    matrix::IntersectionMatrix,
    dimension::TopologicalDimension,
    row::TopologicalLocation,
    col::TopologicalLocation,
)
    matrix.cells[_matrix_index(row, col)] = dimension
    return matrix
end

function set_at_least!(
    matrix::IntersectionMatrix,
    row::TopologicalLocation,
    col::TopologicalLocation,
    dimension::TopologicalDimension,
)
    idx = _matrix_index(row, col)
    matrix.cells[idx] = max_dimension(matrix.cells[idx], dimension)
    return matrix
end

function set_at_least!(matrix::IntersectionMatrix, other::IntersectionMatrix)
    for i in eachindex(matrix.cells)
        matrix.cells[i] = max_dimension(matrix.cells[i], other.cells[i])
    end
    return matrix
end

function de9im_string(matrix::IntersectionMatrix)
    return String(collect(dimension_char.(Tuple(matrix.cells))))
end

Base.copy(matrix::IntersectionMatrix) =
    IntersectionMatrix(StaticArrays.MVector{_DE9IM_SIZE, TopologicalDimension}(Tuple(matrix.cells)))

Base.show(io::IO, matrix::IntersectionMatrix) = print(io, de9im_string(matrix))

function _matches_dimension(dimension::TopologicalDimension, pattern::Char)
    c = uppercase(pattern)
    if c == '*'
        return true
    elseif c == 'T'
        return is_true_dimension(dimension)
    elseif c == 'F'
        return is_false_dimension(dimension)
    elseif c in ('0', '1', '2')
        return dimension == dimension_from_char(c)
    else
        throw(ArgumentError("Invalid DE-9IM pattern character '$pattern'."))
    end
end

function matches(matrix::IntersectionMatrix, pattern::AbstractString)
    length(pattern) == _DE9IM_SIZE ||
        throw(ArgumentError("A DE-9IM pattern must have 9 characters."))

    return all(_matches_dimension(dim, pat) for (dim, pat) in zip(matrix.cells, pattern))
end
