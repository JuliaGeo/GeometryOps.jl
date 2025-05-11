abstract type AbstractPreparation end

abstract type AbstractPreparationTrait end

function preptrait end

const PrepTuple = Tuple{Vararg{<: AbstractPreparation}}

struct SpatialIndex{IndexType} <: AbstractPreparation
    index::IndexType
end

struct SpatialEdgeIndex{IndexType} <: AbstractPreparation
    index::IndexType
end

