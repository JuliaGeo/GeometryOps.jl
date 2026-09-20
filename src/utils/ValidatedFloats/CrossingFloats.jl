module CrossingFloats

import LinearAlgebra: cross, dot
import StaticArrays: StaticVector, SVector

export CrossingFloat, center, radius, isbounded, certify, certified_sign,
       scale_pow2, PowerOfTwo

include("src/crossing_float.jl")
include("src/crossing_arithmetic.jl")
include("src/crossing_vectors.jl")
include("src/exact_vectors.jl")

end
