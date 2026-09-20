module ValidatedFloats

export ValidatedFloat, center, radius, isbounded, certify, certified_sign, CrossingFloats

include("CrossingFloats.jl")
include("src/validated_float.jl")
include("src/bounds.jl")
include("src/arithmetic.jl")

end
