module UnitSpherical

using CoordinateTransformations
using StaticArrays, LinearAlgebra
import GeoInterface as GI, GeoFormatTypes as GFT

import Random

# using TestItems # this is a thin package that allows TestItems.@testitem to be parsed.

include("point.jl")
include("coordinate_transforms.jl")
include("slerp.jl")
include("cap.jl")

export UnitSphericalPoint, UnitSphereFromGeographic, GeographicFromUnitSphere, 
       slerp, SphericalCap
export spherical_distance

end