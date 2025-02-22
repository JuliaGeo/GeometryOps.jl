"""
    UnitSphericalPoint(x, y, z)

Represents a point on the unit 2-sphere in 3-D cartesian space, i.e., (x, y, z).
"""
struct UnitSphericalPoint{T <: Number}
	data::NTuple{3, T}
end
UnitSphericalPoint(x, y, z) = UnitSphericalPoint((x, y, z))

# Define GeoInterface, but with lat long coordinates
# This is so that we can write geometries that have USPs in them 
# and quickly at that!

GI.isgeometry(::UnitSphericalPoint) = true
GI.trait(::UnitSphericalPoint) = GI.PointTrait()
GI.geomtrait(::UnitSphericalPoint) = GI.PointTrait()

# TODO this may be controversial
# discuss with collaborators
GI.ncoord(::UnitSphericalPoint) = 2
function GI.getcoord(p::UnitSphericalPoint, i::Integer)
    x, y, z = p.data
    if i == 1
        return atand(y, x)
    elseif i == 2
        return asind(z)
    end
end

# Define ancillary GeoInterface functions like crs, coordtype, crstrait, etc.
GI.crs(::UnitSphericalPoint) = UNIT_SPHERICAL_CRS
GI.crstrait(::UnitSphericalPoint) = UnitSphericalTrait()
# GI.coordtype(::UnitSphericalPoint{T}) where T = T # assume T is a float always, for sanity's sake


# define the 4 basic mathematical operators elementwise on the data tuple
Base.:+(p::UnitSphericalPoint, q::UnitSphericalPoint) = UnitSphericalPoint(p.data .+ q.data)
Base.:-(p::UnitSphericalPoint, q::UnitSphericalPoint) = UnitSphericalPoint(p.data .- q.data)
Base.:*(p::UnitSphericalPoint, q::UnitSphericalPoint) = UnitSphericalPoint(p.data .* q.data)
Base.:/(p::UnitSphericalPoint, q::UnitSphericalPoint) = UnitSphericalPoint(p.data ./ q.data)
# Define sum on a UnitSphericalPoint to sum across its data
Base.sum(p::UnitSphericalPoint) = sum(p.data)

# define dot and cross products
LinearAlgebra.dot(p::UnitSphericalPoint, q::UnitSphericalPoint) = sum(p * q)
function LinearAlgebra.cross(a::UnitSphericalPoint, b::UnitSphericalPoint)
	a1, a2, a3 = a.data
    b1, b2, b3 = b.data
	UnitSphericalPoint((a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1))
end