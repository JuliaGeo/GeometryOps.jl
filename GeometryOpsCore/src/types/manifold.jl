
#=
# `Manifold`s

A manifold is mathematically defined as a topological space that resembles Euclidean space locally.

In GeometryOps, this represents the domain or space on which your geometries live.  
Manifolds can be accessible via the crs info of the geometry - OR can be specified explicitly.
For example you may pass planar geometry around using GeoJSON, but because the spec says GeoJSON is only geographic,
GeometryOps will interpret GeoJSON geometries as geographic on WGS84, unless told otherwise.

In GeometryOps (and geodesy more generally), there are three manifolds we care about:
- [`Planar`](@ref): the 2d plane, a completely Euclidean manifold
- [`Spherical`](@ref): the unit sphere, but one where areas are multiplied by the radius of the Earth.  This is not Euclidean globally, but all map projections attempt to represent the sphere on the Euclidean 2D plane to varying degrees of success.
- [`Geodesic`](@ref): the ellipsoid, the closest we can come to representing the Earth by a simple geometric shape.  Parametrized by `semimajor_axis` and `inv_flattening`.
- [`AutoManifold`](@ref): a special manifold that automatically selects the best manifold for the operation when it's executed.  Resolves to [`Planar`](@ref), [`Spherical`](@ref), or [`Geodesic`](@ref) depending on the input geometry.

Generally, we aim to have `Linear` and `Spherical` be operable everywhere, whereas `Geodesic` will only apply in specific circumstances.
Currently, those circumstances are [`area`](@ref), `arclength`, and [`segmentize`](@ref), but this could be extended with time and https://github.com/JuliaGeo/SphericalGeodesics.jl.
=#

export Manifold, AutoManifold, Planar, Spherical, Geodesic

"""
    abstract type Manifold

A manifold is mathematically defined as a topological space that resembles Euclidean space locally.

We use the manifold definition to define the space in which an operation should be performed, or where a geometry lies.

Currently we have [`Planar`](@ref), [`Spherical`](@ref), and [`Geodesic`](@ref) manifolds.
"""
abstract type Manifold end

"""
    AutoManifold()

The `AutoManifold` is a special manifold that automatically selects the best manifold for the operation.
It does not carry any parameters, nor does it indicate anything about the nature of the space.

This gets resolved to a specific manifold when an operation is applied, using the `format` method.  
"""
struct AutoManifold <: Manifold end

"""
    Planar()

A planar manifold refers to the 2D Euclidean plane.  

Z coordinates may be accepted but will not influence geometry calculations, which 
are done purely on 2D geometry.  This is the standard "2.5D" model used by e.g. GEOS.
"""
struct Planar <: Manifold
end

"""
    Spherical(; radius)

A spherical manifold means that the geometry is on the 3-sphere (but is represented by 2-D longitude and latitude).  

## Extended help

!!! note
    The traditional definition of spherical coordinates in physics and mathematics, 
    ``r, \\theta, \\phi``, uses the _colatitude_, that measures angular displacement from the `z`-axis.  
    
    Here, we use the geographic definition of longitude and latitude, meaning
    that `lon` is longitude between -180 and 180, and `lat` is latitude between 
    `-90` (south pole) and `90` (north pole).
"""
Base.@kwdef struct Spherical{T} <: Manifold
    radius::T = WGS84_EARTH_MEAN_RADIUS # this should be theWGS84 defined mean radius
end

"""
    Geodesic(; semimajor_axis, inv_flattening)

A geodesic manifold means that the geometry is on a 3-dimensional ellipsoid, parameterized by `semimajor_axis` (``a`` in mathematical parlance)
and `inv_flattening` (``1/f``).

Usually, this is only relevant for area and segmentization calculations.  It becomes more relevant as one grows closer to the poles (or equator).
"""
Base.@kwdef struct Geodesic{T} <: Manifold
    semimajor_axis::T = WGS84_EARTH_SEMI_MAJOR_RADIUS     # WGS84 by default
    inv_flattening::T = WGS84_EARTH_INV_FLATTENING # WGS84 by default
end


# specifically for manifolds
# not used now but will be used later
abstract type EllipsoidParametrization end

struct SemimajorAxisInvFlattening{T} <: EllipsoidParametrization
    semimajor_axis::T
    inv_flattening::T
end

# this should be the full Ellipsoid parametrization from 
struct FullEllipsoidParametrization{T} <: EllipsoidParametrization
    semimajor_axis::T
    semiminor_axis::T
    inv_flattening::T
end

