# # Equals

export equals

#=
## What is equals?

The equals function checks if two geometries are equal. They are equal if they
share the same set of points and edges.

To provide an example, consider these two lines:
```@example cshape
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie

l1 = GI.LineString([(0.0, 0.0), (0.0, 10.0)])
l2 = GI.LineString([(0.0, -10.0), (0.0, 3.0)])
f, a, p = lines(GI.getpoint(l1), color = :blue)
scatter!(GI.getpoint(l1), color = :blue)
lines!(GI.getpoint(l2), color = :orange)
scatter!(GI.getpoint(l2), color = :orange)
```
We can see that the two lines do not share a commen set of points and edges in
the plot, so they are not equal:
```@example cshape
equals(l1, l2)  # returns false
```

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait. This is also used in the
implementation, since it's a lot less work! 

Note that while we need the same set of points and edges, they don't need to be
provided in the same order for polygons. For for example, we need the same set
points for two multipoints to be equal, but they don't have to be saved in the
same order. This requires checking every point against every other point in the
two geometries we are comparing.  
=#

"""
    equals(geom1, geom2)::Bool

Compare two Geometries return true if they are the same geometry.

## Examples
```jldoctest
import GeometryOps as GO, GeoInterface as GI
poly1 = GI.Polygon([[(0,0), (0,5), (5,5), (5,0), (0,0)]])
poly2 = GI.Polygon([[(0,0), (0,5), (5,5), (5,0), (0,0)]])

GO.equals(poly1, poly2)
# output
true
```
"""
equals(geom_a, geom_b) = equals(
    GI.trait(geom_a), geom_a,
    GI.trait(geom_b), geom_b,
)

"""
    equals(::T, geom_a, ::T, geom_b)::Bool

Two geometries of the same type, which don't have a equals function to dispatch
off of should throw an error.
"""
equals(::T, geom_a, ::T, geom_b) where T = error("Cant compare $T yet")

"""
    equals(trait_a, geom_a, trait_b, geom_b)

Two geometries which are not of the same type cannot be equal so they always
return false.
"""
equals(trait_a, geom_a, trait_b, geom_b) = false

"""
    equals(::GI.PointTrait, p1, ::GI.PointTrait, p2)::Bool

Two points are the same if they have the same x and y (and z if 3D) coordinates.
"""
function equals(::GI.PointTrait, p1, ::GI.PointTrait, p2)
    GI.ncoord(p1) == GI.ncoord(p2) || return false
    GI.x(p1) == GI.x(p2) || return false
    GI.y(p1) == GI.y(p2) || return false
    if GI.is3d(p1)
        GI.z(p1) == GI.z(p2) || return false 
    end
    return true
end

"""
    equals(::GI.MultiPointTrait, mp1, ::GI.MultiPointTrait, mp2)::Bool

Two multipoints are equal if they share the same set of points.
"""
function equals(::GI.MultiPointTrait, mp1, ::GI.MultiPointTrait, mp2)
    GI.npoint(mp1) == GI.npoint(mp2) || return false
    for p1 in GI.getpoint(mp1)
        has_match = false  # if point has a matching point in other multipoint
        for p2 in GI.getpoint(mp2)
            if equals(p1, p2)
                has_match = true
                break
            end
        end
        has_match || return false  # if no matching point, can't be equal
    end
    return true  # all points had a match
end

"""
    equals(::T, l1, ::T, l2) where {T<:GI.AbstractCurveTrait} ::Bool

Two curves are equal if they share the same set of points going around the
curve. 
"""
function equals(::T, l1, ::T, l2) where {T<:GI.AbstractCurveTrait}
    # Check line lengths match
    n1 = GI.npoint(l1)
    n2 = GI.npoint(l2)
    # TODO: do we need to account for repeated last point??
    n1 == n2 || return false
    
    # Find first matching point if it exists
    p1 = GI.getpoint(l1, 1)
    offset = nothing
    for i in 1:n2
        if equals(p1, GI.getpoint(l2, i))
            offset = i - 1
            break
        end
    end
    isnothing(offset) && return false

    # Then check all points are the same wrapping around line
    for i in 1:n1
        pi = GI.getpoint(l1, i)
        j = i + offset
        j = j <= n1 ? j : (j - n1)
        pj = GI.getpoint(l2, j)
        equals(pi, pj) || return false
    end
    return true
end

"""
    equals(::GI.PolygonTrait, geom_a, ::GI.PolygonTrait, geom_b)::Bool

Two polygons are equal if they share the same exterior edge and holes.
"""
function equals(::GI.PolygonTrait, geom_a, ::GI.PolygonTrait, geom_b)
    # Check if exterior is equal
    equals(GI.getexterior(geom_a), GI.getexterior(geom_b)) || return false
    # Check if number of holes are equal
    GI.nhole(geom_a) == GI.nhole(geom_b) || return false
    # Check if holes are equal
    for ihole in GI.gethole(geom_a)
        has_match = false
        for jhole in GI.gethole(geom_b)
            if equals(ihole, jhole)
                has_match = true
                break
            end
        end
        has_match || return false
    end
    return true
end

"""
    equals(::GI.PolygonTrait, geom_a, ::GI.PolygonTrait, geom_b)::Bool

Two multipolygons are equal if they share the same set of polygons.
"""
function equals(::GI.MultiPolygonTrait, geom_a, ::GI.MultiPolygonTrait, geom_b)
    # Check if same number of polygons
    GI.npolygon(geom_a) == GI.npolygon(geom_b) || return false
    # Check if each polygon has a matching polygon
    for poly_a in GI.getpolygon(geom_a)
        has_match = false
        for poly_b in GI.getpolygon(geom_b)
            if equals(poly_a, poly_b)
                has_match = true
                break
            end
        end
        has_match || return false
    end
    return true
end