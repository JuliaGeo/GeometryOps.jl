#=
# Force dimensions (xy, xyz)

These functions force the geometry to be 2D or 3D.  They work on any geometry, vector of geometries, feature collection, or table!

They're implemented by `apply` pretty simply.
=#

export forcexy, forcexyz

"""
    forcexy(geom)

Force the geometry to be 2D.  Works on any geometry, vector of geometries, feature collection, or table!
"""
function forcexy(geom)
    return apply(GI.PointTrait(), geom) do point
        (GI.x(point), GI.y(point))
    end
end

"""
    forcexyz(geom, z = 0)

Force the geometry to be 3D.  Works on any geometry, vector of geometries, feature collection, or table!

The `z` parameter is the default z value - if a point has no z value, it will be set to this value.  
If it does, then the z value will be kept.
"""
function forcexyz(geom, z = 0)
    return apply(GI.PointTrait(), geom) do point
        x, y = GI.x(point), GI.y(point)
        z = GI.is3d(geom) ? GI.z(point) : z
        (x, y, z)
    end
end
