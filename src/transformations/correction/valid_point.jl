#=
# Valid Points

We define a point as valid if it has no NaN or Inf values.

Other corrections are the domain of the user and can be defined relatively easily using this template.

=#

struct ValidPoint <: GeometryCorrection end

application_level(::ValidPoint) = Union{GI.LineStringTrait, GI.LinearRingTrait}

function (::ValidPoint)(::Union{GI.LineStringTrait, GI.LinearRingTrait}, geom)
    new_coords = NTuple{2, Float64}[]
    sizehint!(new_coords, GI.npoint(geom))
    for coord in GI.getpoint(geom)
        x, y = GI.x(coord), GI.y(coord)
        if isfinite(x) && isfinite(y)
            push!(new_coords, (x, y))
        end
    end
    return rebuild(geom, new_coords)
end