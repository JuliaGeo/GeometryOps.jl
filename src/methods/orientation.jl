# """
#     isparallel(line1::LineString, line2::LineString)::Bool

# Return `true` if each segment of `line1` is parallel to the correspondent segment of `line2`

# ## Examples
# ```julia
# import GeoInterface as GI, GeometryOps as GO
# julia> line1 = GI.LineString([(9.170356, 45.477985), (9.164434, 45.482551), (9.166644, 45.484003)])
# GeoInterface.Wrappers.LineString{false, false, Vector{Tuple{Float64, Float64}}, Nothing, Nothing}([(9.170356, 45.477985), (9.164434, 45.482551), (9.166644, 45.484003)], nothing, nothing)

# julia> line2 = GI.LineString([(9.169356, 45.477985), (9.163434, 45.482551), (9.165644, 45.484003)])
# GeoInterface.Wrappers.LineString{false, false, Vector{Tuple{Float64, Float64}}, Nothing, Nothing}([(9.169356, 45.477985), (9.163434, 45.482551), (9.165644, 45.484003)], nothing, nothing)

# julia> 
# GO.isparallel(line1, line2)
# true
# ```
# """
# function isparallel(line1, line2)::Bool
#     seg1 = linesegment(line1)
#     seg2 = linesegment(line2)

#     for i in eachindex(seg1)
#         coors2 = nothing
#         coors1 = seg1[i]
#         coors2 = seg2[i]
#         _isparallel(coors1, coors2) == false && return false
#     end
#     return true
# end

# @inline function _isparallel(p1, p2)
#     slope1 = bearing_to_azimuth(rhumb_bearing(GI.x(p1), GI.x(p2)))
#     slope2 = bearing_to_azimuth(rhumb_bearing(GI.y(p1), GI.y(p2)))

#     return slope1 === slope2
# end

_isparallel(((ax, ay), (bx, by)), ((cx, cy), (dx, dy))) = 
    _isparallel(bx - ax, by - ay, dx - cx, dy - cy)

_isparallel(Δx1, Δy1, Δx2, Δy2) = (Δx1 * Δy2 == Δy1 * Δx2)  