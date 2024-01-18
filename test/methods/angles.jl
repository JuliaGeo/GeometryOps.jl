
angles(geom; ::Type{T} = Float64) where T = _angles(T, GI.trait(geom), geom)

_angles(::Type{T}, ::Union{GI.PointTrait, GI.LineTrait}, geom) where T = T[]

function _angles(::Type{T}, ::Union{GI.LineStringTrait, GI.LinearRingTrait}, geom)
    angle_list = Vector{T}(undef, GI.npoints(geom) - 2)

    return
end

function _angles(::Type{T}, ::GI.Polygon, geom)
    angle_list = Vector{T}(undef, GI.npoints(geom) - 1)  # TODO add in check for repeted last coordinate
    local p1, last_diff
    start_diff = (zero(T), zero(T))
    for (i, p2) in enumerate(GI.getpoint(geom))
        if i > 1
            current_diff = (GI.x(p2) - GI.x(p1), GI.y(p2) - GI.y(p1))
            if i == 2
                start_diff = current_diff
            else
                
            end
            last_diff = current_diff
        end
        p1 = p2
    end
    # Calculate needed vectors
    pdiff = diff(coords[1])
    npoints = length(pdiff)
    v1 = -pdiff
    v2 = vcat(pdiff[2:end], pdiff[1:1])
    v1_dot_v2 = [sum(v1[i] .* v2[i]) for i in collect(1:npoints)]
    mag_v1 = sqrt.([sum(v1[i].^2) for i in collect(1:npoints)])
    mag_v2 = sqrt.([sum(v2[i].^2) for i in collect(1:npoints)])
    # Protect against division by 0 caused by very close points
    replace!(mag_v1, 0=>eps(FT))
    replace!(mag_v2, 0=>eps(FT))
    angles = real.(
        acos.(
            clamp!(v1_dot_v2 ./ mag_v1 ./ mag_v2, FT(-1), FT(1))
        ) * 180 / pi
    )

    #= The first angle computed was for the second vertex, and the last was for
    the first vertex. Scroll one position down to make the last vertex be the
    first. =#
    sangles = circshift(angles, 1)
    # Now determine if any vertices are concave and adjust angles accordingly.
    sgn = convex_angle_test(coords[1])
    for i in eachindex(sangles)
        sangles[i] = (sgn[i] == -1) ? (-sangles[i] + 360) : sangles[i]
    end
    return sangles

end

"""
    polyedge(p1, p2)

Outputs the coefficients of the line passing through p1 and p2.
The line is of the form w1x + w2y + w3 = 0. 
Inputs:
    p1 <Vector{Float}> [x, y] point
    p2 <Vector{Float}> [x, y] point
Outputs:
    Three-element vector for coefficents of line passing through p1 and p2
Note:
    See note on calc_poly_angles for credit for this function.
"""
function polyedge(p1::Vector{<:FT}, p2) where FT 
    x1 = p1[1]
    y1 = p1[2]
    x2 = p2[1]
    y2 = p2[2]
    w = if x1 == x2
            [-1/x1, 0, 1]
        elseif y1 == y2
            [0, -1/y1, 1]
        elseif x1 == y1 && x2 == y2
            [1, 1, 0]
        else
            v = (y1 - y2)/(x1*(y2 - y1) - y1*(x2 - x1) + eps(FT))
            [v, -v*(x2 - x1)/(y2 - y1), 1]
        end
    return w
end

"""
    orient_coords(coords)

Take given coordinates and make it so that the first point has the smallest
x-coordiante and so that the coordinates are ordered in a clockwise sequence.
Duplicates vertices will be removed and the coordiantes will be closed (first
and last point are the same).

Input:
    coords  <RingVec> vector of points [x, y]
Output:
    coords  <RingVec> oriented clockwise with smallest x-coordinate first
"""
function orient_coords(coords::RingVec)
    # extreem_idx is point with smallest x-value - if tie, choose lowest y-value
    extreem_idx = 1
    for i in eachindex(coords)
        ipoint = coords[i]
        epoint = coords[extreem_idx]
        if ipoint[1] < epoint[1]
            extreem_idx = i
        elseif ipoint[1] == epoint[1] && ipoint[2] < epoint[2]
            extreem_idx = i
        end
    end
    # extreem point must be first point in list
    new_coords = similar(coords)
    circshift!(new_coords, coords, -extreem_idx + 1)
    valid_ringvec!(new_coords)

    # if coords are counterclockwise, switch to clockwise
    orient_matrix = hcat(
        ones(3),
        vcat(new_coords[1]', new_coords[2]', new_coords[end-1]') # extreem/adjacent points
    )
    if det(orient_matrix) > 0
        reverse!(new_coords)
    end
    return new_coords
end

"""
    convex_angle_test(coords::RingVec{T})

Determine which angles in the polygon are convex, with the assumption that the
first angle is convex, no other vertex has a smaller x-coordinate, and the
vertices are assumed to be ordered in a clockwise sequence. The test is based on
the fact that every convex vertex is on the positive side of the line passing
through the two vertices immediately following each vertex being considered. 
Inputs:
    coords <RingVec{Float}> Vector of [x, y] vectors that make up the exterior
        of a polygon
Outputs:
        sgn <Vector of 1s and -1s> One element for each [x,y] pair - if 1 then
            the angle at that vertex is convex, if it is -1 then the angle is
            concave.
"""
function convex_angle_test(coords::RingVec{T}) where T
    L = 10^25
    # Extreme points used in following loop, apended by a 1 for dot product
    top_left = [-L, -L, 1]
    top_right = [-L, L, 1]
    bottom_left = [L, -L, 1]
    bottom_right = [L, L, 1]
    sgn = [1]  # First vertex is convex

    for k in collect(2:length(coords)-1)
        p1 = coords[k-1]
        p2 = coords[k]  # Testing this point for concavity
        p3 = coords[k+1]
        # Coefficents of polygon edge passing through p1 and p2
        w = polyedge(p1, p2)

        #= Establish the positive side of the line w1x + w2y + w3 = 0.
        The positive side of the line should be in the right side of the vector
        (p2- p3).Δx and Δy give the direction of travel, establishing which of
        the extreme points (see above) should be on the + side. If that point is
        on the negative side of the line, then w is replaced by -w. =#
        Δx = p2[1] - p1[1]
        Δy = p2[2] - p1[2]
        if Δx == Δy == 0
            throw(ArgumentError("Data into convextiy test is 0 or duplicated"))
        end
        vector_product =
            if Δx <= 0 && Δy >= 0  # Bottom_right should be on + side.
                dot(w, bottom_right)
            elseif Δx <= 0 && Δy <=0  # Top_right should be on + side.
                dot(w, top_right)
            elseif Δx>=0 && Δy<=0  # Top_left should be on + side.
                dot(w, top_left)
            else  # Bottom_left should be on + side.
                dot(w, bottom_left)
            end
            w *= sign(vector_product)

            # For vertex at p2 to be convex, p3 has to be on + side of line
            if (w[1]*p3[1] + w[2]*p3[2] + w[3]) < 0
                push!(sgn, -1)
            else
                push!(sgn, 1)
            end
    end
    return sgn
end

"""
    calc_poly_angles(coords::PolyVec{T})

Computes internal polygon angles (in degrees) of an arbitrary simple polygon.
The program eliminates duplicate points, except that the first row must equal
the last, so that the polygon is closed.
Inputs:
    coords  <PolyVec{Float}> coordinates from a polygon
Outputs:
    Vector of polygon's interior angles in degrees

Note - Translated into Julia from the following program (including helper
    functions convex_angle_test and polyedge):
    Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
    Digital Image Processing Using MATLAB, Prentice-Hall, 2004
    Revision: 1.6 Date: 2003/11/21 14:44:06
Warning - Assumes polygon has clockwise winding order. Use orient_coords! to
    update coordinates prior to use
"""
function calc_poly_angles(coords::PolyVec{FT}) where {FT<:AbstractFloat}
    # Calculate needed vectors
    pdiff = diff(coords[1])
    npoints = length(pdiff)
    v1 = -pdiff
    v2 = vcat(pdiff[2:end], pdiff[1:1])
    v1_dot_v2 = [sum(v1[i] .* v2[i]) for i in collect(1:npoints)]
    mag_v1 = sqrt.([sum(v1[i].^2) for i in collect(1:npoints)])
    mag_v2 = sqrt.([sum(v2[i].^2) for i in collect(1:npoints)])
    # Protect against division by 0 caused by very close points
    replace!(mag_v1, 0=>eps(FT))
    replace!(mag_v2, 0=>eps(FT))
    angles = real.(
        acos.(
            clamp!(v1_dot_v2 ./ mag_v1 ./ mag_v2, FT(-1), FT(1))
        ) * 180 / pi
    )

    #= The first angle computed was for the second vertex, and the last was for
    the first vertex. Scroll one position down to make the last vertex be the
    first. =#
    sangles = circshift(angles, 1)
    # Now determine if any vertices are concave and adjust angles accordingly.
    sgn = convex_angle_test(coords[1])
    for i in eachindex(sangles)
        sangles[i] = (sgn[i] == -1) ? (-sangles[i] + 360) : sangles[i]
    end
    return sangles
end