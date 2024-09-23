using Random, Distributions

"""
    generate_random_poly(
        x,
        y,
        nverts,
        avg_radius,
        irregularity,
        spikiness,
    )
Generate a random polygon given a center point, number of vertices, and measures
of the polygon's irregularity.
Inputs:
    x               <Real> x-coordinate for center point
    y               <Real> y-coordinate for center point
    nverts          <Int> number of vertices
    avg_radius      <Real> average radius for each point to center point
    irregularity    <Real> measure of randomness for difference in angle between
                        points (between 0 and 1)
    spikiness       <Real> measure of randomness for difference in radius
                        between points (between 0 and 1)
    rng             <RNG> random number generator for polygon creation
Output:
    Vector{Vector{Vector{T}}} representing polygon coordinates
Note:
    Check your outputs! No guarantee that the polygon's aren't self-intersecting
"""
function generate_random_poly(
    x,
    y,
    nverts,
    avg_radius,
    irregularity::T,
    spikiness::T,
    rng = Xoshiro()
) where T <: Real
    # Check argument validity
    @assert 0 <= irregularity <= 1 "Irregularity must be between 0 and 1"
    @assert 0 <= spikiness <= 1 "Spikiness must be between 0 and 1"
    # Setup basic parameters
    avg_angle = 2π / nverts
    ϵ_angle = irregularity * avg_angle
    # ϵ_rad = spikiness * avg_radius
    smallest_angle_step = avg_angle - ϵ_angle
    # smallest_rad = avg_radius - ϵ_rad
    current_angle = rand(rng) * 2π
    rad_distribution = Distributions.Normal(avg_radius, spikiness)
    points = [zeros(T, 2) for _ in 1:nverts]
    # Generate angle steps around polygon
    angle_steps = zeros(T, nverts)
    cumsum = T(0)
    for i in 1:nverts
        step = smallest_angle_step + 2ϵ_angle * rand(rng)
        angle_steps[i] = step
        cumsum += step
    end
    angle_steps ./= (cumsum / 2π)
    # Generate polygon points at given angles and radii
    for i in 1:nverts
        rad = clamp(rand(rad_distribution), 0, 2avg_radius) #smallest_rad + 2ϵ_rad * rand(rng)
        points[i][1] = x + rad * cos(current_angle)
        points[i][2] = y + rad * sin(current_angle)
        current_angle -= angle_steps[i]
    end
    points[end] .= points[1]
    return [points]
end
