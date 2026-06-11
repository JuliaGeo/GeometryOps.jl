# # RelateNG point location
#
# Point-location machinery for RelateNG. This file holds the ports of three
# small, tightly coupled JTS classes, in this order (JTS file boundaries
# preserved as clearly marked sections):
#
# 1. `LinearBoundary`        (JTS LinearBoundary.java)
# 2. `AdjacentEdgeLocator`   (JTS AdjacentEdgeLocator.java) — Task 11
# 3. `RelatePointLocator`    (JTS RelatePointLocator.java)  — Task 12

#==========================================================================
# LinearBoundary (port of JTS LinearBoundary.java)
==========================================================================#

"""
    LinearBoundary(lines, rule::BoundaryNodeRule)

Determines the boundary points of a linear geometry, using a
[`BoundaryNodeRule`](@ref). `lines` is an iterable of linestrings
(any GeoInterface linestring-like geometries); the endpoint degree of
every line endpoint is counted and the rule decides which degrees are
boundary points.

Coordinate keys are normalized via `_node_point` (kernel.jl): exact
`(Float64, Float64)` tuples with signed zeros normalized (`-0.0 → +0.0`),
so lookups here agree with the `NodeKey` vertex-node identity from the
kernel (Task 7) under Dict bit-pattern hashing.

Faithful to Java: only *empty* lines are skipped. Closed lines are NOT
special-cased — a closed line contributes degree 2 to its closure vertex
(both endpoints coincide), which is never a boundary under the Mod-2 or
monovalent rules but would be under e.g. the endpoint rule.
"""
struct LinearBoundary{BR <: BoundaryNodeRule, P}
    vertex_degree::Dict{P, Int}
    has_boundary::Bool
    rule::BR
end

function LinearBoundary(lines, rule::BoundaryNodeRule)
    # assert: dim(geom) == 1
    vertex_degree = _compute_boundary_points(lines)
    has_boundary = _check_boundary(vertex_degree, rule)
    return LinearBoundary(vertex_degree, has_boundary, rule)
end

function _check_boundary(vertex_degree::Dict, rule::BoundaryNodeRule)
    for degree in values(vertex_degree)
        if is_in_boundary(rule, degree)
            return true
        end
    end
    return false
end

has_boundary(lb::LinearBoundary) = lb.has_boundary

function is_boundary(lb::LinearBoundary, pt)
    key = _node_point(pt)
    haskey(lb.vertex_degree, key) || return false
    degree = lb.vertex_degree[key]
    return is_in_boundary(lb.rule, degree)
end

function _compute_boundary_points(lines)
    vertex_degree = Dict{Tuple{Float64, Float64}, Int}()
    for line in lines
        n = GI.npoint(line)
        n == 0 && continue
        _add_endpoint!(_node_point(GI.getpoint(line, 1)), vertex_degree)
        _add_endpoint!(_node_point(GI.getpoint(line, n)), vertex_degree)
    end
    return vertex_degree
end

function _add_endpoint!(p, degree::Dict)
    dim = get(degree, p, 0)
    dim += 1
    degree[p] = dim
    return nothing
end
