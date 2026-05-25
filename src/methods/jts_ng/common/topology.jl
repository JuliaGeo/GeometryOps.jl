# # Shared JTS NG topology vocabulary

"""
    TopologicalLocation

The DE-9IM location of a point relative to a geometry.
"""
@enum TopologicalLocation::Int8 loc_interior=0 loc_boundary=1 loc_exterior=2

"""
    TopologicalDimension

The topological dimension values stored in a DE-9IM matrix.
"""
@enum TopologicalDimension::Int8 dim_false=-1 dim_point=0 dim_line=1 dim_area=2

"""
    SidePosition

The side of a directed edge used by NG topology code.
"""
@enum SidePosition::Int8 side_right=-1 side_on=0 side_left=1

const _TOPOLOGICAL_LOCATIONS = (loc_interior, loc_boundary, loc_exterior)
const _TOPOLOGICAL_DIMENSIONS = (dim_false, dim_point, dim_line, dim_area)

location_index(loc::TopologicalLocation) = Int(loc) + 1
dimension_value(dim::TopologicalDimension) = Int(dim)

is_false_dimension(dim::TopologicalDimension) = dim == dim_false
is_true_dimension(dim::TopologicalDimension) = dim != dim_false

function max_dimension(a::TopologicalDimension, b::TopologicalDimension)
    return dimension_value(a) >= dimension_value(b) ? a : b
end

function dimension_char(dim::TopologicalDimension)
    if dim == dim_false
        return 'F'
    elseif dim == dim_point
        return '0'
    elseif dim == dim_line
        return '1'
    elseif dim == dim_area
        return '2'
    else
        error("Unknown topological dimension: $dim")
    end
end

function dimension_from_char(c::Char)
    uc = uppercase(c)
    if uc == 'F'
        return dim_false
    elseif uc == '0'
        return dim_point
    elseif uc == '1'
        return dim_line
    elseif uc == '2'
        return dim_area
    else
        throw(ArgumentError("Expected one of F, 0, 1, or 2, got '$c'."))
    end
end

"""
    DimensionLocation

Packed dimension/location vocabulary shared by the NG engines.  This does not
carry OverlayNG label semantics such as collapse state or effective boundary
state.
"""
struct DimensionLocation
    dimension::TopologicalDimension
    location::TopologicalLocation
end

"""
    BoundaryNodeRule

Abstract supertype for JTS-style line boundary node rules.
"""
abstract type BoundaryNodeRule end

struct Mod2BoundaryNodeRule <: BoundaryNodeRule end
struct EndpointBoundaryNodeRule <: BoundaryNodeRule end
struct MultivalentEndpointBoundaryNodeRule <: BoundaryNodeRule end
struct MonovalentEndpointBoundaryNodeRule <: BoundaryNodeRule end

is_in_boundary(::Mod2BoundaryNodeRule, boundary_count::Integer) = isodd(boundary_count)
is_in_boundary(::EndpointBoundaryNodeRule, boundary_count::Integer) = boundary_count > 0
is_in_boundary(::MultivalentEndpointBoundaryNodeRule, boundary_count::Integer) = boundary_count > 1
is_in_boundary(::MonovalentEndpointBoundaryNodeRule, boundary_count::Integer) = boundary_count == 1
