# # RelateNG polygon node converter
#
# Port of JTS `PolygonNodeConverter.java`. Method order parallels the Java
# file (`convert`, `convertShellAndHoles`, `convertHoles`, `createSection`,
# `extractUnique`, `next`, `findShell`), so this file diffs against its Java
# counterpart. Indices are 1-based throughout (Java's `findShell` returns -1
# for "no shell"; here `_find_shell` returns 0).

"""
    polygon_node_convert(m::Manifold, poly_sections::Vector{<:NodeSection}; exact)

Converts the node sections at a polygon node where a shell and one or more
holes touch, or two or more holes touch. This converts the node topological
structure from the OGC "touching-rings" (AKA "minimal-ring") model to the
equivalent "self-touch" (AKA "inverted/exverted ring" or "maximal ring")
model. In the "self-touch" model the converted [`NodeSection`](@ref) corners
enclose areas which all lie inside the polygon (i.e. they do not enclose
hole edges). This allows `RelateNode` (Task 17) to use simple area-additive
semantics for adding edges and propagating edge locations.

The input node sections are assumed to have canonical orientation (CW shells
and CCW holes). The arrangement of shells and holes must be topologically
valid. Specifically, the node sections must not cross or be collinear.

This supports multiple shell-shell touches (including ones containing
holes), and hole-hole touches. This generalizes the relate algorithm to
support both the OGC model and the self-touch model.

Converts a list of sections of valid polygon rings to have "self-touching"
structure. There are the same number of output sections as input ones.
Sorts (and thereby mutates) `poly_sections`; returns the converted sections.

Port of `PolygonNodeConverter.convert`. The angle sort goes through
[`edge_angle_compare`](@ref) (the `NodeSection.EdgeAngleComparator` port),
which is why the manifold and the `exact` flag are threaded in (the Java
method is geometry-context-free).
"""
function polygon_node_convert(m::Manifold, poly_sections::Vector{<:NodeSection}; exact)
    # Stable, like Java List.sort: equal-angle sections (duplicates) keep
    # their input (prepareSections) order, so extractUnique sees them adjacent.
    sort!(poly_sections; alg = MergeSort,
        lt = (a, b) -> edge_angle_compare(m, a, b; exact) < 0)

    #TODO: move uniquing up to caller
    sections = _extract_unique(poly_sections)
    length(sections) == 1 && return sections

    #-- find shell section index
    shell_index = _find_shell(sections)
    if shell_index == 0
        return _convert_holes(sections)
    end
    #-- at least one shell is present.  Handle multiple ones if present
    converted_sections = NodeSection[]
    next_shell_index = shell_index
    while true  # Java do-while
        next_shell_index = _convert_shell_and_holes(sections, next_shell_index, converted_sections)
        next_shell_index == shell_index && break
    end
    return converted_sections
end

# Port of PolygonNodeConverter.convertShellAndHoles: walk CCW from the shell
# section at `shell_index`, closing each shell-hole / hole-hole / hole-shell
# corner into a self-touch shell corner. Returns the index of the next shell
# section (the do-while cursor in `polygon_node_convert`).
function _convert_shell_and_holes(sections::Vector{<:NodeSection}, shell_index::Integer,
        converted_sections::Vector{NodeSection})
    shell_section = sections[shell_index]
    in_vertex = get_vertex(shell_section, 0)
    i = _next(sections, shell_index)
    while !is_shell(sections[i])
        hole_section = sections[i]
        # Assert: is_shell(hole_section) == false
        out_vertex = get_vertex(hole_section, 1)
        ns = _create_section(shell_section, in_vertex, out_vertex)
        push!(converted_sections, ns)

        in_vertex = get_vertex(hole_section, 0)
        i = _next(sections, i)
    end
    #-- create final section for corner from last hole to shell
    out_vertex = get_vertex(shell_section, 1)
    ns = _create_section(shell_section, in_vertex, out_vertex)
    push!(converted_sections, ns)
    return i
end

# Port of PolygonNodeConverter.convertHoles: no shell at the node, so each
# pair of angularly-adjacent hole sections contributes one self-touch corner.
function _convert_holes(sections::Vector{<:NodeSection})
    converted_sections = NodeSection[]
    copy_section = sections[1]
    for i in eachindex(sections)
        inext = _next(sections, i)
        in_vertex = get_vertex(sections[i], 0)
        out_vertex = get_vertex(sections[inext], 1)
        ns = _create_section(copy_section, in_vertex, out_vertex)
        push!(converted_sections, ns)
    end
    return converted_sections
end

# Port of PolygonNodeConverter.createSection: a shell (ring id 0) area
# section of the same polygon as `ns`, with the given edge vertices.
_create_section(ns::NodeSection, v0, v1) = NodeSection(
    is_a(ns), DIM_A, section_id(ns), Int32(0), get_polygonal(ns),
    is_node_at_vertex(ns), v0, node_pt(ns), v1)

# Port of PolygonNodeConverter.extractUnique: drop consecutive duplicate
# sections (the list is sorted, so duplicates are adjacent).
function _extract_unique(sections::Vector{S}) where {S <: NodeSection}
    unique_sections = S[]
    last_unique = sections[1]
    push!(unique_sections, last_unique)
    for ns in sections
        if compare_to(last_unique, ns) != 0
            push!(unique_sections, ns)
            last_unique = ns
        end
    end
    return unique_sections
end

# Port of PolygonNodeConverter.next: circular increment (1-based).
function _next(ns::Vector{<:NodeSection}, i::Integer)
    next = i + 1
    if next > length(ns)
        next = 1
    end
    return next
end

# Port of PolygonNodeConverter.findShell: index of the first shell section,
# or 0 if there is none (Java returns -1).
function _find_shell(poly_sections::Vector{<:NodeSection})
    for i in eachindex(poly_sections)
        is_shell(poly_sections[i]) && return i
    end
    return 0
end
