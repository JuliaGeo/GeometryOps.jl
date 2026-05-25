using XML

import WellKnownGeometry
import GeoFormatTypes as GFT
import GeometryOps as GO
import GeoInterface as GI

export JTSOperation,
    JTSCase,
    JTSTestSet,
    JTSEmptyGeometry,
    JTSRawGeometry,
    jts_wkt_to_geom,
    load_test_set,
    load_test_cases,
    fixture_family,
    geometry_category,
    case_category,
    is_overlay_operation,
    is_relate_operation

const OVERLAY_OPERATION_NAMES = Set([
    "intersection",
    "intersectionng",
    "union",
    "unionng",
    "difference",
    "differenceng",
    "symdifference",
    "symdifferenceng",
])

const RELATE_OPERATION_NAMES = Set([
    "relate",
    "intersects",
    "contains",
    "covers",
    "coveredby",
    "within",
    "touches",
    "crosses",
    "overlaps",
    "disjoint",
    "equalstopo",
])

const _WKT_PREFIXES = (
    "POINT",
    "LINESTRING",
    "LINEARRING",
    "POLYGON",
    "MULTIPOINT",
    "MULTILINESTRING",
    "MULTIPOLYGON",
    "GEOMETRYCOLLECTION",
)

struct JTSEmptyGeometry
    wkt::String
end

struct JTSRawGeometry
    wkt::String
    parse_error::String
end

"""
    jts_wkt_to_geom(wkt::AbstractString)

Convert a JTS WKT string to a GeometryOps geometry via WellKnownGeometry.jl and
`GO.tuples`.  JTS fixtures often put newlines inside WKT, so this normalizes
whitespace before parsing.
"""
function jts_wkt_to_geom(wkt::AbstractString)
    sanitized_wkt = replace(strip(wkt), r"\s+" => " ")
    sanitized_wkt = replace(sanitized_wkt, r"\(\s+" => "(")
    sanitized_wkt = replace(sanitized_wkt, r"\s+\)" => ")")
    sanitized_wkt = replace(sanitized_wkt, r",\s+" => ",")
    isempty(sanitized_wkt) && return nothing
    _is_simple_empty_wkt(sanitized_wkt) && return JTSEmptyGeometry(sanitized_wkt)
    geom = GFT.WellKnownText(GFT.Geom(), sanitized_wkt)
    try
        return GO.tuples(geom)
    catch err
        occursin(r"\bEMPTY\b"i, sanitized_wkt) || rethrow()
        return JTSRawGeometry(sanitized_wkt, sprint(showerror, err))
    end
end

struct JTSOperation
    name::String
    argument_refs::Vector{String}
    arguments::Vector{Any}
    expected::Any
    expected_text::String
    attributes::Dict{String,String}
end

function Base.show(io::IO, op::JTSOperation)
    print(io, "JTSOperation(name = $(op.name), expected = $(typeof(op.expected)))")
end

struct JTSCase
    description::String
    geom_a::Any
    geom_b::Any
    operations::Vector{JTSOperation}
end

function Base.show(io::IO, case::JTSCase)
    print(io, "JTSCase(description = $(repr(case.description)), operations = $(length(case.operations)))")
end

struct JTSTestSet
    filepath::String
    cases::Vector{JTSCase}
end

load_test_cases(filepath::AbstractString; kwargs...) = load_test_set(filepath; kwargs...).cases

function load_test_set(filepath::AbstractString; operations = nothing)
    doc = read(filepath, XML.Node)
    run = _first_element(doc, "run")
    isnothing(run) && throw(ArgumentError("Expected a <run> root in $filepath."))

    cases = JTSCase[]
    for case_node in _element_children(run, "case")
        case = parse_case(case_node; operations)
        isnothing(case) || push!(cases, case)
    end
    return JTSTestSet(String(filepath), cases)
end

function parse_case(case_node::XML.Node; operations = nothing)
    desc_node = _first_element(case_node, "desc")
    a_node = _first_element(case_node, "a")
    b_node = _first_element(case_node, "b")

    description = isnothing(desc_node) ? "" : _node_text(desc_node)
    geom_a = isnothing(a_node) ? nothing : jts_wkt_to_geom(_node_text(a_node))
    geom_b = isnothing(b_node) ? nothing : jts_wkt_to_geom(_node_text(b_node))

    parsed_operations = JTSOperation[]
    for test_node in _element_children(case_node, "test")
        for op_node in _element_children(test_node, "op")
            op = parse_operation(op_node, geom_a, geom_b)
            _operation_allowed(op.name, operations) && push!(parsed_operations, op)
        end
    end

    isempty(parsed_operations) && !isnothing(operations) && return nothing
    return JTSCase(description, geom_a, geom_b, parsed_operations)
end

function parse_operation(op_node::XML.Node, geom_a, geom_b)
    attrs = Dict(String(k) => String(v) for (k, v) in XML.attributes(op_node))
    name = attrs["name"]
    argument_refs = _argument_refs(attrs)
    arguments = Any[_resolve_argument(ref, geom_a, geom_b) for ref in argument_refs]
    expected_text = _node_text(op_node)
    expected = _parse_expected(expected_text)
    return JTSOperation(name, argument_refs, arguments, expected, expected_text, attrs)
end

is_overlay_operation(name::AbstractString) = lowercase(name) in OVERLAY_OPERATION_NAMES
is_overlay_operation(op::JTSOperation) = is_overlay_operation(op.name)

is_relate_operation(name::AbstractString) = lowercase(name) in RELATE_OPERATION_NAMES
is_relate_operation(op::JTSOperation) = is_relate_operation(op.name)

function fixture_family(filepath::AbstractString)
    filename = basename(filepath)
    if occursin("Relate", filename) || occursin("Prepared", filename)
        return :relate
    elseif occursin("Overlay", filename) || occursin("UnaryUnion", filename)
        return :overlay
    else
        return :other
    end
end

function geometry_category(geom)
    isnothing(geom) && return :missing
    try
        GI.isempty(geom) && return :empty
    catch
    end
    return _geometry_category(GI.trait(geom))
end

geometry_category(::JTSEmptyGeometry) = :empty
geometry_category(geom::JTSRawGeometry) = _wkt_category(geom.wkt)

_geometry_category(::GI.PointTrait) = :point
_geometry_category(::GI.MultiPointTrait) = :point
_geometry_category(::GI.LineTrait) = :line
_geometry_category(::GI.LineStringTrait) = :line
_geometry_category(::GI.LinearRingTrait) = :line
_geometry_category(::GI.MultiLineStringTrait) = :line
_geometry_category(::GI.PolygonTrait) = :area
_geometry_category(::GI.MultiPolygonTrait) = :area
_geometry_category(::GI.GeometryCollectionTrait) = :collection
_geometry_category(::GI.FeatureTrait) = :feature
_geometry_category(::GI.AbstractGeometryTrait) = :geometry
_geometry_category(_) = :unknown

function case_category(case::JTSCase)
    a = geometry_category(case.geom_a)
    b = geometry_category(case.geom_b)
    if a == :empty || b == :empty
        return :empty
    elseif a == :collection || b == :collection
        return :collection
    else
        return Symbol(string(a), "_", string(b))
    end
end

function _argument_refs(attrs::Dict{String,String})
    arg_keys = filter(k -> occursin(r"^arg\d+$", k), collect(keys(attrs)))
    sort!(arg_keys, by = k -> parse(Int, k[4:end]))
    return [attrs[k] for k in arg_keys]
end

function _resolve_argument(ref::AbstractString, geom_a, geom_b)
    lower_ref = lowercase(strip(ref))
    if lower_ref == "a"
        return geom_a
    elseif lower_ref == "b"
        return geom_b
    else
        return _parse_literal(ref)
    end
end

function _parse_expected(text::AbstractString)
    stripped = strip(text)
    isempty(stripped) && return nothing
    _is_wkt(stripped) && return jts_wkt_to_geom(stripped)
    return _parse_literal(stripped)
end

function _parse_literal(text::AbstractString)
    stripped = strip(text)
    lower_text = lowercase(stripped)
    if lower_text == "true"
        return true
    elseif lower_text == "false"
        return false
    end

    int_value = tryparse(Int, stripped)
    !isnothing(int_value) && return int_value

    float_value = tryparse(Float64, stripped)
    !isnothing(float_value) && return float_value

    return stripped
end

function _is_wkt(text::AbstractString)
    upper_text = uppercase(strip(text))
    return any(startswith(upper_text, prefix) for prefix in _WKT_PREFIXES)
end

function _is_simple_empty_wkt(text::AbstractString)
    return occursin(
        r"^(POINT|LINESTRING|LINEARRING|POLYGON|MULTIPOINT|MULTILINESTRING|MULTIPOLYGON|GEOMETRYCOLLECTION) EMPTY$"i,
        strip(text),
    )
end

function _wkt_category(text::AbstractString)
    upper_text = uppercase(strip(text))
    _is_simple_empty_wkt(upper_text) && return :empty
    startswith(upper_text, "GEOMETRYCOLLECTION") && return :collection
    startswith(upper_text, "POINT") && return :point
    startswith(upper_text, "MULTIPOINT") && return :point
    startswith(upper_text, "LINESTRING") && return :line
    startswith(upper_text, "LINEARRING") && return :line
    startswith(upper_text, "MULTILINESTRING") && return :line
    startswith(upper_text, "POLYGON") && return :area
    startswith(upper_text, "MULTIPOLYGON") && return :area
    return :unknown
end

function _element_children(node::XML.Node, tag_name::AbstractString)
    return [child for child in children(node) if tag(child) == tag_name]
end

function _first_element(node::XML.Node, tag_name::AbstractString)
    for child in children(node)
        tag(child) == tag_name && return child
        nested = _first_element(child, tag_name)
        isnothing(nested) || return nested
    end
    return nothing
end

function _node_text(node::XML.Node)
    if isempty(children(node))
        return strip(value(node))
    end
    return strip(join((_node_text(child) for child in children(node)), ""))
end

function _operation_allowed(name::AbstractString, operations::Nothing)
    return true
end

function _operation_allowed(name::AbstractString, operation::Union{AbstractString,Symbol,Regex})
    return _operation_allowed(name, (operation,))
end

function _operation_allowed(name::AbstractString, operations)
    return any(operations) do op
        if op isa Regex
            return occursin(op, name)
        else
            return lowercase(String(op)) == lowercase(name)
        end
    end
end
