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
    JTSPrecisionModel,
    JTSFixtureRule,
    jts_wkt_to_geom,
    find_jts_test_files,
    load_test_set,
    load_test_sets,
    load_test_cases,
    fixture_family,
    geometry_category,
    case_category,
    is_overlay_operation,
    is_relate_operation,
    is_runnable,
    is_skipped,
    is_broken,
    is_unimplemented

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

"""
    JTSEmptyGeometry(wkt)

Placeholder for simple `... EMPTY` WKT forms that fixture tests still need to classify.
"""
struct JTSEmptyGeometry
    wkt::String
end

"""
    JTSRawGeometry(wkt, parse_error)

Lossless placeholder for fixture geometry text that cannot yet be lowered to `GO.tuples`.
"""
struct JTSRawGeometry
    wkt::String
    parse_error::String
end

"""
    JTSPrecisionModel

Parsed `<precisionModel>` metadata from a JTS XML test run.
"""
struct JTSPrecisionModel
    model_type::String
    scale::Union{Nothing,Float64}
    offsetx::Union{Nothing,Float64}
    offsety::Union{Nothing,Float64}
    attributes::Dict{String,String}
end

JTSPrecisionModel(attributes::Dict{String,String}) = JTSPrecisionModel(
    get(attributes, "type", "FLOATING"),
    _parse_optional_float(get(attributes, "scale", nothing)),
    _parse_optional_float(get(attributes, "offsetx", nothing)),
    _parse_optional_float(get(attributes, "offsety", nothing)),
    attributes,
)

"""
    JTSFixtureRule(; file = nothing, case = nothing, operation = nothing, reason = "")

Conformance-status rule for fixture operations.  Matching ops can be tagged as
skipped, known broken, or known unimplemented without editing upstream XML.
"""
struct JTSFixtureRule
    file::Any
    case::Any
    operation::Any
    reason::String
end

JTSFixtureRule(; file = nothing, case = nothing, operation = nothing, reason = "") =
    JTSFixtureRule(file, case, operation, String(reason))

"""
    jts_wkt_to_geom(wkt::AbstractString)

Convert JTS fixture WKT to a GeometryOps tuple geometry.  Unsupported payloads
return `JTSRawGeometry` so broad fixture loading stays non-throwing.
"""
function jts_wkt_to_geom(wkt::AbstractString)
    sanitized_wkt = replace(strip(wkt), r"\s+" => " ")
    sanitized_wkt = replace(sanitized_wkt, r"\(\s+" => "(")
    sanitized_wkt = replace(sanitized_wkt, r"\s+\)" => ")")
    sanitized_wkt = replace(sanitized_wkt, r",\s+" => ",")
    isempty(sanitized_wkt) && return nothing
    _is_wkt(sanitized_wkt) ||
        return JTSRawGeometry(sanitized_wkt, "Unsupported non-WKT geometry payload.")
    _is_simple_empty_wkt(sanitized_wkt) && return JTSEmptyGeometry(sanitized_wkt)
    geom = GFT.WellKnownText(GFT.Geom(), sanitized_wkt)
    try
        return GO.tuples(geom)
    catch err
        return JTSRawGeometry(sanitized_wkt, sprint(showerror, err))
    end
end

"""
    JTSOperation

One parsed `<op>` expectation, including resolved arguments and conformance status.
"""
struct JTSOperation
    name::String
    argument_refs::Vector{String}
    arguments::Vector{Any}
    expected::Any
    expected_text::String
    attributes::Dict{String,String}
    status::Symbol
    status_reason::String
end

JTSOperation(
    name::AbstractString,
    argument_refs::Vector{String},
    arguments::Vector{Any},
    expected,
    expected_text::AbstractString,
    attributes::Dict{String,String},
) = JTSOperation(
    String(name),
    argument_refs,
    arguments,
    expected,
    String(expected_text),
    attributes,
    :ok,
    "",
)

function Base.show(io::IO, op::JTSOperation)
    status = op.status == :ok ? "" : ", status = $(op.status)"
    print(io, "JTSOperation(name = $(op.name), expected = $(typeof(op.expected))$status)")
end

"""
    JTSCase

One parsed `<case>` with A/B geometries and the operations asserted for them.
"""
struct JTSCase
    description::String
    geom_a::Any
    geom_b::Any
    operations::Vector{JTSOperation}
end

function Base.show(io::IO, case::JTSCase)
    print(io, "JTSCase(description = $(repr(case.description)), operations = $(length(case.operations)))")
end

"""
    JTSTestSet

One parsed JTS XML file, including run metadata and all selected cases.
"""
struct JTSTestSet
    filepath::String
    description::String
    precision_model::Union{Nothing,JTSPrecisionModel}
    cases::Vector{JTSCase}
end

JTSTestSet(filepath::String, cases::Vector{JTSCase}) = JTSTestSet(filepath, "", nothing, cases)

"""
    load_test_cases(filepath; kwargs...)

Load only the parsed case vector from a JTS XML fixture.
"""
load_test_cases(filepath::AbstractString; kwargs...) = load_test_set(filepath; kwargs...).cases

"""
    load_test_set(filepath; operations = nothing, skip = (), broken = (), unimplemented = ())

Parse one JTS XML fixture file into a reusable `JTSTestSet`.
"""
function load_test_set(
    filepath::AbstractString;
    operations = nothing,
    skip = (),
    broken = (),
    unimplemented = (),
)
    doc = read(filepath, XML.Node)
    run = _first_element(doc, "run")
    isnothing(run) && throw(ArgumentError("Expected a <run> root in $filepath."))

    description_node = _first_child_element(run, "desc")
    description = isnothing(description_node) ? "" : String(_node_text(description_node))
    precision_model = _parse_precision_model(_first_child_element(run, "precisionModel"))

    cases = JTSCase[]
    for case_node in _element_children(run, "case")
        case = parse_case(
            case_node,
            String(filepath);
            operations,
            skip,
            broken,
            unimplemented,
        )
        isnothing(case) || push!(cases, case)
    end
    return JTSTestSet(String(filepath), description, precision_model, cases)
end

"""
    load_test_sets(path_or_paths; family = nothing, filename = nothing, kwargs...)

Discover or load multiple JTS XML fixtures and parse each with `load_test_set`.
"""
function load_test_sets(path_or_paths; family = nothing, filename = nothing, kwargs...)
    files = if path_or_paths isa AbstractString
        find_jts_test_files(path_or_paths; family, filename)
    else
        files = String.(path_or_paths)
        filter!(files) do file
            _fixture_family_allowed(file, family) && _filename_allowed(basename(file), filename)
        end
        sort!(files)
        files
    end
    return [load_test_set(path; kwargs...) for path in files]
end

"""
    find_jts_test_files(path; family = nothing, filename = nothing)

Return sorted JTS XML fixture paths from a file or directory tree.
"""
function find_jts_test_files(path::AbstractString; family = nothing, filename = nothing)
    files = if isfile(path)
        [String(path)]
    elseif isdir(path)
        found = String[]
        for (dirpath, _, filenames) in walkdir(path)
            for name in filenames
                endswith(lowercase(name), ".xml") || continue
                push!(found, joinpath(dirpath, name))
            end
        end
        found
    else
        throw(ArgumentError("No JTS XML file or directory found at $path."))
    end

    filter!(files) do file
        _fixture_family_allowed(file, family) && _filename_allowed(basename(file), filename)
    end
    sort!(files)
    return files
end

function parse_case(
    case_node::XML.Node,
    filepath::String;
    operations = nothing,
    skip = (),
    broken = (),
    unimplemented = (),
)
    desc_node = _first_child_element(case_node, "desc")
    a_node = _first_child_element(case_node, "a")
    b_node = _first_child_element(case_node, "b")

    description = isnothing(desc_node) ? "" : String(_node_text(desc_node))
    geom_a = isnothing(a_node) ? nothing : jts_wkt_to_geom(_node_text(a_node))
    geom_b = isnothing(b_node) ? nothing : jts_wkt_to_geom(_node_text(b_node))

    parsed_operations = JTSOperation[]
    for test_node in _element_children(case_node, "test")
        for op_node in _element_children(test_node, "op")
            op = parse_operation(op_node, geom_a, geom_b)
            if _operation_allowed(op.name, operations)
                push!(
                    parsed_operations,
                    _tag_operation(op, filepath, description; skip, broken, unimplemented),
                )
            end
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

"""
    is_runnable(op)

Return true when fixture rules left this operation as a normal test target.
"""
is_runnable(op::JTSOperation) = op.status == :ok

"""
    is_skipped(op)

Return true when a fixture rule marked this operation as intentionally skipped.
"""
is_skipped(op::JTSOperation) = op.status == :skip

"""
    is_broken(op)

Return true when a fixture rule marked this operation as a known failure.
"""
is_broken(op::JTSOperation) = op.status == :broken

"""
    is_unimplemented(op)

Return true when a fixture rule marked this operation as planned future work.
"""
is_unimplemented(op::JTSOperation) = op.status == :unimplemented

"""
    is_overlay_operation(name_or_op)

Identify JTS fixture operations that belong to overlay conformance.
"""
is_overlay_operation(name::AbstractString) = lowercase(name) in OVERLAY_OPERATION_NAMES
is_overlay_operation(op::JTSOperation) = is_overlay_operation(op.name)

"""
    is_relate_operation(name_or_op)

Identify JTS fixture operations that belong to relate/predicate conformance.
"""
is_relate_operation(name::AbstractString) = lowercase(name) in RELATE_OPERATION_NAMES
is_relate_operation(op::JTSOperation) = is_relate_operation(op.name)

"""
    fixture_family(filepath)

Classify a fixture filename as `:relate`, `:overlay`, or `:other`.
"""
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

"""
    geometry_category(geom)

Classify parsed fixture geometry by coarse topological dimension.
"""
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

"""
    case_category(case)

Classify a fixture case from its A/B geometry categories.
"""
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

_parse_optional_float(::Nothing) = nothing
_parse_optional_float(value::AbstractString) = tryparse(Float64, strip(value))

_parse_precision_model(::Nothing) = nothing

function _parse_precision_model(node::XML.Node)
    attrs = Dict(String(k) => String(v) for (k, v) in XML.attributes(node))
    return JTSPrecisionModel(attrs)
end

function _tag_operation(
    op::JTSOperation,
    filepath::String,
    case_description::String;
    skip = (),
    broken = (),
    unimplemented = (),
)
    for (status, rules) in (
        (:skip, skip),
        (:unimplemented, unimplemented),
        (:broken, broken),
    )
        reason = _matching_rule_reason(rules, filepath, case_description, op.name)
        isnothing(reason) || return JTSOperation(
            op.name,
            op.argument_refs,
            op.arguments,
            op.expected,
            op.expected_text,
            op.attributes,
            status,
            reason,
        )
    end
    return op
end

function _matching_rule_reason(rules, filepath::String, case_description::String, op_name::String)
    for rule in _as_rules(rules)
        _rule_matches(rule, filepath, case_description, op_name) && return rule.reason
    end
    return nothing
end

_as_rules(::Nothing) = ()
_as_rules(rule::JTSFixtureRule) = (rule,)
_as_rules(rules) = rules

function _rule_matches(
    rule::JTSFixtureRule,
    filepath::String,
    case_description::String,
    op_name::String,
)
    return _matcher_matches(rule.file, filepath) &&
           _matcher_matches(rule.case, case_description) &&
           _matcher_matches(rule.operation, op_name; exact_string = true)
end

_matcher_matches(::Nothing, value; exact_string = false) = true

function _matcher_matches(matcher::Function, value; exact_string = false)
    return matcher(value)
end

function _matcher_matches(matcher::Regex, value; exact_string = false)
    return occursin(matcher, String(value))
end

function _matcher_matches(matcher::Union{AbstractString,Symbol}, value; exact_string = false)
    matcher_text = lowercase(String(matcher))
    value_text = lowercase(String(value))
    exact_string && return matcher_text == value_text
    return occursin(matcher_text, value_text)
end

function _matcher_matches(matchers, value; exact_string = false)
    return any(matcher -> _matcher_matches(matcher, value; exact_string), matchers)
end

function _fixture_family_allowed(filepath::AbstractString, family::Nothing)
    return true
end

function _fixture_family_allowed(filepath::AbstractString, family)
    return _matcher_matches(family, fixture_family(filepath); exact_string = true)
end

function _filename_allowed(filename::AbstractString, pattern::Nothing)
    return true
end

function _filename_allowed(filename::AbstractString, pattern)
    return _matcher_matches(pattern, filename)
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

function _first_child_element(node::XML.Node, tag_name::AbstractString)
    for child in children(node)
        tag(child) == tag_name && return child
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
