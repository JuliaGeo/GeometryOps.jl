using XML

import WellKnownGeometry
import GeoFormatTypes as GFT
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

"""
    jts_wkt_to_geom(wkt::String)

Convert a JTS WKT string to a GeometryOps geometry, via WellKnownGeometry.jl and GO.tuples.

The reason this exists is because WellKnownGeometry doesn't work well with newlines in subsidiary geometries,
so this sanitizes the input before parsing and converting.

Some WKT is parsed via LibGEOS instead, because WellKnownGeometry can't handle it:

- WKT containing `EMPTY` (including nested, e.g. `GEOMETRYCOLLECTION(POLYGON EMPTY, ...)`),
  because GI wrapper geometries cannot be empty.
- `GEOMETRYCOLLECTION`, because WellKnownGeometry mis-splits subgeometries
  preceded by whitespace (e.g. `GEOMETRYCOLLECTION( LINESTRING (1 2, 1 1))`).
- `LINEARRING`, which WellKnownGeometry does not know at all.
"""
function jts_wkt_to_geom(wkt::String)
    sanitized_wkt = join(strip.(split(wkt, "\n")), "")
    upper_wkt = uppercase(lstrip(sanitized_wkt))
    if occursin("EMPTY", upper_wkt) ||
            startswith(upper_wkt, "GEOMETRYCOLLECTION") ||
            startswith(upper_wkt, "LINEARRING")
        return LG.readgeom(sanitized_wkt)
    end
    geom = GFT.WellKnownText(GFT.Geom(), sanitized_wkt)
    return GO.tuples(geom)
end

# Operations whose expected result is a boolean, not a geometry.
const BOOLEAN_OPS = Set(["relate", "intersects", "disjoint", "contains", "within",
    "covers", "coveredby", "crosses", "touches", "overlaps", "equalstopo", "equals"])

function parse_expected(operation, raw::String)
    lowercase(operation) in BOOLEAN_OPS && return parse(Bool, lowercase(strip(raw)))
    return jts_wkt_to_geom(raw)
end

# Note: geometry fields are untyped because `GO.tuples` returns a bare
# coordinate tuple for POINT WKT, not a wrapper geometry.
struct TestItem{T}
    operation::String
    arg1::Any
    arg2::Any                        # `nothing` for unary ops like getboundary
    pattern::Union{Nothing, String}  # DE-9IM pattern from `arg3` (relate ops)
    expected_result::T
end

Base.show(io::IO, ::MIME"text/plain", item::TestItem) = print(io, "TestItem(operation = $(item.operation), expects $(item.expected_result isa Bool ? item.expected_result : GI.trait(item.expected_result)))")
Base.show(io::IO, item::TestItem) = show(io, MIME"text/plain"(), item)

struct Case
    description::String
    geom_a::Any
    geom_b::Any  # `nothing` for cases (e.g. TestBoundary) with no <b>
    items::Vector{TestItem}
end

"""
    Run

The parsed contents of one JTS XML test file: run-level metadata plus all cases.

- `precision_model` is `"FLOATING"` (also the default when no `<precisionModel>`
  element is present) or `"FIXED"` (a `<precisionModel>` with a `scale` attribute).
  Used to skip FIXED-precision cases, which assume snapped coordinates.
- `n_skipped_ops` counts `<op>` elements that could not be represented as a
  `TestItem` (e.g. `arg1`/`arg2` referring to something other than the case's
  A/B geometries); they are dropped from `items` but never silently — the count
  is recorded here.
"""
struct Run
    filepath::String
    description::Union{Nothing, String}
    precision_model::String
    cases::Vector{Case}
    n_skipped_ops::Int
end

function load_test_run(filepath::String)
    doc = read(filepath, XML.Node) # lazy parsing
    run = only(children(doc))
    description = nothing
    precision_model = "FLOATING"
    test_cases = Case[]
    n_skipped_ops = Ref(0)
    for child in children(run)
        t = tag(child)
        if t == "case"
            push!(test_cases, parse_case(child, n_skipped_ops))
        elseif t == "desc"
            description = strip(value(only(children(child))))
        elseif t == "precisionModel"
            pm_attrs = something(XML.attributes(child), Dict{String, String}())
            if haskey(pm_attrs, "type")
                precision_model = uppercase(pm_attrs["type"])
            elseif haskey(pm_attrs, "scale")
                precision_model = "FIXED"
            end
        end
        # other tags (comments, <resultMatcher>, ...) are run-level config we don't use
    end
    return Run(filepath, description, precision_model, test_cases, n_skipped_ops[])
end

"""
    load_test_cases(filepath)

Parse a JTS XML test file and return its `Vector{Case}`.
Convenience wrapper around [`load_test_run`](@ref) for consumers that
don't need the run-level metadata (e.g. `overlay_runner.jl`).
"""
load_test_cases(filepath::String) = load_test_run(filepath).cases

# Extract the text content of an XML element, tolerating missing/empty content.
function _element_text(node::XML.Node)
    isnothing(children(node)) && return ""
    kids = children(node)
    isempty(kids) && return ""
    return join((something(value(kid), "") for kid in kids), "\n")
end

function parse_case(case::XML.Node, n_skipped_ops::Ref{Int} = Ref(0))
    description = ""
    a = nothing
    b = nothing
    test_elements = XML.Node[]
    for child in children(case)
        t = tag(child)
        if t == "desc"
            description = strip(_element_text(child))
        elseif t == "a"
            a = jts_wkt_to_geom(_element_text(child))
        elseif t == "b"
            b = jts_wkt_to_geom(_element_text(child))
        elseif t == "test"
            push!(test_elements, child)
        end
    end
    isnothing(a) && error("case \"$description\" has no <a> geometry")

    # Resolve an argN attribute value to the case's A/B geometry, or `nothing`
    # if it refers to anything else (e.g. a previous op's result).
    resolve_arg(arg) = lowercase(arg) == "a" ? a : (lowercase(arg) == "b" ? b : nothing)

    items = TestItem[]
    for item in test_elements
        for op in children(item)
            tag(op) == "op" || continue
            op_attrs = XML.attributes(op)
            operation = op_attrs["name"]
            arg1 = resolve_arg(op_attrs["arg1"])
            if isnothing(arg1)
                # arg1 refers to something other than A/B; we can't run this op.
                n_skipped_ops[] += 1
                continue
            end
            arg2_raw = get(op_attrs, "arg2", nothing)
            arg2 = isnothing(arg2_raw) ? nothing : resolve_arg(arg2_raw)
            if !isnothing(arg2_raw) && isnothing(arg2)
                n_skipped_ops[] += 1
                continue
            end
            pattern = get(op_attrs, "arg3", nothing)
            expected_result = parse_expected(operation, _element_text(op))
            push!(items, TestItem(operation, arg1, arg2, pattern, expected_result))
        end
    end
    return Case(description, a, b, items)
end
