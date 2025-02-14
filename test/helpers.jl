module TestHelpers

using Test, GeoInterface, ArchGDAL, GeometryBasics, LibGEOS

export @test_implementations, @testset_implementations

const TEST_MODULES = [GeoInterface, ArchGDAL, GeometryBasics, LibGEOS]

# Monkey-patch GeometryBasics to have correct methods.
# TODO: push this up to GB!

@eval GeometryBasics begin
    # MultiGeometry ncoord implementations
    GeoInterface.ncoord(::GeoInterface.MultiPolygonTrait, ::GeometryBasics.MultiPolygon{N}) where N = N
    GeoInterface.ncoord(::GeoInterface.MultiPointTrait, ::GeometryBasics.MultiPoint{N}) where N = N
    # LinearRing and LineString confusion
    GeometryBasics.geointerface_geomtype(::GeoInterface.LinearRingTrait) = LineString
    function GeoInterface.convert(::Type{LineString}, ::GeoInterface.LinearRingTrait, geom)
        return GeoInterface.convert(LineString, GeoInterface.LineStringTrait(), geom) # forward to the linestring conversion method
    end
    # Line interface
    function GeoInterface.convert(::Type{Line}, ::GeoInterface.LineTrait, geom)
        p1, p2 = GeoInterface.getpoint(geom)
        return Line(GeoInterface.convert(Point, GeoInterface.PointTrait(), p1), GeoInterface.convert(Point, GeoInterface.PointTrait(), p2))
    end
    # GeometryCollection interface - currently just a large Union
    const _ALL_GB_GEOM_TYPES = Union{Point, Line, LineString, Polygon, MultiPolygon, MultiLineString, MultiPoint}
    GeometryBasics.geointerface_geomtype(::GeoInterface.GeometryCollectionTrait) = Vector{_ALL_GB_GEOM_TYPES}
    function GeoInterface.convert(::Type{Vector{_ALL_GB_GEOM_TYPES}}, ::GeoInterface.GeometryCollectionTrait, geoms)
        return _ALL_GB_GEOM_TYPES[GeoInterface.convert(GeometryBasics, g) for g in GeoInterface.getgeom(geoms)]
    end
end


# Macro to run a block of `code` for multiple modules,
# using GeoInterface.convert for each var in `args`
macro test_implementations(code::Expr)
    _test_implementations_inner(TEST_MODULES, code)
end
macro test_implementations(modules::Union{Expr,Vector}, code::Expr)
    _test_implementations_inner(modules, code)
end

function _test_implementations_inner(modules::Union{Expr,Vector}, code::Expr)
    vars = Dict{Symbol,Symbol}()
    code1 = _quasiquote!(code, vars)
    modules1 = modules isa Expr ? modules.args : modules
    tests = Expr(:block)

    for mod in modules1
        expr = Expr(:block)
        for (var, genkey) in pairs(vars)
            push!(expr.args, :($genkey = $GeoInterface.convert($mod, $var)))
        end
        push!(expr.args, :(@test $code1))
        push!(tests.args, expr)
    end

    return esc(tests)
end

# Macro to run a block of `code` for multiple modules,
# using GeoInterface.convert for each var in `args`
macro testset_implementations(code::Expr)
    _testset_implementations_inner("", TEST_MODULES, code)
end
macro testset_implementations(arg, code::Expr)
    if arg isa String || arg isa Expr && arg.head == :string
        _testset_implementations_inner(arg, TEST_MODULES, code)
    else
        _testset_implementations_inner("", arg, code)
    end
end
macro testset_implementations(title, modules::Union{Expr,Vector}, code::Expr)
    _testset_implementations_inner(title, modules, code)
end

function _testset_implementations_inner(title, modules::Union{Expr,Vector}, code::Expr)
    vars = Dict{Symbol,Symbol}()
    code1 = _quasiquote!(code, vars)
    modules1 = modules isa Expr ? modules.args : modules
    testsets = Expr(:block)

    for mod in modules1
        expr = Expr(:block)
        for (var, genkey) in pairs(vars)
            push!(expr.args, :($genkey = $GeoInterface.convert($mod, $var)))
        end
        # Manually define the testset macrocall and all string interpolation
        testset = Expr(
            :macrocall, 
            Symbol("@testset"), 
            LineNumberNode(@__LINE__, @__FILE__), 
            Expr(:string, mod, " ", title), 
            code1
        )
        push!(expr.args, testset)
        push!(testsets.args, expr)
    end

    return esc(testsets)
end

# Taken from BenchmarkTools.jl
_quasiquote!(ex, vars) = ex
function _quasiquote!(ex::Expr, vars::Dict)
    if ex.head === :($)
        v = ex.args[1]
        gen = if v isa Symbol 
            haskey(vars, v) ? vars[v] : gensym(v)
        else
            gensym()
        end
        vars[v] = gen
        return v 
    elseif ex.head !== :quote
        for i in 1:length(ex.args)
            ex.args[i] = _quasiquote!(ex.args[i], vars)
        end
    end
    return ex
end

end
