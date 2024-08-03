using Test, GeoInterface, ArchGDAL, GeometryBasics, LibGEOS

const TEST_MODULES = [GeoInterface, ArchGDAL, GeometryBasics, LibGEOS]

# Monkey-patch GeometryBasics to have correct methods.
# TODO: push this up to GB!

@eval GeometryBasics begin
    # MultiGeometry ncoord implementations
    GeoInterface.ncoord(::GeoInterface.MultiPolygonTrait, ::GeometryBasics.MultiPolygon{N}) where N = N
    GeoInterface.ncoord(::GeoInterface.MultiLineStringTrait, ::GeometryBasics.MultiLineString{N}) where N = N
    GeoInterface.ncoord(::GeoInterface.MultiPointTrait, ::GeometryBasics.MultiPoint{N}) where N = N
    # LinearRing and LineString confusion
    GeometryBasics.geointerface_geomtype(::GeoInterface.LinearRingTrait) = LineString
    function GeoInterface.convert(::Type{LineString}, ::GeoInterface.LinearRingTrait, geom)
        return GeoInterface.convert(LineString, GeoInterface.LineStringTrait(), geom) # forward to the linestring conversion method
    end
    # Line interface
    GeometryBasics.geointerface_geomtype(::GeoInterface.LineTrait) = Line
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
macro test_implementations(modules, code::Expr)
    _test_implementations_inner(modules, code)
end

function _test_implementations_inner(modules, code)
    vars = Dict{Symbol,Symbol}()
    code1 = esc(_quasiquote!(code, vars))
    modules1 = modules isa Expr ? modules.args : modules
    tests = Expr(:block)
    @show vars

    for mod in modules1[1:end]
        expr = Expr(:block, :(mod = $mod))
        for (var, genkey) in pairs(vars)
            push!(expr.args, :($genkey = GeoInterface.convert(mod, $var)))
        end
        push!(expr.args, code1)
        dump(expr)
        push!(tests.args, expr)
    end

    return testsets
end

# Macro to run a block of `code` for multiple modules,
# using GeoInterface.convert for each var in `args`
macro testset_implementations(code::Expr)
    _testset_implementations_inner("", TEST_MODULES, code)
end
macro testset_implementations(title::String, code::Expr)
    _testset_implementations_inner(title::String, TEST_MODULES, code)
end
macro testset_implementations(modules, code::Expr)
    _testset_implementations_inner("", modules, code)
end
macro testset_implementations(title::String, modules, code::Expr)
    _testset_implementations_inner(title, modules, code)
end

function _testset_implementations_inner(title, modules, code)
    vars = Dict{Symbol,Symbol}()
    code1 = _quasiquote!(code, vars)
    modules1 = modules isa Expr ? modules.args : modules
    testsets = Expr(:block)

    for mod in modules1[1:end]
        expr = Expr(:block, :(mod = $mod))
        for (var, genkey) in pairs(vars)
            push!(expr.args, :($genkey = GeoInterface.convert(mod, $var)))
        end
        push!(expr.args, code1)
        dump(expr)
        push!(testsets.args, expr)
    end

    return testsets
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

macro f(x)
    dump(x)
end
