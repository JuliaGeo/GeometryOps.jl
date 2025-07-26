module TestHelpers

import GeometryOps as GO

using Test, GeoInterface, ArchGDAL, GeometryBasics, LibGEOS

export @test_implementations, @testset_implementations

const TEST_MODULES = [GeoInterface, ArchGDAL, GeometryBasics, LibGEOS]

function conversion_expr(mod, var, genkey)
    quote
        $genkey = if $var isa $(GeoInterface.Extents.Extent)
            if $mod in ($GeoInterface, $LibGEOS)
                $var
            else
                $GeoInterface.convert($mod, $(GO.extent_to_polygon)($var))
            end
        # GeometryBasics does not have a Line geometry type.
        # elseif $mod in ($GeometryBasics,) && $GeoInterface.trait($var) isa $GeoInterface.LineTrait
        #     $var
        elseif $mod in ($GeoInterface, $ArchGDAL, $GeometryBasics) && $GeoInterface.isempty($var)
            $var
        else
            $GeoInterface.convert($mod, $var)
        end
    end
end
# Monkey-patch GeometryBasics to have correct methods.
# TODO: push this up to GB!

    # TODO: remove when GB GI pr lands
    @static if hasmethod(GeometryBasics.convert, (Type{GeometryBasics.LineString}, GeoInterface.LinearRingTrait, Any))
        function GeoInterface.convert(
            ::Type{GeometryBasics.LineString}, 
            ::GeoInterface.LinearRingTrait, 
            geom
            )
            return GeoInterface.convert(GeometryBasics.LineString, GeoInterface.LineStringTrait(), geom)
        end
        GeometryBasics.geointerface_geomtype(::GeoInterface.LinearRingTrait) = GeometryBasics.LineString
    end
    # end todo
    # GeometryCollection interface - currently just a large Union
    const _ALL_GB_GEOM_TYPES = Union{GeometryBasics.Point, GeometryBasics.LineString, GeometryBasics.Polygon, GeometryBasics.MultiPolygon, GeometryBasics.MultiLineString, GeometryBasics.MultiPoint}
    GeometryBasics.geointerface_geomtype(::GeoInterface.GeometryCollectionTrait) = Vector{_ALL_GB_GEOM_TYPES}
    function GeoInterface.convert(::Type{Vector{<: _ALL_GB_GEOM_TYPES}}, ::GeoInterface.GeometryCollectionTrait, geoms)
        return _ALL_GB_GEOM_TYPES[GeoInterface.convert(GeometryBasics, g) for g in GeoInterface.getgeom(geoms)]
    end

    function GeoInterface.convert(
        ::Type{GeometryBasics.LineString}, 
        type::GeoInterface.LineStringTrait, 
        geom::GeoInterface.Wrappers.LinearRing{false, false, GO.StaticArrays.SVector{N, Tuple{Float64, Float64}}, Nothing, Nothing} where N
        )
        return GeometryBasics.LineString(GeometryBasics.Point2{Float64}.(collect(geom.geom)))
    end


    @eval ArchGDAL begin
        function GeoInterface.convert(
            ::Type{T},
            type::GeoInterface.PolygonTrait,
            geom,
        ) where {T<:IGeometry}
            f = get(lookup_method, typeof(type), nothing)
            isnothing(f) && error(
                "Cannot convert an object of $(typeof(geom)) with the $(typeof(type)) trait (yet). Please report an issue.",
            )
            poly = createpolygon()
            foreach(GeoInterface.getring(geom)) do ring
                xs = GeoInterface.x.(GeoInterface.getpoint(ring)) |> collect
                ys = GeoInterface.y.(GeoInterface.getpoint(ring)) |> collect
                subgeom = unsafe_createlinearring(xs, ys)
                result = GDAL.ogr_g_addgeometrydirectly(poly, subgeom)
                @ogrerr result "Failed to add linearring."
            end
            return poly
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
            push!(expr.args, conversion_expr(mod, var, genkey))
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
            push!(expr.args, conversion_expr(mod, var, genkey))
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

end # module
