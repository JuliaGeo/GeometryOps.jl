module TestHelpers

using Test, GeoInterface, ArchGDAL, GeometryBasics, LibGEOS

export @test_implementations, @testset_implementations

"Segmentizes geometries to a max distance of 1% of the max dimension of the bounding box when converting.  Useful to test densified geoms in tests!"
module DensifiedGeometries 

    using GeoInterface
    import GeoInterface as GI

    import GeometryOps as GO

    struct DensifiedWrapperGeometry{T, Trait}
        trait::Trait
        geom::T
    end

    geointerface_geomtype(x) = DensifiedWrapperGeometry

    function GeoInterface.convert(::Type{DensifiedWrapperGeometry}, t::GI.AbstractGeometryTrait, geom)
        ext = GI.extent(geom)
        xrange = ext.X[2] - ext.X[1]
        yrange = ext.Y[2] - ext.Y[1]
        max_range = max(max(xrange, yrange), 1.0)
        densified_geom = GO.segmentize(geom, max_distance = max_range / 100)
        return DensifiedWrapperGeometry(t, densified_geom)
    end
    # you can't densify a point
    # NO, MULTIPOINT DOES NOT COUNT
    GeoInterface.convert(::Type{DensifiedWrapperGeometry}, t::Union{GI.AbstractPointTrait, GI.MultiPointTrait}, geom) = DensifiedWrapperGeometry(t, geom)

    GeoInterface.trait(geom::DensifiedWrapperGeometry) = geom.trait
    GeoInterface.geomtrait(geom::DensifiedWrapperGeometry) = geom.trait
    GeoInterface.isgeometry(geom::DensifiedWrapperGeometry) = geom.trait <: GI.AbstractGeometryTrait

    GeoInterface.ngeom(geom::DensifiedWrapperGeometry) = GeoInterface.ngeom(geom.trait, geom.geom)
    GeoInterface.getgeom(geom::DensifiedWrapperGeometry, i) = GeoInterface.getgeom(geom.trait, geom.geom, i)

    GeoInterface.getexterior(geom::DensifiedWrapperGeometry) = GeoInterface.getexterior(geom.trait, geom.geom)
    GeoInterface.nhole(geom::DensifiedWrapperGeometry) = GeoInterface.nhole(geom.trait, geom.geom)
    GeoInterface.gethole(geom::DensifiedWrapperGeometry, i) = GeoInterface.gethole(geom.trait, geom.geom, i)

    GeoInterface.npoint(geom::DensifiedWrapperGeometry) = GeoInterface.npoint(geom.trait, geom.geom)
    GeoInterface.getpoint(geom::DensifiedWrapperGeometry, i) = GeoInterface.getpoint(geom.trait, geom.geom, i)

    GeoInterface.ncoord(geom::DensifiedWrapperGeometry) = GeoInterface.ncoord(geom.trait, geom.geom)
    GeoInterface.getcoord(geom::DensifiedWrapperGeometry, i) = GeoInterface.getcoord(geom.trait, geom.geom, i)
    GeoInterface.coordnames(geom::DensifiedWrapperGeometry) = GeoInterface.coordnames(geom.trait, geom.geom)

end # module DensifiedGeometries

const TEST_MODULES = [GeoInterface, ArchGDAL, GeometryBasics, LibGEOS, DensifiedGeometries]

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
        # This instantiates a ContextTestSet that holds the contents of the let block as "context"
        # and displays it as:
        # ```
        # Expression: compare_GO_LG_clipping(GO_f, LG_f, p1, p2)
        # Context: Testing geometry from module = GeoInterface
        # ```
        # This is a bit of a hack to get around the 90000 different testsets you see if we did this
        # per testset.
        # Ideally, we'd have a custom testset that can (a) display whether the test failed for _every_ geom
        # or just some, and (b) display the context in a more readable way.
        # But for now this is good enough.
        testset = Expr(
            :macrocall, 
            Symbol("@testset"), 
            LineNumberNode(@__LINE__, @__FILE__), 
            Expr(:let, :(var"Testing geometry from module"=$(mod)), code1),
        )
        push!(expr.args, testset)
        push!(testsets.args, expr)
    end

    # Construct a toplevel testset that displays the title and contains all the context testsets
    toplevel = Expr(
        :macrocall,
        Symbol("@testset"),
        LineNumberNode(@__LINE__, @__FILE__),
        Expr(:string, title),
        testsets
    )

    return esc(toplevel)
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
