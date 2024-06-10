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

@eval LibGEOS begin
    function GI.convert(
        ::Type{GeometryCollection},
        ::GeometryCollectionTrait,
        geom;
        context = get_global_context(),
    )
        return GeometryCollection(GI.convert.((LibGEOS,), GI.getgeom(geom)))
    end
end

# Macro to run a block of `code` for multiple modules, 
# using GeoInterface.convert for each var in `args`
macro test_all_implementations(args, code::Expr)
    _test_all_implementations_inner("", args, TEST_MODULES, code)
end
macro test_all_implementations(title::String, args, code::Expr)
    _test_all_implementations_inner(title::String, args, TEST_MODULES, code)
end
macro test_all_implementations(args, modules, code::Expr)
    _test_all_implementations_inner("", args, modules, code)
end
macro test_all_implementations(title::String, args, modules, code::Expr)
    _test_all_implementations_inner(title, args, modules, code)
end

function _test_all_implementations_inner(title, args, modules, code)
    args1 = esc(args)
    code1 = esc(code)
    modules1 = modules isa Expr ? modules.args : modules

    let_expr = if args isa Symbol # Handle a single variable name
        quote
            let $args1 = GeoInterface.convert(mod, $args1)
                $code1
            end
        end
    else # Handle a tuple of variable names
        quote
            let ($args1 = map(g -> GeoInterface.convert(mod, g), $args1))
                $code1
            end
        end
    end
    testsets = Expr(:block)

    for mod in modules1
        expr = quote
            mod = $mod
            @testset "$mod" begin
                $let_expr
            end
        end
        push!(testsets.args, expr)
    end

    quote 
        @testset "$($title)" begin
            $testsets
        end
    end
end
