module GeometryOpsTestHelpers

import GeometryOps as GO
using GeoInterface
using Test

export @test_implementations, @testset_implementations

# List of test modules - will be populated when extensions load
# GeoInterface is always available as it's a regular dependency
const TEST_MODULES = Module[GeoInterface]

"""
    conversion_expr(mod, var, genkey)

Generate an expression to convert a geometry variable to a specific module's type.
Handles special cases for Extents and empty geometries.
"""
function conversion_expr(mod, var, genkey)
    quote
        $genkey = if $var isa $(GeoInterface.Extents.Extent)
            # GeoInterface and LibGEOS support Extents directly
            if string(nameof($mod)) in ("GeoInterface", "LibGEOS")
                $var
            else
                $GeoInterface.convert($mod, $(GO.extent_to_polygon)($var))
            end
        # These modules do not support empty geometries.
        # GDAL does but AG does not
        elseif string(nameof($mod)) in ("GeoInterface", "ArchGDAL", "GeometryBasics") && $GeoInterface.isempty($var)
            $var
        else
            $GeoInterface.convert($mod, $var)
        end
    end
end

"""
    @test_implementations(code)
    @test_implementations(modules, code)

Macro to run a block of `code` for multiple modules, using GeoInterface.convert
for each variable prefixed with `\$` in the code block.

# Examples
```julia
point = GI.Point(1.0, 2.0)
@test_implementations begin
    \$point isa GeoInterface.AbstractGeometry
end
```
"""
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

"""
    @testset_implementations(code)
    @testset_implementations(title, code)
    @testset_implementations(modules, code)
    @testset_implementations(title, modules, code)

Macro to run a block of `code` for multiple modules within separate testsets,
using GeoInterface.convert for each variable prefixed with `\$` in the code block.

# Examples
```julia
point = GI.Point(1.0, 2.0)
@testset_implementations "Point tests" begin
    @test GeoInterface.x(\$point) == 1.0
end
```
"""
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

end # module GeometryOpsTestHelpers
