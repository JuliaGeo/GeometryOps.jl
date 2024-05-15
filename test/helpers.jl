using Test, GeoInterface, ArchGDAL, GeometryBasics, LibGEOS

const TEST_MODULES = [GeoInterface, ArchGDAL, GeometryBasics, LibGEOS]

# Macro to run a block of `code` for multiple modules, 
# using GeoInterface.convert for each var in `args`
macro test_all_implementations(args, code)
    _test_all_implementations_inner("", args, code)
end
macro test_all_implementations(title::String, args, code)
    _test_all_implementations_inner(string(title, " "), args, code)
end

function _test_all_implementations_inner(title, args, code)
    args1 = esc(args)
    code1 = esc(code)

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
    quote
        for mod in TEST_MODULES 
            @testset "$($(isempty(title) ? "" : "$title : " ))$mod" begin
                $let_expr
            end
        end
    end
end
