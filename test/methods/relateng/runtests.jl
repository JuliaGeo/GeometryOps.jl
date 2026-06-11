using SafeTestsets

@safetestset "DE9IM" begin include("de9im.jl") end
@safetestset "Predicates" begin include("predicates.jl") end
@safetestset "Kernel" begin include("kernel.jl") end
@safetestset "Kernel conformance" begin include("kernel_conformance.jl") end
@safetestset "Point locator" begin include("point_locator.jl") end
# Further files appended here as tasks land:
# ...
