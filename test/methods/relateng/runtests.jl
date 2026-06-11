using SafeTestsets

@safetestset "DE9IM" begin include("de9im.jl") end
@safetestset "Predicates" begin include("predicates.jl") end
@safetestset "Kernel" begin include("kernel.jl") end
@safetestset "Kernel conformance" begin include("kernel_conformance.jl") end
@safetestset "Point locator" begin include("point_locator.jl") end
@safetestset "RelateGeometry" begin include("relate_geometry.jl") end
@safetestset "Node topology" begin include("node_topology.jl") end
@safetestset "TopologyComputer" begin include("topology_computer.jl") end
@safetestset "XML harness" begin include("xml_harness.jl") end
# Further files appended here as tasks land:
# ...
