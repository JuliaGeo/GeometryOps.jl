using SafeTestsets

@safetestset "DE9IM" begin include("de9im.jl") end
@safetestset "Predicates" begin include("predicates.jl") end
@safetestset "Kernel" begin include("kernel.jl") end
@safetestset "Kernel conformance" begin include("kernel_conformance.jl") end
@safetestset "Spherical kernel conformance" begin include("kernel_conformance_spherical.jl") end
@safetestset "S2 edge-crossing conformance" begin include("s2_crossings_conformance.jl") end
@safetestset "Point locator" begin include("point_locator.jl") end
@safetestset "Indexed point-in-area locator" begin include("indexed_point_in_area.jl") end
@safetestset "RelateGeometry" begin include("relate_geometry.jl") end
@safetestset "Node topology" begin include("node_topology.jl") end
@safetestset "TopologyComputer" begin include("topology_computer.jl") end
@safetestset "Edge intersector" begin include("edge_intersector.jl") end
@safetestset "XML harness" begin include("xml_harness.jl") end
@safetestset "RelateNG engine" begin include("relate_ng.jl") end
@safetestset "JTS XML suite" begin include("xml_suite.jl") end
@safetestset "Spherical relate end-to-end" begin include("spherical_end_to_end.jl") end
@safetestset "LibGEOS differential fuzz" begin include("fuzz.jl") end
@safetestset "Allocations and type stability" begin include("allocations.jl") end
# Further files appended here as tasks land:
# ...
