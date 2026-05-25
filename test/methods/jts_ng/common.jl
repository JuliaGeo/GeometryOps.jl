using Test
import GeometryOps as GO

@testset "Algorithm markers" begin
    relate = GO.RelateNG()
    overlay = GO.OverlayNG()

    @test GO.manifold(relate) isa GO.Planar
    @test GO.manifold(overlay) isa GO.Planar
    @test GO.rebuild(relate, GO.Planar()) === relate
    @test GO.rebuild(overlay, GO.Planar()) === overlay
    @test relate.boundary_node_rule isa GO.Mod2BoundaryNodeRule
    @test !relate.prepared
    @test !overlay.strict
    @test !overlay.area_result_only
    @test overlay.optimized
    @test overlay.precision_model isa GO.NoPrecisionModel
end

@testset "Topology vocabulary" begin
    @test GO.location_index(GO.loc_interior) == 1
    @test GO.location_index(GO.loc_boundary) == 2
    @test GO.location_index(GO.loc_exterior) == 3

    @test GO.dimension_char(GO.dim_false) == 'F'
    @test GO.dimension_char(GO.dim_point) == '0'
    @test GO.dimension_char(GO.dim_line) == '1'
    @test GO.dimension_char(GO.dim_area) == '2'
    @test GO.dimension_from_char('f') == GO.dim_false
    @test GO.max_dimension(GO.dim_point, GO.dim_area) == GO.dim_area

    @test !GO.is_in_boundary(GO.Mod2BoundaryNodeRule(), 2)
    @test GO.is_in_boundary(GO.Mod2BoundaryNodeRule(), 3)
    @test GO.is_in_boundary(GO.EndpointBoundaryNodeRule(), 1)
    @test GO.is_in_boundary(GO.MultivalentEndpointBoundaryNodeRule(), 2)
    @test GO.is_in_boundary(GO.MonovalentEndpointBoundaryNodeRule(), 1)
end

@testset "IntersectionMatrix" begin
    matrix = GO.IntersectionMatrix()
    @test GO.de9im_string(matrix) == "FFFFFFFFF"
    @test matrix[GO.loc_interior, GO.loc_interior] == GO.dim_false

    matrix[GO.loc_interior, GO.loc_interior] = GO.dim_point
    @test GO.de9im_string(matrix) == "0FFFFFFFF"
    @test GO.matches(matrix, "T********")
    @test !GO.matches(matrix, "F********")

    GO.set_at_least!(matrix, GO.loc_interior, GO.loc_interior, GO.dim_line)
    @test matrix[GO.loc_interior, GO.loc_interior] == GO.dim_line
    GO.set_at_least!(matrix, GO.loc_interior, GO.loc_interior, GO.dim_point)
    @test matrix[GO.loc_interior, GO.loc_interior] == GO.dim_line

    other = GO.IntersectionMatrix("F1FFFFFFF")
    GO.set_at_least!(matrix, other)
    @test matrix[GO.loc_interior, GO.loc_boundary] == GO.dim_line
    @test sprint(show, matrix) == GO.de9im_string(matrix)
end
