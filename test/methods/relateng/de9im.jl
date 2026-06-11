using Test
import GeometryOps as GO
import GeometryOps: DE9IM

@testset "codes" begin
    @test GO.LOC_INTERIOR == 0 && GO.LOC_BOUNDARY == 1 && GO.LOC_EXTERIOR == 2
    @test GO.DIM_FALSE == -1 && GO.DIM_P == 0 && GO.DIM_L == 1 && GO.DIM_A == 2
    @test GO.dim_char(GO.DIM_FALSE) == 'F'
    @test GO.dim_char(GO.DIM_A) == '2'
    @test GO.dim_code('T') == GO.DIM_TRUE && GO.dim_code('*') == GO.DIM_DONTCARE
    # JTS Dimension.toDimensionValue upper-cases its input (Dimension.java)
    @test GO.dim_code('t') == GO.DIM_TRUE
    @test GO.dim_code('f') == GO.DIM_FALSE
    @test_throws ArgumentError GO.dim_code('X')
end

@testset "DimensionLocation" begin
    # Constants verbatim from JTS DimensionLocation.java
    @test GO.DL_POINT_INTERIOR == 103
    @test GO.DL_LINE_INTERIOR == 110 && GO.DL_LINE_BOUNDARY == 111
    @test GO.DL_AREA_INTERIOR == 120 && GO.DL_AREA_BOUNDARY == 121
    @test GO.dimloc_location(GO.DL_AREA_BOUNDARY) == GO.LOC_BOUNDARY
    @test GO.dimloc_location(GO.DL_POINT_INTERIOR) == GO.LOC_INTERIOR
    @test GO.dimloc_dimension(GO.DL_LINE_INTERIOR) == GO.DIM_L
    @test GO.dimloc_dimension(GO.DL_EXTERIOR) == GO.DIM_FALSE
    # Two-arg overload (JTS DimensionLocation.dimension(dimLoc, exteriorDim)):
    # exterior returns the passed exterior dimension, others ignore it.
    @test GO.dimloc_dimension(GO.DL_EXTERIOR, GO.DIM_A) == GO.DIM_A
    @test GO.dimloc_dimension(GO.DL_EXTERIOR, GO.DIM_FALSE) == GO.DIM_FALSE
    @test GO.dimloc_dimension(GO.DL_AREA_INTERIOR, GO.DIM_P) == GO.DIM_A
    @test GO.dimloc_dimension(GO.DL_LINE_BOUNDARY, GO.DIM_A) == GO.DIM_L
    @test GO.dimloc_dimension(GO.DL_POINT_INTERIOR, GO.DIM_A) == GO.DIM_P
    @test GO.dimloc_location(GO.DL_EXTERIOR) == GO.LOC_EXTERIOR
    @test GO.dimloc_area(GO.LOC_INTERIOR) == GO.DL_AREA_INTERIOR
    @test GO.dimloc_line(GO.LOC_BOUNDARY) == GO.DL_LINE_BOUNDARY
    @test GO.dimloc_point(GO.LOC_EXTERIOR) == GO.DL_EXTERIOR
end

@testset "DE9IM" begin
    im = DE9IM("212101212")
    @test string(im) == "212101212"
    @test im[GO.LOC_INTERIOR, GO.LOC_INTERIOR] == GO.DIM_A
    @test im[GO.LOC_BOUNDARY, GO.LOC_BOUNDARY] == GO.DIM_P
    @test DE9IM() == DE9IM("FFFFFFFFF")
    im2 = GO.with_entry(DE9IM(), GO.LOC_INTERIOR, GO.LOC_BOUNDARY, GO.DIM_L)
    @test string(im2) == "F1FFFFFFF"
    # pattern matching (JTS IntersectionMatrix.matches semantics)
    @test GO.matches(DE9IM("212101212"), "T*F**FFF*") == false
    @test GO.matches(DE9IM("2FF1FF212"), "T*F**FFF*") == false
    @test GO.matches(DE9IM("2FF1FFFF2"), "T*F**FFF*") == true
    @test GO.matches(DE9IM("0FFFFFFF2"), "0FFFFFFF2") == true
    # 'T' pattern also matches a 'T' matrix entry (JTS IntersectionMatrix.matches)
    @test GO.matches(GO.DE9IM("TFFFFFFFF"), "TFFFFFFFF") == true
    # lowercase pattern codes are accepted (JTS Dimension.toDimensionValue upper-cases)
    @test GO.matches(GO.DE9IM("2FF1FFFF2"), "t*f**fff*") == true
    @test_throws ArgumentError DE9IM("212")        # wrong length
    @test_throws ArgumentError GO.matches(DE9IM(), "T*F**FF")  # wrong length
end

@testset "BoundaryNodeRule" begin
    @test GO.is_in_boundary(GO.Mod2Boundary(), 1) == true
    @test GO.is_in_boundary(GO.Mod2Boundary(), 2) == false
    @test GO.is_in_boundary(GO.Mod2Boundary(), 3) == true
    @test GO.is_in_boundary(GO.EndpointBoundary(), 1) == true
    @test GO.is_in_boundary(GO.EndpointBoundary(), 2) == true
    @test GO.is_in_boundary(GO.MultivalentEndpointBoundary(), 1) == false
    @test GO.is_in_boundary(GO.MultivalentEndpointBoundary(), 2) == true
    @test GO.is_in_boundary(GO.MonovalentEndpointBoundary(), 1) == true
    @test GO.is_in_boundary(GO.MonovalentEndpointBoundary(), 2) == false
end
