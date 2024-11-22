import GeometryOps: _build_ab_list, PolyNode, _diff_delay_cross_f, _diff_delay_bounce_f, _diff_step, _get_poly_type

function _trace_polynodes_for_mutual(::Type{T}, a_list, b_list, a_idx_list, f_step, poly_a, poly_b) where T
    n_a_pts, n_b_pts = length(a_list), length(b_list)
    total_pts = n_a_pts + n_b_pts 
    n_cross_pts = length(a_idx_list)
    return_polys = Vector{_get_poly_type(T)}(undef, 0)
    
    # Keep track of processed intersection points
    visited_pts = 0
    processed_pts = 0
    first_idx = 1
    @info "Starting main loop to process $(n_cross_pts) crossing points"
    while processed_pts < n_cross_pts
        @info "Processing crossing point $(processed_pts + 1) of $(n_cross_pts)"
        curr_list, curr_npoints = a_list, n_a_pts
        on_a_list = true
        
        # Find first unprocessed intersecting point
        visited_pts += 1
        processed_pts += 1
        first_idx = findnext(x -> x != 0, a_idx_list, first_idx)
        idx = a_idx_list[first_idx]
        a_idx_list[first_idx] = 0
        start_pt = a_list[idx]

        # Initialize polygon trace
        curr = curr_list[idx]
        pt_list = [curr.point]
        @info "Starting point" curr.point
        
        # For each entry point, find matching exit point
        if curr.ent_exit # If entry point
            @info "Found entry point at index $(idx)" point=curr.point
            # Find next exit point by scanning forward
            next_idx = idx
            while true
                next_idx = next_idx < curr_npoints ? next_idx + 1 : 1
                next_pt = curr_list[next_idx]
                
                # Add each point to trace
                push!(pt_list, next_pt.point)
                
                if next_pt.crossing && !next_pt.ent_exit
                    @info "Found matching exit point at index $(next_idx)" point=next_pt.point
                    # Found exit point - mark it and break
                    processed_pts += 1
                    a_idx_list[next_pt.idx] = 0
                    break
                end
                visited_pts += 1
                @assert visited_pts < total_pts "Trace error"
            end
            
            @info "Completing polygon by adding remaining points back to start" current_point=curr_list[next_idx].point target_point=curr.point
            # Continue adding remaining points back to start
            curr_idx = next_idx
            while true
                curr_idx = curr_idx < curr_npoints ? curr_idx + 1 : 1
                if curr_idx == idx
                    break
                end
                push!(pt_list, curr_list[curr_idx].point)
                visited_pts += 1
                @assert visited_pts < total_pts "Trace error"
            end
        end

        @info "Created polygon with $(length(pt_list)) points" points=pt_list
        
        push!(return_polys, GI.Polygon([pt_list]))
    end
    #= To trace both polygons:
    1. Need to track entry/exit points for both polygons simultaneously
    2. For each entry point in either polygon:
       - Find matching exit point in same polygon
       - Add segment between entry and exit
       - Mark both points as processed
    3. Could potentially alternate between polygons to build medial axis
    4. Would need careful handling of overlapping regions
    =#
    
    @info "Completed tracing with $(length(return_polys)) polygons found" polygons=return_polys
    return return_polys
end

function _trace_polynodes_for_mutual_2(::Type{T}, a_list, b_list, a_idx_list, f_step, poly_a, poly_b) where T
    n_a_pts, n_b_pts = length(a_list), length(b_list)
    total_pts = n_a_pts + n_b_pts 
    n_cross_pts = length(a_idx_list)
    return_polys = Vector{_get_poly_type(T)}(undef, 0)

    # Keep track of processed intersection points
    visited_pts = 0
    processed_pts = 0
    first_idx = 1

    if n_cross_pts == 0
        return [node.point for node in a_list]
    end

    @info "Starting main loop to process $(n_cross_pts) crossing points"
    # while processed_pts < n_cross_pts
        @info "Processing crossing point $(processed_pts + 1) of $(n_cross_pts)"
        curr_list, curr_npoints = a_list, n_a_pts
        on_a_list = true

        skipping = false
        current_ent_exit = false
        is_first_point = true

        pt_list = Tuple{T, T}[]

        first_intersection_point = findfirst(x -> x.crossing, a_list)
        
        for point in view(a_list, vcat(first_intersection_point:length(a_list), 1:first_intersection_point - 1))
            if point.crossing
                if is_first_point
                    is_first_point = false
                    current_ent_exit = point.ent_exit
                end

                if point.ent_exit == GO.exit
                    skipping = false
                else # entry point
                    skipping = true
                end
                push!(pt_list, point.point)
                current_ent_exit = point.ent_exit
            elseif !skipping
                push!(pt_list, point.point)
            end
        end

    # end

    if pt_list[1] != pt_list[end]
        push!(pt_list, pt_list[1])
    end
    return pt_list
end
    

function mutual_boundary_adjustment(::GI.PolygonTrait, a, ::GI.PolygonTrait, b; exact = true)
    T = Float64
    
    if GI.nhole(a) > 0 || GI.nhole(b) > 0
        throw(ArgumentError("""
        Holes are not supported in `mutual_boundary_adjustment` yet.  
        Please file an issue at JuliaGeo/GeometryOps.jl if you need this.

        Found $(GI.nhole(a)) holes in polygon a and $(GI.nhole(b)) holes in polygon b.
        """))
    end
    # get exterior rings of both polygons
    aext, bext = GI.getexterior(a), GI.getexterior(b)
    # build list of points and indices of intersection points
    a_list, b_list, a_idx_list = _build_ab_list(T, aext, bext, _diff_delay_cross_f, _diff_delay_bounce_f; exact)
    # trace polynodes for mutual boundary adjustment
    pt_list = _trace_polynodes_for_mutual_2(T, a_list, b_list, a_idx_list, _diff_step, a, b)
    
    return GI.Polygon([GI.LinearRing(pt_list)])
end


@testset "Overlapping rectangles" begin
    # Create two rectangles that overlap at corner by 1/6 area
    # First rectangle: 6x6 at origin
    rect1 = GI.Polygon([GI.LinearRing([(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (0.0, 6.0), (0.0, 0.0)])])

    # Second rectangle: 6x6 with 2x3 overlap area
    # Positioned so bottom left corner is at (4,3)
    rect2 = GI.Polygon([GI.LinearRing([(4.0, 3.0), (10.0, 3.0), (10.0, 9.0), (4.0, 9.0), (4.0, 3.0)])])

    expected_rect1 = GI.Polygon([[(6.0, 3.0), (4.0, 6.0), (0.0, 6.0), (0.0, 0.0), (6.0, 0.0), (6.0, 3.0)]])
    expected_rect2 = GI.Polygon([[(6.0, 3.0), (10.0, 3.0), (10.0, 9.0), (4.0, 9.0), (4.0, 6.0), (6.0, 3.0)]])

    received_rect1 = mutual_boundary_adjustment(rect1, rect2)
    received_rect2 = mutual_boundary_adjustment(rect2, rect1)

    @test GO.equals(received_rect1, expected_rect1)
    @test GO.equals(received_rect2, expected_rect2)

end

@testset "Overlapping circles" begin
    p1 = GI.Point(0.0, 0.0)
    p2 = GI.Point(1.0, 0.0)

    c1 = GO.buffer(p1, 0.6)
    c2 = GO.buffer(p2, 0.6)

    received_c1 = mutual_boundary_adjustment(c1, c2)
    received_c2 = mutual_boundary_adjustment(c2, c1)

    @test GO.equals(received_c1, expected_c1)
end

@testset "Edge case: No intersection points" begin     
    # Create two non-overlapping rectangles     
    rect1 = GI.Polygon([GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)])])     
    rect2 = GI.Polygon([GI.LinearRing([(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 3.0), (2.0, 2.0)])])      
    # When there are no intersection points, the function should return the original polygon     
    received_rect1 = mutual_boundary_adjustment(rect1, rect2)     
    received_rect2 = mutual_boundary_adjustment(rect2, rect1)      
    @test GO.equals(received_rect1, rect1)     
    @test GO.equals(received_rect2, rect2) 
end

@testset "Edge case: Multiple intersection points in wrong order" begin
    # Create a concave polygon that intersects with a rectangle in multiple places
    # The intersection points will not be in sequential order along the polygon boundary
    concave = GI.Polygon([GI.LinearRing([
        (0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (3.0, 4.0), 
        (3.0, 1.0), (1.0, 1.0), (1.0, 4.0), (0.0, 4.0), (0.0, 0.0)
    ])])
    
    rect = GI.Polygon([GI.LinearRing([
        (2.0, 0.5), (5.0, 0.5), (5.0, 3.5), (2.0, 3.5), (2.0, 0.5)
    ])])
    
    # This case should break the current implementation because it assumes
    # intersection points will be found in sequential order
    @test_throws AssertionError mutual_boundary_adjustment(concave, rect)
    @test_throws AssertionError mutual_boundary_adjustment(rect, concave)
end


using GeoJSON
alex_fc = GeoJSON.read("/Users/singhvi/Downloads/test.geojson")

mutual_boundary_adjustment(alex_fc.geometry[1], alex_fc.geometry[2])