# converts a vector of vector of vector of tuples to vec of vec of vec of array
"""
convert_tuple_to_array(Vector{Vector{Vector{Tuple{::Real}}}})::Vector{Vector{Vector{Vector{::Real}}}}

    This function converts a vector of vectors of vectors of tuples to a 
    vector of vectors of vectors of vectors.
"""
function convert_tuple_to_array(VVVT)
    new_return_obj = Vector{Vector{Vector{Vector{Float64}}}}(undef, 0)
    for vec1 in VVVT
        new_vec1 = Vector{Vector{Vector{Float64}}}(undef, 0)
        for vec2 in vec1
            new_vec2 = Vector{Vector{Float64}}(undef, 0)
            for tup in vec2
                push!(new_vec2, [tup[1], tup[2]])
            end
            push!(new_vec1, new_vec2)
        end
        push!(new_return_obj, new_vec1)
    end
    return new_return_obj
end