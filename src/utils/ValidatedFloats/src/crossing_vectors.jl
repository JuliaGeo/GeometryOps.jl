# Products can temporarily exceed the scalar domain; only sums consume them.
# Admitted inputs bound these intermediates below 2^805 and above the EFT grid.
@inline function cross3(a::NTuple{3, CrossingFloat}, b::NTuple{3, CrossingFloat})
    check_domain = Val(false)
    x = admit(add(mul(a[2], b[3], check_domain), -mul(a[3], b[2], check_domain), check_domain))
    y = admit(add(mul(a[3], b[1], check_domain), -mul(a[1], b[3], check_domain), check_domain))
    z = admit(add(mul(a[1], b[2], check_domain), -mul(a[2], b[1], check_domain), check_domain))
    return (x, y, z)
end

@inline function dot3(a::NTuple{3, CrossingFloat}, b::NTuple{3, CrossingFloat})
    check_domain = Val(false)
    x_product = mul(a[1], b[1], check_domain)
    y_product = mul(a[2], b[2], check_domain)
    xy_sum = add(x_product, y_product, check_domain)
    z_product = mul(a[3], b[3], check_domain)
    return admit(add(xy_sum, z_product, check_domain))
end

"Return the bounded three-dimensional cross product as an `SVector`."
@inline function cross(a::StaticVector{3, CrossingFloat},
                       b::StaticVector{3, CrossingFloat})
    return SVector{3, CrossingFloat}(cross3(Tuple(a), Tuple(b)))
end

"Return the bounded three-dimensional dot product."
@inline dot(a::StaticVector{3, CrossingFloat}, b::StaticVector{3, CrossingFloat}) =
    dot3(Tuple(a), Tuple(b))

@inline cross(a::NTuple{3, CrossingFloat}, b::NTuple{3, CrossingFloat}) = cross3(a, b)
@inline dot(a::NTuple{3, CrossingFloat}, b::NTuple{3, CrossingFloat}) = dot3(a, b)
