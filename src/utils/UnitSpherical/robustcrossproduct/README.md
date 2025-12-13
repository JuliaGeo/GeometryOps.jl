# RobustCrossProduct.jl

This is an implementation of robust cross products that will give you an orthogonal 
vector on the unit sphere to two input points, in nearly all cases (including degeneracy and 
antipodal points).

This was adapted from Google's s2 library, licensed under the Apache 2.0 license.

The main entry point is `robust_cross_product(a, b)`, which will return a UnitSphericalPoint.
In general it's about 10x slower than LinearAlgebra.cross if no adjustment is required, 
but it is substantially stabler.