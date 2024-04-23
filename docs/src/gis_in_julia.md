# GIS in Julia: a brief introduction

This file is an introduction to GIS concepts in Julia.  We will begin by introducing the set of basic 
geometry types in the simple features model and how you can create them in Julia.

Then, we will look at relations between geometries, from the DE-9IM model and spatial predicates, to polygon set operations and distances between geometries.  
We will touch briefly on the nature of Cartesian versus spherical (Haversine) and ellipsoidal (Geodesic) distances, and where each is applicable.

Following this, we will examine how geometries can be manipulated in Julia via the apply function and similar methods.

Finally, we'll examine file input/output and visualization.  We will briefly touch on rasterized data, but that is better left to Rasters.jl and similar 
packages, which we will introduce.
