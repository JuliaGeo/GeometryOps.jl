# Paper

This page will follow the JOSS template, but I'm just getting some ideas down that we can use later if required.

## Statement of need

!!! note
    These are basically a bunch of ideas/lines that we can use, it's not fleshed out yet!

Many full-featured geospatial libraries exist already - GEOS (and its predecessor JTS) being the main sources 
from which libraries like GDAL, R's sf, Python's shapely and others obtain functionality.

GeometryOps is able to integrate all of Julia's functionality like arbitrary float types and multiple dispatch.  
This, along with being written in a high-level language like Julia, significantly lowers the barrier to entry
for new contributors.  We are aiming to make impossible workflows (whether because of memory constraints, speed,
or lack of good algorithms) possible, not merely for performance improvement!

A guiding philosophy of GeometryOps is the ability to provide choice.  Consider the `segmentize` function - 
users can choose whether to interpolate in linear (handled directly) or geodetic (handled by Karney's GeographicLib) 
space.  This is not a functionality that is available in most other libraries, where one must instead choose the library
(for example, s2 vs sf) or activate a global switch to toggle certain functionality.

Another example of this philosophy is the `fix` interface.  Users can create their own fixes, relative to their needs, and hook
those into an already known syntax.  Error messages can be customized and users can potentially even have them show plots indicating
exactly where the error is.  No other geometry library (to my knowledge) offers this level of flexibility.

GeometryOps also utilizes exact predicates to return geometrically correct answers, which even GEOS does not. 

## Ongoing research projects
- Subzero.jl (ice floe simulation, OOM better performance than Matlab, 3x better than LibGEOS with more accurate results) using polygon intersection
- Alex Gardner's stuff (glacier tracking, statistics and forecasting) using polygonize and generic spatial predicates/set ops
- Anyone else?

## Citations

- Core Julia packages: julia, ExactPredicates, GeoInterface
- Foster 2019 paper (polygon clipping)
- Hao-Sun paper (point-in-polygon)
- Previous efforts in Julia: PolygonOps.jl, ...
- 
