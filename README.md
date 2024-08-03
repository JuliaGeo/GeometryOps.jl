<img width="400" alt="GeometryOps.jl" src="https://github.com/JuliaGeo/GeometryOps.jl/assets/32143268/92c5526d-23a9-4e01-aee0-2fcea99c5001">

![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGeo.github.io/GeometryOps.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGeo.github.io/GeometryOps.jl/dev/)
[![Build Status](https://github.com/JuliaGeo/GeometryOps.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaGeo/GeometryOps.jl/actions/workflows/CI.yml?query=branch%3Amain)

<img src="docs/src/assets/logo.png" alt="GeometryOps logo" width="250">

> [!WARNING]
> This package is still under heavy development!  Use with care.

GeometryOps.jl is a package for geometric calculations on (primarily 2D) geometries.

The driving idea behind this package is to unify all the disparate packages for geometric calculations in Julia, and make them [GeoInterface.jl](https://github.com/JuliaGeo/GeoInterface.jl)-compatible. We are focusing primarily on 2/2.5D geometries for now.  All methods in this package will consume any geometry which is compatible with GeoInterface - see its [integrations page](https://juliageo.org/GeoInterface.jl/stable/reference/integrations/) for more info on that!

Most of the use cases are driven by GIS and similar Earth data workflows, so this might be a bit specialized towards that, but methods should always be general to any coordinate space.

We welcome contributions, either as pull requests or discussion on issues!

## Methods 

GeometryOps tries to offer most of the basic geometry operations you'd need, implemented in pure Julia and accepting any GeoInterface.jl compatible type.

- General geometry methods (OGC methods): `equals`, `extent`, `distance`, `crosses`, `contains`, `intersects`, etc
- Targeted function application over large nested geometries (`apply`) and reduction over geometries (`applyreduce`)
    - Both `apply` and `applyreduce` consume arbitrary tables as well, like DataFrames!
- `signed_area`, `centroid`, `distance`, etc for valid geometries
- Line and polygon simplification (`simplify`)
- Polygon clipping, `intersection`, `difference` and `union`
- Generalized barycentric coordinates in polygons (`barycentric_coordinates`)
- Projection of geometries between coordinate reference systems using [Proj.jl](https://github.com/JuliaGeo/Proj.jl)
- Polygonization of raster images by contour detection (`polygonize`)
- Segmentization/densification of geometry, both linearly and by geodesic paths (`segmentize`)

See the "API" page in the docs for a more complete list!

## How to navigate the docs

GeometryOps' [docs](https://juliageo.org/GeometryOps.jl/stable) are divided into three main sections: tutorials, explanations and source code.  
Documentation and examples for many functions can be found in the source code section, since we use literate programming in GeometryOps.

- Tutorials are meant to teach the fundamental concepts behind GeometryOps, and how to perform certain operations.
- Explanations usually contain little code, and explain in more detail how GeometryOps works.
- Source code usually contains explanations and examples at the top of the page, followed by annotated source code from that file.

## Performance comparison to other packages

From the wonderful [vector-benchmark](https://www.github.com/kadyb/vector-benchmark),

[![download-3](https://github.com/JuliaGeo/GeometryOps.jl/assets/32143268/0be8672c-c90f-4e1d-81c5-8522317c5e29)](https://github.com/kadyb/vector-benchmark/pull/12)

More benchmarks coming soon!

### Planned additions

- Buffering, hulls (convex and otherwise)
- Checks for valid geometries (empty linestrings, null points, etc) ([#14](https://github.com/JuliaGeo/GeometryOps.jl/issues/14))
- Operations on spherical (non-Euclidean) geometry ([#17](https://github.com/JuliaGeo/GeometryOps.jl/issues/17))
