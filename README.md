## GeometryOps.jl


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://asinghvi17.github.io/GeometryOps.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://asinghvi17.github.io/GeometryOps.jl/dev/)
[![Build Status](https://github.com/asinghvi17/GeometryOps.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/asinghvi17/GeometryOps.jl/actions/workflows/CI.yml?query=branch%3Amain)

<img src="docs/src/assets/logo.png" alt="GeometryOps logo" width="250">

GeometryOps.jl is a package for geometric calculations on (primarily 2D) geometries.

The driving idea behind this package is to unify all the disparate packages for geometric calculations in Julia, and make them GeoInterface.jl-compatible. We seem to be focusing primarily on 2/2.5D geometries for now.

Most of the usecases are driven by GIS and similar Earth data workflows, so this might be a bit specialized towards that, but methods should always be general to any coordinate space.

## Methods 

- Signed area, centroid, distance, etc
- Iteration into geometries (`apply`)
- Line and polygon simplification
- Generalized barycentric coordinates in polygons

### Planned additions

- OGC methods (crosses, contains, intersects, etc)
- Polygon union, intersection and clipping
- Arclength interpolation (absolute and relative)
