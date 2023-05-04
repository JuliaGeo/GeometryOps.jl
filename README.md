# GeometryOps

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://asinghvi17.github.io/GeometryOps.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://asinghvi17.github.io/GeometryOps.jl/dev/)
[![Build Status](https://github.com/asinghvi17/GeometryOps.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/asinghvi17/GeometryOps.jl/actions/workflows/CI.yml?query=branch%3Amain)

GeometryOps.jl is a package which aims to frankenstein mishmash all the methods we need from all the existing geometry packages, then make them compatible with the GeoInterface.

Currently, `src/methods/signed_area.jl` is the only method I would call "complete". 

## Contributing

Contributions are welcome!  We're trying to write this package using literate programming, so you should add lots of comments :D and check out the Literate.jl package for how the syntax works!  

Look at `methods/signed_area.jl` as mentioned above to get an idea of how a function should be defined generically and how to use Literate to the best effect!