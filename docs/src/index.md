```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "GeometryOps.jl"
  text: ""
  tagline: Blazing fast geometry operations in pure Julia
  image:
    src: /logo.png
    alt: GeometryOps
  actions:
    - theme: brand
      text: Introduction
      link: /introduction
    - theme: alt
      text: View on Github
      link: https://github.com/JuliaGeo/GeometryOps.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - icon: <img width="64" height="64" src="https://rawcdn.githack.com/JuliaLang/julia-logo-graphics/f3a09eb033b653970c5b8412e7755e3c7d78db9e/images/juliadots.iconset/icon_512x512.png" alt="Julia code"/>
    title: Pure Julia code
    details: Fast, understandable, extensible functions
    link: /introduction
  - icon: <img width="64" height="64" src="https://fredrikekre.github.io/Literate.jl/v2/assets/logo.png" />
    title: Literate programming
    details: Documented source code with examples!
    link: /source/methods/clipping/cut
  - icon: <img width="64" height="64" src="https://rawcdn.githack.com/JuliaGeo/juliageo.github.io/4788480c2a5f7ae36df67a4b142e3a963024ac91/img/juliageo.svg" />
    title: Full integration with GeoInterface
    details: Use any GeoInterface.jl-compatible geometry
    link: https://juliageo.org/GeoInterface.jl/stable
---


<p style="margin-bottom:2cm"></p>

<div class="vp-doc" style="width:80%; margin:auto">

<h1> What is GeometryOps.jl? </h1>

GeometryOps.jl is a package for geometric calculations on (primarily 2D) geometries.

The driving idea behind this package is to unify all the disparate packages for geometric calculations in Julia, and make them [GeoInterface.jl](https://github.com/JuliaGeo/GeoInterface.jl)-compatible. We seem to be focusing primarily on 2/2.5D geometries for now.

Most of the usecases are driven by GIS and similar Earth data workflows, so this might be a bit specialized towards that, but methods should always be general to any coordinate space.

We welcome contributions, either as pull requests or discussion on issues!


</div>

```

