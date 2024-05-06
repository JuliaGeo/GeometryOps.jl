import{_ as e,c as a,o as t,a6 as o}from"./chunks/framework.xqf8KCwY.js";const f=JSON.parse('{"title":"Introduction","description":"","frontmatter":{},"headers":[],"relativePath":"introduction.md","filePath":"introduction.md","lastUpdated":null}'),i={name:"introduction.md"},r=o('<h1 id="Introduction" tabindex="-1">Introduction <a class="header-anchor" href="#Introduction" aria-label="Permalink to &quot;Introduction {#Introduction}&quot;">​</a></h1><p>GeometryOps.jl is a package for geometric calculations on (primarily 2D) geometries.</p><p>The driving idea behind this package is to unify all the disparate packages for geometric calculations in Julia, and make them <a href="https://github.com/JuliaGeo/GeoInterface.jl" target="_blank" rel="noreferrer">GeoInterface.jl</a>-compatible. We seem to be focusing primarily on 2/2.5D geometries for now.</p><p>Most of the usecases are driven by GIS and similar Earth data workflows, so this might be a bit specialized towards that, but methods should always be general to any coordinate space.</p><p>We welcome contributions, either as pull requests or discussion on issues!</p><h2 id="Main-concepts" tabindex="-1">Main concepts <a class="header-anchor" href="#Main-concepts" aria-label="Permalink to &quot;Main concepts {#Main-concepts}&quot;">​</a></h2><h3 id="The-apply-paradigm" tabindex="-1">The <code>apply</code> paradigm <a class="header-anchor" href="#The-apply-paradigm" aria-label="Permalink to &quot;The `apply` paradigm {#The-apply-paradigm}&quot;">​</a></h3><div class="tip custom-block"><p class="custom-block-title">Note</p><p>See the <a href="/GeometryOps.jl/previews/PR135/source/primitives#Primitive-functions">Primitive Functions</a> page for more information on this.</p></div><p>The <code>apply</code> function allows you to decompose a given collection of geometries down to a certain level, and then operate on it.</p><p>Functionally, it&#39;s similar to <code>map</code> in the way you apply it to geometries.</p><p><code>apply</code> and <code>applyreduce</code> take any geometry, vector of geometries, collection of geometries, or table (like <code>Shapefile.Table</code>, <code>DataFrame</code>, or <code>GeoTable</code>)!</p><h3 id="What&#39;s-this-GeoInterface.Wrapper-thing?" tabindex="-1">What&#39;s this <code>GeoInterface.Wrapper</code> thing? <a class="header-anchor" href="#What&#39;s-this-GeoInterface.Wrapper-thing?" aria-label="Permalink to &quot;What&amp;#39;s this `GeoInterface.Wrapper` thing? {#What&#39;s-this-GeoInterface.Wrapper-thing?}&quot;">​</a></h3><p>Write a comment about GeoInterface.Wrapper and why it helps in type stability to guarantee a particular return type.</p>',13),n=[r];function c(s,p,l,d,h,m){return t(),a("div",null,n)}const g=e(i,[["render",c]]);export{f as __pageData,g as default};
