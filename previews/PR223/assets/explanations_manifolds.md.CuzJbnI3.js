import{_ as a,c as o,a5 as t,o as i}from"./chunks/framework.C4-WeekU.js";const f=JSON.parse('{"title":"Manifolds","description":"","frontmatter":{},"headers":[],"relativePath":"explanations/manifolds.md","filePath":"explanations/manifolds.md","lastUpdated":null}'),s={name:"explanations/manifolds.md"};function r(n,e,d,l,h,c){return i(),o("div",null,e[0]||(e[0]=[t('<h1 id="manifolds" tabindex="-1">Manifolds <a class="header-anchor" href="#manifolds" aria-label="Permalink to &quot;Manifolds&quot;">​</a></h1><p>A manifold is, mathematically, a description of some space that is locally Euclidean (i.e., locally flat). All geographic projections, and the surface of the sphere and ellipsoid, fall under this category of space - and these are all the spaces that are relevant to geographic geometry.</p><h2 id="What-manifolds-are-available?" tabindex="-1">What manifolds are available? <a class="header-anchor" href="#What-manifolds-are-available?" aria-label="Permalink to &quot;What manifolds are available? {#What-manifolds-are-available?}&quot;">​</a></h2><p>GeometryOps has three <a href="/GeometryOps.jl/previews/PR223/source/src/types#Manifold"><code>Manifold</code></a> types: <a href="./@ref"><code>Planar</code></a>, <a href="./@ref"><code>Spherical</code></a>, and <a href="./@ref"><code>Geodesic</code></a>.</p><ul><li><p><code>Planar()</code> is, as the name suggests, a perfectly Cartesian, usually 2-dimensional, space. The shortest path from one point to another is a straight line.</p></li><li><p><code>Spherical(; radius)</code> describes points on the surface of a sphere of a given radius. The most convenient sphere for geometry processing is the unit sphere, but one can also use the sphere of the Earth for e.g. projections.</p></li><li><p><code>Geodesic(; semimajor_axis, inv_flattening)</code> describes points on the surface of a flattened ellipsoid, similar to the Earth. The parameters describe the curvature and shape of the ellipsoid, and are equivalent to the flags <code>+a</code> and <code>+f</code> in Proj&#39;s ellipsoid specification. The default values are the values of the WGS84 ellipsoid. For <code>Geodesic</code>, we need an <code>AbstractGeodesic</code> that can wrap representations from Proj.jl and SphericalGeodesics.jl.</p></li></ul><p>The idea here is that the manifold describes how the geometry needs to be treated.</p><h2 id="Why-this-is-needed" tabindex="-1">Why this is needed <a class="header-anchor" href="#Why-this-is-needed" aria-label="Permalink to &quot;Why this is needed {#Why-this-is-needed}&quot;">​</a></h2><p>The classical problem this is intended to solve is that in GIS, latitude and longitude coordinates are often treated as planar coordinates, when they in fact live on the sphere/ellipsoid, and must be treated as such. For example, computing the area of the USA on the lat/long plane yields a result of <code>1116</code>, which is plainly nonsensical.</p><h2 id="How-this-is-done" tabindex="-1">How this is done <a class="header-anchor" href="#How-this-is-done" aria-label="Permalink to &quot;How this is done {#How-this-is-done}&quot;">​</a></h2><p>In order to avoid this, we&#39;ve introduced three complementary CRS-related systems to the JuliaGeo ecosystem.</p><ol><li><p>GeoInterface&#39;s <code>crstrait</code>. This is a method that returns the ideal CRS <em>type</em> of a geometry, either Cartesian or Geographic.</p></li><li><p>Proj&#39;s <code>PreparedCRS</code> type, which extracts ellipsoid parameters and the nature of the projection from a coordinate reference system, and caches the results in a struct. This allows GeometryOps to quickly determine the correct manifold to use for a given geometry.</p></li><li><p>GeometryOps&#39;s <code>Manifold</code> type, which defines the surface on which to perform operations. This is what allows GeometryOps to perform calculations correctly depending on the nature of the geometry.</p></li></ol><p>The way this flow works, is that when you load a geometry using GeoDataFrames, its CRS is extracted and parsed into a <code>PreparedCRS</code> type. This is then used to determine the manifold to use for the geometry, and the geometry is converted to the manifold&#39;s coordinate system.</p><p>There is a table of known geographic coordinate systems in GeoFormatTypes.jl, and anything else is assumed to be a Cartesian or planar coordinate system. CRStrait is used as the cheap determinant, but PreparedCRS is more general and better to use if possible.</p><p>When GeometryOps sees a geometry, it first checks its CRS to see if it is a geographic coordinate system. If it is, it uses the <code>PreparedCRS</code>, or falls back to <code>crstrait</code> and geographic defaults to determine the manifold.</p><h2 id="Algorithms-and-manifolds" tabindex="-1">Algorithms and manifolds <a class="header-anchor" href="#Algorithms-and-manifolds" aria-label="Permalink to &quot;Algorithms and manifolds {#Algorithms-and-manifolds}&quot;">​</a></h2><p>Algorithms define what operation is performed on the geometry; however, the choice of algorithm can also depend on the manifold. L&#39;Huilier&#39;s algorithm for the area of a polygon is not applicable to the plane, but is applicable to either the sphere or ellipsoid, for example.</p>',16)]))}const m=a(s,[["render",r]]);export{f as __pageData,m as default};
