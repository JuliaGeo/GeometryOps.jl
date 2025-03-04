import{_ as o,C as l,c as p,o as h,az as a,j as i,G as t,a as n,w as r}from"./chunks/framework.aIRhVz5C.js";const _=JSON.parse('{"title":"Geometry Corrections","description":"","frontmatter":{},"headers":[],"relativePath":"source/transformations/correction/geometry_correction.md","filePath":"source/transformations/correction/geometry_correction.md","lastUpdated":null}'),k={name:"source/transformations/correction/geometry_correction.md"},c={class:"jldocstring custom-block",open:""},d={class:"jldocstring custom-block",open:""},y={class:"jldocstring custom-block",open:""},g={class:"jldocstring custom-block",open:""},m={class:"jldocstring custom-block",open:""};function E(u,s,C,F,f,b){const e=l("Badge");return h(),p("div",null,[s[20]||(s[20]=a('<h1 id="Geometry-Corrections" tabindex="-1">Geometry Corrections <a class="header-anchor" href="#Geometry-Corrections" aria-label="Permalink to &quot;Geometry Corrections {#Geometry-Corrections}&quot;">​</a></h1><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> fix</span></span></code></pre></div><p>This file simply defines the <code>GeometryCorrection</code> abstract type, and the interface that any <code>GeometryCorrection</code> must implement.</p><p>A geometry correction is a transformation that is applied to a geometry to correct it in some way.</p><p>For example, a <code>ClosedRing</code> correction might be applied to a <code>Polygon</code> to ensure that its exterior ring is closed.</p><h2 id="interface" tabindex="-1">Interface <a class="header-anchor" href="#interface" aria-label="Permalink to &quot;Interface&quot;">​</a></h2><p>All <code>GeometryCorrection</code>s are callable structs which, when called, apply the correction to the given geometry, and return either a copy or the original geometry (if nothing needed to be corrected).</p><p>See below for the full interface specification.</p>',8)),i("details",c,[i("summary",null,[s[0]||(s[0]=i("a",{id:"GeometryOps.GeometryCorrection-source-transformations-correction-geometry_correction",href:"#GeometryOps.GeometryCorrection-source-transformations-correction-geometry_correction"},[i("span",{class:"jlbinding"},"GeometryOps.GeometryCorrection")],-1)),s[1]||(s[1]=n()),t(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),s[3]||(s[3]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">abstract type</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeometryCorrection</span></span></code></pre></div><p>This abstract type represents a geometry correction.</p><p><strong>Interface</strong></p><p>Any <code>GeometryCorrection</code> must implement two functions: * <code>application_level(::GeometryCorrection)::AbstractGeometryTrait</code>: This function should return the <code>GeoInterface</code> trait that the correction is intended to be applied to, like <code>PointTrait</code> or <code>LineStringTrait</code> or <code>PolygonTrait</code>. * <code>(::GeometryCorrection)(::AbstractGeometryTrait, geometry)::(some_geometry)</code>: This function should apply the correction to the given geometry, and return a new geometry.</p>',4)),t(e,{type:"info",class:"source-link",text:"source"},{default:r(()=>s[2]||(s[2]=[i("a",{href:"https://github.com/JuliaGeo/GeometryOps.jl/blob/cb1a04e4a5c685bedf49aed43f13093d0ec3c781/src/transformations/correction/geometry_correction.jl#L28-L38",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),s[21]||(s[21]=a(`<p>Any geometry correction must implement the interface as given above.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    abstract type GeometryCorrection</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">This abstract type represents a geometry correction.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"># Interface</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Any \`GeometryCorrection\` must implement two functions:</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    * \`application_level(::GeometryCorrection)::AbstractGeometryTrait\`: This function should return the \`GeoInterface\` trait that the correction is intended to be applied to, like \`PointTrait\` or \`LineStringTrait\` or \`PolygonTrait\`.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    * \`(::GeometryCorrection)(::AbstractGeometryTrait, geometry)::(some_geometry)\`: This function should apply the correction to the given geometry, and return a new geometry.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">abstract type</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeometryCorrection </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">application_level</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(gc</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GeometryCorrection</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> error</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Not implemented yet for </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(gc)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(gc</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GeometryCorrection</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)(geometry) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> gc</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">trait</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(geometry), geometry)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(gc</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GeometryCorrection</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)(trait</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GI.AbstractGeometryTrait</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, geometry) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> error</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Not implemented yet for </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(gc)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> and </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(trait)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">.&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> fix</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(geometry; corrections </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeometryCorrection[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ClosedRing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(),], kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    traits </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> application_level</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">.(corrections)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    final_geometry </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> geometry</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Trait </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PointTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">MultiPointTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">LineStringTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">LinearRingTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">MultiLineStringTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PolygonTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">MultiPolygonTrait)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        available_corrections </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> findall</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">==</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Trait, traits)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        isempty</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(available_corrections) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;&amp;</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> continue</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        @debug</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> &quot;Correcting for </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(Trait)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        net_function </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> reduce</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">∘</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, corrections[available_corrections])</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">        final_geometry </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> apply</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(net_function, Trait, final_geometry; kwargs</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> final_geometry</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><h2 id="Available-corrections" tabindex="-1">Available corrections <a class="header-anchor" href="#Available-corrections" aria-label="Permalink to &quot;Available corrections {#Available-corrections}&quot;">​</a></h2>`,3)),i("details",d,[i("summary",null,[s[4]||(s[4]=i("a",{id:"GeometryOps.ClosedRing-source-transformations-correction-geometry_correction",href:"#GeometryOps.ClosedRing-source-transformations-correction-geometry_correction"},[i("span",{class:"jlbinding"},"GeometryOps.ClosedRing")],-1)),s[5]||(s[5]=n()),t(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),s[7]||(s[7]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ClosedRing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> GeometryCorrection</span></span></code></pre></div><p>This correction ensures that a polygon&#39;s exterior and interior rings are closed.</p><p>It can be called on any geometry correction as usual.</p><p>See also <a href="/GeometryOps.jl/previews/PR271/api#GeometryOps.GeometryCorrection"><code>GeometryCorrection</code></a>.</p>',4)),t(e,{type:"info",class:"source-link",text:"source"},{default:r(()=>s[6]||(s[6]=[i("a",{href:"https://github.com/JuliaGeo/GeometryOps.jl/blob/cb1a04e4a5c685bedf49aed43f13093d0ec3c781/src/transformations/correction/closed_ring.jl#L38-L46",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),i("details",y,[i("summary",null,[s[8]||(s[8]=i("a",{id:"GeometryOps.DiffIntersectingPolygons-source-transformations-correction-geometry_correction",href:"#GeometryOps.DiffIntersectingPolygons-source-transformations-correction-geometry_correction"},[i("span",{class:"jlbinding"},"GeometryOps.DiffIntersectingPolygons")],-1)),s[9]||(s[9]=n()),t(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),s[11]||(s[11]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">DiffIntersectingPolygons</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> GeometryCorrection</span></span></code></pre></div><p>This correction ensures that the polygons included in a multipolygon aren&#39;t intersecting. If any polygon&#39;s are intersecting, they will be made nonintersecting through the <a href="/GeometryOps.jl/previews/PR271/api#GeometryOps.difference-Union{Tuple{T}, Tuple{Any, Any}, Tuple{Any, Any, Type{T}}} where T&lt;:AbstractFloat"><code>difference</code></a> operation to create a unique set of disjoint (other than potentially connections by a single point) polygons covering the same area. See also <a href="/GeometryOps.jl/previews/PR271/api#GeometryOps.GeometryCorrection"><code>GeometryCorrection</code></a>, <a href="/GeometryOps.jl/previews/PR271/api#GeometryOps.UnionIntersectingPolygons"><code>UnionIntersectingPolygons</code></a>.</p>',2)),t(e,{type:"info",class:"source-link",text:"source"},{default:r(()=>s[10]||(s[10]=[i("a",{href:"https://github.com/JuliaGeo/GeometryOps.jl/blob/cb1a04e4a5c685bedf49aed43f13093d0ec3c781/src/transformations/correction/intersecting_polygons.jl#L92-L99",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),i("details",g,[i("summary",null,[s[12]||(s[12]=i("a",{id:"GeometryOps.GeometryCorrection-source-transformations-correction-geometry_correction-2",href:"#GeometryOps.GeometryCorrection-source-transformations-correction-geometry_correction-2"},[i("span",{class:"jlbinding"},"GeometryOps.GeometryCorrection")],-1)),s[13]||(s[13]=n()),t(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),s[15]||(s[15]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">abstract type</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeometryCorrection</span></span></code></pre></div><p>This abstract type represents a geometry correction.</p><p><strong>Interface</strong></p><p>Any <code>GeometryCorrection</code> must implement two functions: * <code>application_level(::GeometryCorrection)::AbstractGeometryTrait</code>: This function should return the <code>GeoInterface</code> trait that the correction is intended to be applied to, like <code>PointTrait</code> or <code>LineStringTrait</code> or <code>PolygonTrait</code>. * <code>(::GeometryCorrection)(::AbstractGeometryTrait, geometry)::(some_geometry)</code>: This function should apply the correction to the given geometry, and return a new geometry.</p>',4)),t(e,{type:"info",class:"source-link",text:"source"},{default:r(()=>s[14]||(s[14]=[i("a",{href:"https://github.com/JuliaGeo/GeometryOps.jl/blob/cb1a04e4a5c685bedf49aed43f13093d0ec3c781/src/transformations/correction/geometry_correction.jl#L28-L38",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),i("details",m,[i("summary",null,[s[16]||(s[16]=i("a",{id:"GeometryOps.UnionIntersectingPolygons-source-transformations-correction-geometry_correction",href:"#GeometryOps.UnionIntersectingPolygons-source-transformations-correction-geometry_correction"},[i("span",{class:"jlbinding"},"GeometryOps.UnionIntersectingPolygons")],-1)),s[17]||(s[17]=n()),t(e,{type:"info",class:"jlObjectType jlType",text:"Type"})]),s[19]||(s[19]=a('<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">UnionIntersectingPolygons</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> GeometryCorrection</span></span></code></pre></div><p>This correction ensures that the polygon&#39;s included in a multipolygon aren&#39;t intersecting. If any polygon&#39;s are intersecting, they will be combined through the union operation to create a unique set of disjoint (other than potentially connections by a single point) polygons covering the same area.</p><p>See also <a href="/GeometryOps.jl/previews/PR271/api#GeometryOps.GeometryCorrection"><code>GeometryCorrection</code></a>.</p>',3)),t(e,{type:"info",class:"source-link",text:"source"},{default:r(()=>s[18]||(s[18]=[i("a",{href:"https://github.com/JuliaGeo/GeometryOps.jl/blob/cb1a04e4a5c685bedf49aed43f13093d0ec3c781/src/transformations/correction/intersecting_polygons.jl#L47-L56",target:"_blank",rel:"noreferrer"},"source",-1)])),_:1})]),s[22]||(s[22]=i("hr",null,null,-1)),s[23]||(s[23]=i("p",null,[i("em",null,[n("This page was generated using "),i("a",{href:"https://github.com/fredrikekre/Literate.jl",target:"_blank",rel:"noreferrer"},"Literate.jl"),n(".")])],-1))])}const G=o(k,[["render",E]]);export{_ as __pageData,G as default};
