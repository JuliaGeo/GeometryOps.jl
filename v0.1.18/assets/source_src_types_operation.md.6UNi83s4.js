import{_ as i,c as a,o as n,az as e}from"./chunks/framework.BfKjXwMd.js";const g=JSON.parse('{"title":"Operations","description":"","frontmatter":{},"headers":[],"relativePath":"source/src/types/operation.md","filePath":"source/src/types/operation.md","lastUpdated":null}'),t={name:"source/src/types/operation.md"};function l(p,s,h,r,o,k){return n(),a("div",null,s[0]||(s[0]=[e(`<h1 id="operations" tabindex="-1">Operations <a class="header-anchor" href="#operations" aria-label="Permalink to &quot;Operations&quot;">​</a></h1><blockquote><p><strong>Warning</strong></p><p>Operations are not yet implemented. If you want to implement them then you may do so at your own risk - or file a PR!</p></blockquote><p>Operations are callable structs, that contain the entire specification for what the algorithm will do.</p><p>Sometimes they may be underspecified and only materialized fully when you see the geometry, so you can extract the best manifold for those geometries.</p><ul><li><p>Some conceptual thing that you do to a geometry</p></li><li><p>Overloads on abstract type to decompose user input to have materialized algorithm, manifold, and geoms</p></li><li><p>Run Operation{Alg{Manifold}}(trait, geom) at the lowest level</p></li><li><p>Some indication on whether to use apply or applyreduce? Or are we going too far here</p><ul><li>if we do this, then we also need <code>operation_level</code> to return a geometry trait or traittarget</li></ul></li></ul><p>Operations may look like:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Arclength</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()(geoms)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Arclength</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Geodesic</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())(geoms)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Arclength</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Proj</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())(geoms)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Arclength</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Proj</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Geodesic</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(; </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)))(geoms)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Arclength</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Ericsson</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())(geoms) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># more precise, goes wonky if any points in a triangle are antipodal</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Arclength</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">LHuilier</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">())(geoms) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># less precise, does not go wonky on antipodal points</span></span></code></pre></div><p>Two argument operations, like polygon set operations, may look like:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Union</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">intersection_alg</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(manifold); exact, target)(geom1, geom2)</span></span></code></pre></div><p>Here intersection_alg can be Foster, which we already have in GeometryOps, or GEOS but if we ever implement e.g. RelateNG in GeometryOps, we can add that in.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    abstract type Operation{Alg &lt;: Algorithm} end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Operations are callable structs, that contain the entire specification for what the algorithm will do.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Sometimes they may be underspecified and only materialized fully when you see the geometry, so you can extract</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">the best manifold for those geometries.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">abstract type</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Operation{Alg </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Algorithm</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">} </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#=</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">Here&#39;s an example of how this might work:</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">\`\`\`julia</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">struct XPlusOneOperation{M &lt;: Manifold} &lt;: Operation{NoAlgorithm{M}}</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">    m::M # the manifold always needs to be stored, since it&#39;s not a singleton</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">    x::Int</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">end</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">\`\`\`</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">\`\`\`julia</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">struct Area{Alg&lt;: Algorithm} &lt;: Operation{Alg}</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">    alg::Alg</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">Area() = Area(NoAlgorithm(Planar()))</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">function (op::Area{Alg})(data; threaded = _False(), init = 0.0) where {Alg &lt;: Algorithm}</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">    return GO.area(op.alg, data; threaded, init)</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">end</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">\`\`\`</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">=#</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,13)]))}const c=i(t,[["render",l]]);export{g as __pageData,c as default};
