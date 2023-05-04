import{_ as a,c as i,o as n,az as e}from"./chunks/framework.u4vpEvHV.js";const c=JSON.parse('{"title":"Types","description":"","frontmatter":{},"headers":[],"relativePath":"source/src/types.md","filePath":"source/src/types.md","lastUpdated":null}'),t={name:"source/src/types.md"};function l(p,s,h,k,r,o){return n(),i("div",null,s[0]||(s[0]=[e(`<h1 id="types" tabindex="-1">Types <a class="header-anchor" href="#types" aria-label="Permalink to &quot;Types&quot;">​</a></h1><p>This defines core types that the GeometryOps ecosystem uses, and that are usable in more than just GeometryOps.</p><h2 id="Manifold" tabindex="-1"><code>Manifold</code> <a class="header-anchor" href="#Manifold" aria-label="Permalink to &quot;\`Manifold\` {#Manifold}&quot;">​</a></h2><p>A manifold is mathematically defined as a topological space that resembles Euclidean space locally.</p><p>In GeometryOps (and geodesy more generally), there are three manifolds we care about:</p><ul><li><p><code>Planar</code>: the 2d plane, a completely Euclidean manifold</p></li><li><p><code>Spherical</code>: the unit sphere, but one where areas are multiplied by the radius of the Earth. This is not Euclidean globally, but all map projections attempt to represent the sphere on the Euclidean 2D plane to varying degrees of success.</p></li><li><p><code>Geodesic</code>: the ellipsoid, the closest we can come to representing the Earth by a simple geometric shape. Parametrized by <code>semimajor_axis</code> and <code>inv_flattening</code>.</p></li></ul><p>Generally, we aim to have <code>Linear</code> and <code>Spherical</code> be operable everywhere, whereas <code>Geodesic</code> will only apply in specific circumstances. Currently, those circumstances are <code>area</code> and <code>segmentize</code>, but this could be extended with time and <a href="https://github.com/JuliaGeo/SphericalGeodesics.jl" target="_blank" rel="noreferrer">https://github.com/JuliaGeo/SphericalGeodesics.jl</a>.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Planar, Spherical, Geodesic</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> TraitTarget</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> BoolsAsTypes, True, False, booltype</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    abstract type Manifold</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">A manifold is mathematically defined as a topological space that resembles Euclidean space locally.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">We use the manifold definition to define the space in which an operation should be performed, or where a geometry lies.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Currently we have \`Planar\`, \`Spherical\`, and \`Geodesic\` manifolds.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">abstract type</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Manifold </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    Planar()</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">A planar manifold refers to the 2D Euclidean plane.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Z coordinates may be accepted but will not influence geometry calculations, which</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">are done purely on 2D geometry.  This is the standard &quot;2.5D&quot; model used by e.g. GEOS.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Planar </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Manifold</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    Spherical(; radius)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">A spherical manifold means that the geometry is on the 3-sphere (but is represented by 2-D longitude and latitude).</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">By default, the radius is defined as the mean radius of the Earth, \`6371008.8 m\`.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"># Extended help</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">!!! note</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    The traditional definition of spherical coordinates in physics and mathematics,</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    \`\`r, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\\\</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">theta, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\\\</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">phi\`\`, uses the _colatitude_, that measures angular displacement from the \`z\`-axis.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    Here, we use the geographic definition of longitude and latitude, meaning</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    that \`lon\` is longitude between -180 and 180, and \`lat\` is latitude between</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    \`-90\` (south pole) and \`90\` (north pole).</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Base</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@kwdef</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Spherical{T} </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Manifold</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    radius</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">T</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 6371008.8</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    Geodesic(; semimajor_axis, inv_flattening)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">A geodesic manifold means that the geometry is on a 3-dimensional ellipsoid, parameterized by \`semimajor_axis\` (\`\`a\`\` in mathematical parlance)</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">and \`inv_flattening\` (\`\`1/f\`\`).</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Usually, this is only relevant for area and segmentization calculations.  It becomes more relevant as one grows closer to the poles (or equator).</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Base</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@kwdef</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Geodesic{T} </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Manifold</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    semimajor_axis</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">T</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 6378137.0</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    inv_flattening</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">T</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 298.257223563</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><h2 id="TraitTarget" tabindex="-1"><code>TraitTarget</code> <a class="header-anchor" href="#TraitTarget" aria-label="Permalink to &quot;\`TraitTarget\` {#TraitTarget}&quot;">​</a></h2><p>This struct holds a trait parameter or a union of trait parameters. It&#39;s essentially a way to construct unions.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    TraitTarget{T}</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">This struct holds a trait parameter or a union of trait parameters.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">It is primarily used for dispatch into methods which select trait levels,</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">like \`apply\`, or as a parameter to \`target\`.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"># Constructors</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`julia</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">TraitTarget(GI.PointTrait())</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">TraitTarget(GI.LineStringTrait(), GI.LinearRingTrait()) # and other traits as you may like</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">TraitTarget(TraitTarget(...))</span></span></code></pre></div><p>There are also type based constructors available, but that&#39;s not advised.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">TraitTarget</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">PointTrait)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">TraitTarget</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(Union{GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">LineStringTrait, GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">LinearRingTrait})</span></span></code></pre></div><p>etc.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">struct TraitTarget{T} end</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">TraitTarget(::Type{T}) where T = TraitTarget{T}()</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">TraitTarget(::T) where T&lt;:GI.AbstractTrait = TraitTarget{T}()</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">TraitTarget(::TraitTarget{T}) where T = TraitTarget{T}()</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">TraitTarget(::Type{&lt;:TraitTarget{T}}) where T = TraitTarget{T}()</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">TraitTarget(traits::GI.AbstractTrait...) = TraitTarget{Union{map(typeof, traits)...}}()</span></span>
<span class="line"></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Base.in(::Trait, ::TraitTarget{Target}) where {Trait &lt;: GI.AbstractTrait, Target} = Trait &lt;: Target</span></span></code></pre></div><h2 id="BoolsAsTypes" tabindex="-1"><code>BoolsAsTypes</code> <a class="header-anchor" href="#BoolsAsTypes" aria-label="Permalink to &quot;\`BoolsAsTypes\` {#BoolsAsTypes}&quot;">​</a></h2><p>In <code>apply</code> and <code>applyreduce</code>, we pass <code>threading</code> and <code>calc_extent</code> as types, not simple boolean values.</p><p>This is to help compilation - with a type to hold on to, it&#39;s easier for the compiler to separate threaded and non-threaded code paths.</p><p>Note that if we didn&#39;t include the parent abstract type, this would have been really type unstable, since the compiler couldn&#39;t tell what would be returned!</p><p>We had to add the type annotation on the <code>booltype(::Bool)</code> method for this reason as well.</p><p>TODO: should we switch to <code>Static.jl</code>?</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    abstract type BoolsAsTypes</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">abstract type</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> BoolsAsTypes </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    struct True &lt;: BoolsAsTypes</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">A struct that means \`true\`.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> True </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> BoolsAsTypes</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    struct False &lt;: BoolsAsTypes</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">A struct that means \`false\`.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> False </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> BoolsAsTypes</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    booltype(x)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Returns a \`BoolsAsTypes\` from \`x\`, whether it&#39;s a boolean or a BoolsAsTypes.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> booltype </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@inline</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> booltype</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Bool</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">BoolsAsTypes</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">?</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> True</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> False</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@inline</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> booltype</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">BoolsAsTypes</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">BoolsAsTypes</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,24)]))}const F=a(t,[["render",l]]);export{c as __pageData,F as default};
