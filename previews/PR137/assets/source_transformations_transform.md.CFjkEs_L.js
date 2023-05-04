import{_ as s,c as i,o as a,a6 as n}from"./chunks/framework.ICD8RgvV.js";const d=JSON.parse('{"title":"Pointwise transformation","description":"","frontmatter":{},"headers":[],"relativePath":"source/transformations/transform.md","filePath":"source/transformations/transform.md","lastUpdated":null}'),t={name:"source/transformations/transform.md"},l=n(`<h1 id="Pointwise-transformation" tabindex="-1">Pointwise transformation <a class="header-anchor" href="#Pointwise-transformation" aria-label="Permalink to &quot;Pointwise transformation {#Pointwise-transformation}&quot;">​</a></h1><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    transform(f, obj)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Apply a function \`f\` to all the points in \`obj\`.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Points will be passed to \`f\` as an \`SVector\` to allow</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">using CoordinateTransformations.jl and Rotations.jl</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">without hassle.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`SVector\` is also a valid GeoInterface.jl point, so will</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">work in all GeoInterface.jl methods.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"># Example</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`julia</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">julia&gt; import GeoInterface as GI</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">julia&gt; import GeometryOps as GO</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">julia&gt; geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])]);</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">julia&gt; f = CoordinateTransformations.Translation(3.5, 1.5)</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Translation(3.5, 1.5)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">julia&gt; GO.transform(f, geom)</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">GeoInterface.Wrappers.Polygon{false, false, Vector{GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Float64}}, Nothing, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.Linea</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">rRing{false, false, Vector{StaticArraysCore.SVector{2, Float64}}, Nothing, Nothing}[GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Float64}}, Nothing, Nothing}(StaticArraysCo</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">re.SVector{2, Float64}[[4.5, 3.5], [6.5, 5.5], [8.5, 7.5], [4.5, 3.5]], nothing, nothing), GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Float64}}, Nothing, Nothing}(StaticA</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">rraysCore.SVector{2, Float64}[[6.5, 5.5], [8.5, 7.5], [9.5, 8.5], [6.5, 5.5]], nothing, nothing)], nothing, nothing)</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">With Rotations.jl you need to actuall multiply the Rotation</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">by the \`SVector\` point, which is easy using an anonymous function.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`julia</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">julia&gt; using Rotations</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">julia&gt; GO.transform(p -&gt; one(RotMatrix{2}) * p, geom)</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">GeoInterface.Wrappers.Polygon{false, false, Vector{GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Int64}}, Nothing, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.LinearR</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ing{false, false, Vector{StaticArraysCore.SVector{2, Int64}}, Nothing, Nothing}[GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Int64}}, Nothing, Nothing}(StaticArraysCore.SVe</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">ctor{2, Int64}[[2, 1], [4, 3], [6, 5], [2, 1]], nothing, nothing), GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Int64}}, Nothing, Nothing}(StaticArraysCore.SVector{2, Int64</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">}[[4, 3], [6, 5], [7, 6], [4, 3]], nothing, nothing)], nothing, nothing)</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> transform</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(f, geom; kw</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> _is3d</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(geom)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        return</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> apply</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">PointTrait</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), geom; kw</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">do</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> p</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">            f</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(StaticArrays</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SVector{3}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">((GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p), GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">y</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p), GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">z</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p))))</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    else</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        return</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> apply</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">PointTrait</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(), geom; kw</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">do</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> p</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">            f</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(StaticArrays</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">SVector{2}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">((GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">x</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p), GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">y</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p))))</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,4),e=[l];function p(h,r,k,o,g,F){return a(),i("div",null,e)}const y=s(t,[["render",p]]);export{d as __pageData,y as default};
