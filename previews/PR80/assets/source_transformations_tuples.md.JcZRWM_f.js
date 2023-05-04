import{_ as s,c as a,o as n,a7 as i}from"./chunks/framework.DqBJNxpO.js";const u=JSON.parse('{"title":"Tuple conversion","description":"","frontmatter":{},"headers":[],"relativePath":"source/transformations/tuples.md","filePath":"source/transformations/tuples.md","lastUpdated":null}'),e={name:"source/transformations/tuples.md"},t=i(`<h1 id="Tuple-conversion" tabindex="-1">Tuple conversion <a class="header-anchor" href="#Tuple-conversion" aria-label="Permalink to &quot;Tuple conversion {#Tuple-conversion}&quot;">​</a></h1><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    tuples(obj)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Convert all points in \`obj\` to \`Tuple\`s, wherever the are nested.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Returns a similar object or collection of objects using GeoInterface.jl</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">geometries wrapping \`Tuple\` points.</span></span></code></pre></div><p>Keywords</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">$</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">APPLY_KEYWORDS</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">function tuples(geom; kw...)</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    if _is3d(geom)</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        return apply(PointTrait(), geom; kw...) do p</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">            (Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p)))</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        end</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    else</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        return apply(PointTrait(), geom; kw...) do p</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">            (Float64(GI.x(p)), Float64(GI.y(p)))</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        end</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    end</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">end</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,6),l=[t];function p(o,r,h,c,k,d){return n(),a("div",null,l)}const g=s(e,[["render",p]]);export{u as __pageData,g as default};
