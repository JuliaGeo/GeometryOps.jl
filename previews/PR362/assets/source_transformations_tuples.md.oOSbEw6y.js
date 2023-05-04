import{_ as a,c as n,o as i,az as e}from"./chunks/framework.e25qLPaf.js";const d=JSON.parse('{"title":"Tuple conversion","description":"","frontmatter":{},"headers":[],"relativePath":"source/transformations/tuples.md","filePath":"source/transformations/tuples.md","lastUpdated":null}'),t={name:"source/transformations/tuples.md"};function p(l,s,r,o,h,k){return i(),n("div",null,[...s[0]||(s[0]=[e(`<h1 id="Tuple-conversion" tabindex="-1">Tuple conversion <a class="header-anchor" href="#Tuple-conversion" aria-label="Permalink to &quot;Tuple conversion {#Tuple-conversion}&quot;">​</a></h1><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    tuples(obj)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Convert all points in \`obj\` to \`Tuple\`s, wherever the are nested.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Returns a similar object or collection of objects using GeoInterface.jl</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">geometries wrapping \`Tuple\` points.</span></span></code></pre></div><p>Keywords</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">$</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">APPLY_KEYWORDS</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">function tuples(geom, ::Type{T} = Float64; kw...) where T</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    if _is3d(geom)</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        return apply(PointTrait(), geom; kw...) do p</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">            (T(GI.x(p)), T(GI.y(p)), T(GI.z(p)))</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        end</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    else</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        return apply(PointTrait(), geom; kw...) do p</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">            (T(GI.x(p)), T(GI.y(p)))</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        end</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    end</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">end</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,6)])])}const F=a(t,[["render",p]]);export{d as __pageData,F as default};
