import{_ as s,c as i,o as t,az as n}from"./chunks/framework.D4Qd0f0J.js";const c=JSON.parse('{"title":"","description":"","frontmatter":{},"headers":[],"relativePath":"source/src/types/traittarget.md","filePath":"source/src/types/traittarget.md","lastUpdated":null}'),e={name:"source/src/types/traittarget.md"};function r(l,a,p,h,k,T){return t(),i("div",null,a[0]||(a[0]=[n(`<h2 id="TraitTarget" tabindex="-1"><code>TraitTarget</code> <a class="header-anchor" href="#TraitTarget" aria-label="Permalink to &quot;\`TraitTarget\` {#TraitTarget}&quot;">​</a></h2><p>This struct holds a trait parameter or a union of trait parameters. It&#39;s essentially a way to construct unions.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> TraitTarget</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
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
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Base.in(::Trait, ::TraitTarget{Target}) where {Trait &lt;: GI.AbstractTrait, Target} = Trait &lt;: Target</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,9)]))}const o=s(e,[["render",r]]);export{c as __pageData,o as default};
