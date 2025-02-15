import{_ as i,c as a,o as n,az as e}from"./chunks/framework.BNXt_2zd.js";const E=JSON.parse('{"title":"Types","description":"","frontmatter":{},"headers":[],"relativePath":"source/types.md","filePath":"source/types.md","lastUpdated":null}'),t={name:"source/types.md"};function l(h,s,p,k,r,d){return n(),a("div",null,s[0]||(s[0]=[e(`<h1 id="types" tabindex="-1">Types <a class="header-anchor" href="#types" aria-label="Permalink to &quot;Types&quot;">​</a></h1><p>This file defines some fundamental types used in GeometryOps.</p><blockquote><p><strong>Warning</strong></p><p>Unlike in other Julia packages, only some types are defined in this file, not all. This is because we define types in the files where they are used, to make it easier to understand the code.</p></blockquote><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GEOS</span></span></code></pre></div><h2 id="GEOS" tabindex="-1"><code>GEOS</code> <a class="header-anchor" href="#GEOS" aria-label="Permalink to &quot;\`GEOS\` {#GEOS}&quot;">​</a></h2><p><code>GEOS</code> is a struct which instructs the method it&#39;s passed to as an algorithm to use the appropriate GEOS function via <code>LibGEOS.jl</code> for the operation.</p><p>It&#39;s generally a lot slower than the native Julia implementations, but it&#39;s useful for two reasons:</p><ol><li><p>Functionality which doesn&#39;t exist in GeometryOps can be accessed through the GeometryOps API, but use GEOS in the backend until someone implements a native Julia version.</p></li><li><p>It&#39;s a good way to test the correctness of the native implementations.</p></li></ol><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    GEOS(; params...)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">A struct which instructs the method it&#39;s passed to as an algorithm</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">to use the appropriate GEOS function via \`LibGEOS.jl\` for the operation.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Dispatch is generally carried out using the names of the keyword arguments.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">For example, \`segmentize\` will only accept a \`GEOS\` struct with only a</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`max_distance\` keyword, and no other.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">It&#39;s generally a lot slower than the native Julia implementations, since</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">it must convert to the LibGEOS implementation and back - so be warned!</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GEOS</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    params</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">NamedTuple</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> GEOS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(; params</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    nt </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> NamedTuple</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(params)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> GEOS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(nt)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><p>These are definitions for convenience, so we don&#39;t have to type out <code>alg.params</code> every time.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Base</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">get</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(alg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GEOS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, key, value) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Base</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(alg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">params, key, value)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Base</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">get</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(f</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, alg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GEOS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, key) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Base</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">get</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(f, alg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">params, key)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    enforce(alg::GO.GEOS, kw::Symbol, f)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Enforce the presence of a keyword argument in a \`GEOS\` algorithm, and return \`alg.params[kw]\`.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Throws an error if the key is not present, and mentions \`f\` in the error message (since there isn&#39;t</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">a good way to get the name of the function that called this method).</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> enforce</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(alg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">GEOS</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, kw</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Symbol</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, f)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> haskey</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(alg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">params, kw)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> alg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">params[kw]</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    else</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        error</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(f)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> requires a \`</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$(kw)</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\` keyword argument to the \`GEOS\` algorithm, which was not provided.&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,13)]))}const g=i(t,[["render",l]]);export{E as __pageData,g as default};
