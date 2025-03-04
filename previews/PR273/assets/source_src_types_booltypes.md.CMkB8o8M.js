import{_ as a,c as i,o as n,az as t}from"./chunks/framework.Ct5lodKo.js";const c=JSON.parse('{"title":"BoolsAsTypes","description":"","frontmatter":{},"headers":[],"relativePath":"source/src/types/booltypes.md","filePath":"source/src/types/booltypes.md","lastUpdated":null}'),e={name:"source/src/types/booltypes.md"};function l(p,s,h,k,o,r){return n(),i("div",null,s[0]||(s[0]=[t(`<h1 id="BoolsAsTypes" tabindex="-1"><code>BoolsAsTypes</code> <a class="header-anchor" href="#BoolsAsTypes" aria-label="Permalink to &quot;\`BoolsAsTypes\` {#BoolsAsTypes}&quot;">​</a></h1><p>In <code>apply</code> and <code>applyreduce</code>, we pass <code>threading</code> and <code>calc_extent</code> as types, not simple boolean values.</p><p>This is to help compilation - with a type to hold on to, it&#39;s easier for the compiler to separate threaded and non-threaded code paths.</p><p>Note that if we didn&#39;t include the parent abstract type, this would have been really type unstable, since the compiler couldn&#39;t tell what would be returned!</p><p>We had to add the type annotation on the <code>booltype(::Bool)</code> method for this reason as well.</p><p>!!! note Static.jl</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Static.jl is a package that provides a way to store and manipulate static values.</span></span>
<span class="line"><span>But it creates a lot of invalidations since it breaks the assumption that operations</span></span>
<span class="line"><span>like \`&lt;\`, \`&gt;\` and \`==\` can only return booleans.  So we don&#39;t use it here.</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> BoolsAsTypes, True, False, booltype</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
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
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@inline</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> booltype</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">BoolsAsTypes</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">BoolsAsTypes</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> =</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> x</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,10)]))}const y=a(e,[["render",l]]);export{c as __pageData,y as default};
