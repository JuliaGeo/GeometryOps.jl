import{_ as a,c as t,o as n,az as e}from"./chunks/framework.Cj-4H1V0.js";const k=JSON.parse('{"title":"TaskFunctors","description":"","frontmatter":{},"headers":[],"relativePath":"source/src/types/threading.md","filePath":"source/src/types/threading.md","lastUpdated":null}'),i={name:"source/src/types/threading.md"};function o(r,s,c,p,l,h){return n(),t("div",null,s[0]||(s[0]=[e(`<h1 id="taskfunctors" tabindex="-1">TaskFunctors <a class="header-anchor" href="#taskfunctors" aria-label="Permalink to &quot;TaskFunctors&quot;">​</a></h1><p>Many functions a user might run through <code>apply</code> call into C. But what do you do when the C function you are calling is not threadsafe, or uses a reentrant/context based API?</p><p>The answer is to wrap the function in TaskFunctors, which will use an approximation of task local storage to ensure that each task calls its own copy of the function or C object you are invoking.</p><p>The primary application of this is Proj <code>reproject</code>, and Proj has a reentrant API based on context objects. So, to make this work, we need to create a new context object per task. We do this and pass the vector of Proj transformation objects to <code>TaskFunctors</code>, which <code>apply</code> and <code>_maptasks</code> dispatch on to behave correctly in this circumstance.</p><div class="warning custom-block"><p class="custom-block-title">Missing docstring.</p><p>Missing docstring for <code>TaskFunctors</code>. Check Documenter&#39;s build log for details.</p></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    TaskFunctors(functors)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">A struct to hold the functors,</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">for internal use with \`_maptasks\` where functions have state</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">that cannot be shared accross threads, such as \`Proj.Transformation\`.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`functors\` must be an array or tuple of functors, one per task.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">The number of tasks is the number of elements in \`functors\`.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> TaskFunctors{F}</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    functors</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">F</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,8)]))}const u=a(i,[["render",o]]);export{k as __pageData,u as default};
