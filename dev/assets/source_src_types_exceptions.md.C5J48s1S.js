import{_ as i,c as a,o as n,az as t}from"./chunks/framework.paeFYVqm.js";const E=JSON.parse('{"title":"Errors and exceptions","description":"","frontmatter":{},"headers":[],"relativePath":"source/src/types/exceptions.md","filePath":"source/src/types/exceptions.md","lastUpdated":null}'),h={name:"source/src/types/exceptions.md"};function l(p,s,k,e,r,d){return n(),a("div",null,s[0]||(s[0]=[t(`<h1 id="Errors-and-exceptions" tabindex="-1">Errors and exceptions <a class="header-anchor" href="#Errors-and-exceptions" aria-label="Permalink to &quot;Errors and exceptions {#Errors-and-exceptions}&quot;">​</a></h1><p>We create a few custom exception types in this file, that have nice show methods that we can use for certain errors.</p><p>This makes it substantially easier to catch specific kinds of errors and show them. For example, we can catch <code>WrongManifoldException</code> and show a nice error message, and error hinters can be specialized to that as well.</p><div class="warning custom-block"><p class="custom-block-title">Missing docstring.</p><p>Missing docstring for <code>GeometryOpsCore.WrongManifoldException</code>. Check Documenter&#39;s build log for details.</p></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> WrongManifoldException</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    WrongManifoldException{InputManifold, DesiredManifold, Algorithm} &lt;: Exception</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">This error is thrown when an \`Algorithm\` is called with a manifold that it was not designed for.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">It&#39;s mainly thrown when constructing \`SingleManifoldAlgorithm\` types.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">struct</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> WrongManifoldException{InputManifold, DesiredManifold, Algorithm} </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&lt;:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Base.Exception</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    description</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">String</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">WrongManifoldException</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">{I, D, A}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">where</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> {I, D, A} </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> WrongManifoldException{I, D, A}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Base</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;">showerror</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">IO</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, e</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">WrongManifoldException{I,D,A}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">where</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> {I,D,A}</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    print</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Algorithm &quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    printstyled</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, A; bold </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :green</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    print</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot; is only compatible with manifold &quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    printstyled</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, D; bold </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :blue</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    print</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;,</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\n</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;"> but it was called with manifold &quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    printstyled</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, I; bold </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, color </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> :red</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    print</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;.&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    println</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    \\n</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    To fix this issue, you can specify the manifold explicitly,</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    e.g. \`</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$A</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">$D</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">(); kwargs...)\`, when constructing the algorithm.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    &quot;&quot;&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> !</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">isempty</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(e</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">description)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        print</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, </span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">\\n\\n</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        print</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(io, e</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">description)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span></code></pre></div><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,7)]))}const g=i(h,[["render",l]]);export{E as __pageData,g as default};
