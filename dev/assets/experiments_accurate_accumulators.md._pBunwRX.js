import{_ as a,c as e,o as t,a6 as s}from"./chunks/framework.B9djhgVS.js";const g=JSON.parse('{"title":"Accurate accumulation","description":"","frontmatter":{},"headers":[],"relativePath":"experiments/accurate_accumulators.md","filePath":"experiments/accurate_accumulators.md","lastUpdated":null}'),n={name:"experiments/accurate_accumulators.md"},c=s(`<h1 id="Accurate-accumulation" tabindex="-1">Accurate accumulation <a class="header-anchor" href="#Accurate-accumulation" aria-label="Permalink to &quot;Accurate accumulation {#Accurate-accumulation}&quot;">​</a></h1><p>Accurate arithmetic is a technique which allows you to calculate using more precision than the provided numeric type.</p><p>We will use the accurate sum routines from <a href="https://github.com/JuliaMath/AccurateArithmetic.jl" target="_blank" rel="noreferrer">AccurateArithmetic.jl</a> to show the difference!</p><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span>import GeometryOps as GO, GeoInterface as GI</span></span>
<span class="line"><span>using GeoJSON</span></span>
<span class="line"><span>using AccurateArithmetic</span></span>
<span class="line"><span>using NaturalEarth</span></span>
<span class="line"><span></span></span>
<span class="line"><span>all_adm0 = naturalearth(&quot;admin_0_countries&quot;, 10)</span></span></code></pre></div><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span>GO.area(all_adm0)</span></span></code></pre></div><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span>AccurateArithmetic.sum_oro(GO.area.(all_adm0.geometry))</span></span></code></pre></div><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span>AccurateArithmetic.sum_kbn(GO.area.(all_adm0.geometry))</span></span></code></pre></div><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span>GI.Polygon.(GO.flatten(Union{GI.LineStringTrait, GI.LinearRingTrait}, all_adm0) |&gt; collect .|&gt; x -&gt; [x]) .|&gt; GO.signed_area |&gt; sum</span></span></code></pre></div><div class="language-@example vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">@example</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span>GI.Polygon.(GO.flatten(Union{GI.LineStringTrait, GI.LinearRingTrait}, all_adm0) |&gt; collect .|&gt; x -&gt; [x]) .|&gt; GO.signed_area |&gt; sum_oro</span></span></code></pre></div><p>@example accurate GI.Polygon.(GO.flatten(Union{GI.LineStringTrait, GI.LinearRingTrait}, all_adm0) |&gt; collect .|&gt; x -&gt; [x]) .|&gt; GO.signed_area |&gt; sum_kbn \`\`\`</p>`,10),i=[c];function l(p,o,r,u,d,m){return t(),e("div",null,i)}const _=a(n,[["render",l]]);export{g as __pageData,_ as default};
