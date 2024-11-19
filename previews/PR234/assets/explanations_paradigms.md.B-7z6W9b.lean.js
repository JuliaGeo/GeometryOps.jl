import{_ as a,c as t,a5 as o,o as i}from"./chunks/framework.CzXzEn_U.js";const u=JSON.parse('{"title":"Paradigms","description":"","frontmatter":{},"headers":[],"relativePath":"explanations/paradigms.md","filePath":"explanations/paradigms.md","lastUpdated":null}'),r={name:"explanations/paradigms.md"};function s(n,e,d,p,c,l){return i(),t("div",null,e[0]||(e[0]=[o('<h1 id="paradigms" tabindex="-1">Paradigms <a class="header-anchor" href="#paradigms" aria-label="Permalink to &quot;Paradigms&quot;">​</a></h1><p>GeometryOps exposes functions like <code>apply</code> and <code>applyreduce</code>, as well as the <code>fix</code> and <code>prepare</code> APIs, that represent <em>paradigms</em> of programming, by which we mean the ability to program in a certain way, and in so doing, fit neatly into the tools we&#39;ve built without needing to re-implement the wheel.</p><p>Below, we&#39;ll describe some of the foundational paradigms of GeometryOps, and why you should care!</p><h2 id="apply" tabindex="-1"><code>apply</code> <a class="header-anchor" href="#apply" aria-label="Permalink to &quot;`apply` {#apply}&quot;">​</a></h2><p>The <code>apply</code> function allows you to decompose a given collection of geometries down to a certain level, operate on it, and reconstruct it back to the same nested form as the original. In general, its invocation is:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">apply</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(f, trait</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Trait</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, geom)</span></span></code></pre></div><p>Functionally, it&#39;s similar to <code>map</code> in the way you apply it to geometries - except that you tell it at which level it should stop, by passing a <code>trait</code> to it.</p><p><code>apply</code> will start by decomposing the geometry, feature, featurecollection, iterable, or table that you pass to it, and stop when it encounters a geometry for which <code>GI.trait(geom) isa Trait</code>. This encompasses unions of traits especially, but beware that any geometry which is not explicitly handled, and hits <code>GI.PointTrait</code>, will cause an error.</p><p><code>apply</code> is unlike <code>map</code> in that it returns reconstructed geometries, instead of the raw output of the function. If you want a purely map-like behaviour, like calculating the length of each linestring in your feature collection, then call <code>GO.flatten(f, trait, geom)</code>, which will decompose each geometry to the given <code>trait</code> and apply <code>f</code> to it, returning the decomposition as a flattened vector.</p><h3 id="applyreduce" tabindex="-1"><code>applyreduce</code> <a class="header-anchor" href="#applyreduce" aria-label="Permalink to &quot;`applyreduce` {#applyreduce}&quot;">​</a></h3><p><code>applyreduce</code> is like the previous <code>map</code>-based approach that we mentioned, except that it <code>reduce</code>s the result of <code>f</code> by <code>op</code>. Note that <code>applyreduce</code> does not guarantee associativity, so it&#39;s best to have <code>typeof(init) == returntype(op)</code>.</p><h2 id="fix-and-prepare" tabindex="-1"><code>fix</code> and <code>prepare</code> <a class="header-anchor" href="#fix-and-prepare" aria-label="Permalink to &quot;`fix` and `prepare` {#fix-and-prepare}&quot;">​</a></h2><p>The <code>fix</code> and <code>prepare</code> paradigms are different from <code>apply</code>, though they are built on top of it. They involve the use of structs as &quot;actions&quot;, where a constructed object indicates an action that should be taken. A trait like interface prescribes the level (polygon, linestring, point, etc) at which each action should be applied.</p><p>In general, the idea here is to be able to invoke several actions efficiently and simultaneously, for example when correcting invalid geometries, or instantiating a <code>Prepared</code> geometry with several preparations (sorted edge lists, rtrees, monotone chains, etc.)</p>',14)]))}const m=a(r,[["render",s]]);export{u as __pageData,m as default};
