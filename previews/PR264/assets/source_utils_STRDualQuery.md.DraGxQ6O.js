import{_ as a,c as n,o as i,az as e}from"./chunks/framework.DgVVAcWa.js";const c=JSON.parse('{"title":"","description":"","frontmatter":{},"headers":[],"relativePath":"source/utils/STRDualQuery.md","filePath":"source/utils/STRDualQuery.md","lastUpdated":null}'),p={name:"source/utils/STRDualQuery.md"};function l(t,s,h,r,k,d){return i(),n("div",null,s[0]||(s[0]=[e(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    STRDualQuery</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">A module for performing dual-tree traversals on STRtrees to find potentially overlapping geometry pairs.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">The main entry point is \`maybe_overlapping_geoms_and_query_lists_in_order\`.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">module</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> STRDualQuery</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> SortTileRecursiveTree</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> SortTileRecursiveTree</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> STRtree, STRNode, STRLeafNode</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeoInterface</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Extents</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GeoInterface </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">as</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GI</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;helper function to get the extent of any STR node, since leaf nodes don&#39;t store global extent.&quot;</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    maybe_overlapping_geoms_and_query_lists_in_order(tree_a::STRtree, tree_b::STRtree, edges_a::Vector{&lt;: GI.Line}, edges_b::Vector{&lt;: GI.Line})</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Performs an efficient dual-tree traversal to find potentially overlapping geometry pairs.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Returns a vector of pairs, where each pair contains an index from tree_a and a sorted vector</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">of indices from tree_b that might overlap.</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">The result looks like this:</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">[</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    a1 =&gt; [b1, b2, b3],</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    a2 =&gt; [b4, b5],</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    a3 =&gt; [b6, b7, b8, b9],</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    ...</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">]</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">\`\`\`</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">in which the overlap map is sorted by the tree_a indices, and within each group, the tree_b indices are sorted.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> maybe_overlapping_geoms_and_query_lists_in_order</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(tree_a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">STRtree</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, tree_b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">STRtree</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, edges_a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Vector{&lt;: GI.Line}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, edges_b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Vector{&lt;: GI.Line}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Use DefaultDict to automatically create empty vectors for new keys</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    overlap_map </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Dict{Int, Vector{Int}}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">()</span></span></code></pre></div><p>Start the recursive traversal from the root nodes</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    _dual_tree_traverse!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(overlap_map, tree_a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rootnode, tree_b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">rootnode, edges_a, edges_b)</span></span></code></pre></div><p>Convert to the required output format and sort</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    result </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [(k, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">sort!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(v)) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (k, v) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> pairs</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(overlap_map)]</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">    sort!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(result, by</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">first)  </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># Sort by tree_a indices</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    return</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> result</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    _dual_tree_traverse!(</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        overlap_map::Dict{Int,Vector{Int}},</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        node_a::Union{STRNode,STRLeafNode}, node_b::Union{STRNode,STRLeafNode},</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">        edges_a::Vector{&lt;: GI.Line}, edges_b::Vector{&lt;: GI.Line}</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">    )</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">Recursive helper function that performs the dual-tree traversal and stores results in \`overlap_map\`.</span></span>
<span class="line"><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;&quot;&quot;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">function</span><span style="--shiki-light:#6F42C1;--shiki-dark:#B392F0;"> _dual_tree_traverse!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    overlap_map</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Dict{Int,Vector{Int}}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    node_a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Union{STRNode,STRLeafNode}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    node_b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Union{STRNode,STRLeafNode}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    edges_a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Vector{&lt;: GI.Line}</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    edges_b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">::</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Vector{&lt;: GI.Line}</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><p>Early exit if bounding boxes don&#39;t overlap</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> !</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Extents</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">intersects</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">extent</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(node_a), GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">extent</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(node_b))</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        return</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span></code></pre></div><p>Case 1: Both nodes are leaves</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_a </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">isa</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> STRLeafNode </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&amp;&amp;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_b </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">isa</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> STRLeafNode</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> idx_a </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">indices</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">            dict_vec </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> get!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Int[], overlap_map, idx_a)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">            for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> idx_b </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">indices</span></span></code></pre></div><p>Final extent rejection, this is cheaper than the allocation for <code>push!</code></p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">                if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Extents</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">intersects</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">extent</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(edges_a[idx_a]), GI</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">extent</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(edges_b[idx_b]))</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">                    push!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(dict_vec, idx_b)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">                end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">            end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        return</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span></code></pre></div><p>Case 2: node_a is a leaf, node_b is internal</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_a </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">isa</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> STRLeafNode</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> child_b </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">children</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">            _dual_tree_traverse!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(overlap_map, node_a, child_b, edges_a, edges_b)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        return</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span></code></pre></div><p>Case 3: node_b is a leaf, node_a is internal</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    if</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_b </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">isa</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> STRLeafNode</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> child_a </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">children</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">            _dual_tree_traverse!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(overlap_map, child_a, node_b, edges_a, edges_b)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        return</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span></code></pre></div><p>Case 4: Both nodes are internal</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> child_a </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_a</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">children</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> child_b </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> node_b</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">children</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">            _dual_tree_traverse!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(overlap_map, child_a, child_b, edges_a, edges_b)</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">    end</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">export</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> maybe_overlapping_geoms_and_query_lists_in_order</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">end</span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"> # module</span></span></code></pre></div><p>The code to use this is: elseif accelerator isa DoubleSTRtree # If both of the polygons are quite large, # then we do a dual-tree traversal of the STRtrees # and find all potentially overlapping edges. # This is kind of like an adjacency list of a graph or a sparse matrix that we&#39;re constructing.</p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span># First, we materialize an edge list for each polygon as a vector of \`GI.Line\` objects.</span></span>
<span class="line"><span>ext_a = GI.extent(poly_a)</span></span>
<span class="line"><span>ext_b = GI.extent(poly_b)</span></span>
<span class="line"><span>edges_a, indices_a = to_edgelist(ext_b, poly_a, T) # ::Vector{GI.Line} # with precalculated extents</span></span>
<span class="line"><span>edges_b, indices_b = to_edgelist(ext_a, poly_b, T) # ::Vector{GI.Line} # with precalculated extents</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Now, construct STRtrees from the edge lists.</span></span>
<span class="line"><span># TODO: we can optimize the strtrees by passing in only edges that reside within the area of extent overlap</span></span>
<span class="line"><span># between poly_a and poly_b.</span></span>
<span class="line"><span># this would greatly help for e.g. coverage transactions</span></span>
<span class="line"><span># BUT this needs fixes in SortTileRecursiveTree.jl to allow you to pass in a vector of indices</span></span>
<span class="line"><span># as well as extents / geometries.</span></span>
<span class="line"><span>tree_a = STRtree(edges_a)</span></span>
<span class="line"><span>tree_b = STRtree(edges_b)</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Now we perform a dual-tree traversal to find</span></span>
<span class="line"><span># all potentially overlapping pairs of edges.</span></span>
<span class="line"><span># The \`STRDualQuery\` module is defined in this folder,</span></span>
<span class="line"><span># in \`strtree_dual_query.jl\`.</span></span>
<span class="line"><span>index_list = STRDualQuery.maybe_overlapping_geoms_and_query_lists_in_order(</span></span>
<span class="line"><span>    tree_a, tree_b, edges_a, edges_b</span></span>
<span class="line"><span>)</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Track the last index in \`poly_a\` that we&#39;ve processed.</span></span>
<span class="line"><span>last_a_idx_orig = indices_a[1]</span></span>
<span class="line"><span></span></span>
<span class="line"><span># Iterate over the indices in \`poly_a\` that may potentially overlap with any edge in \`poly_b\`.</span></span>
<span class="line"><span>for (a_idx, b_idxs) in index_list</span></span>
<span class="line"><span>    a_idx_orig = indices_a[a_idx]</span></span>
<span class="line"><span>    # If we&#39;ve skipped any indices in \`poly_a\`,</span></span>
<span class="line"><span>    # then we need to process those indices.</span></span>
<span class="line"><span>    # But we know that they will never intersect, so we can skip calling \`f_intersect\`.</span></span>
<span class="line"><span>    if a_idx_orig &gt; (last_a_idx_orig + 1)</span></span>
<span class="line"><span>        for i in (last_a_idx_orig + 1:a_idx_orig - 1)</span></span>
<span class="line"><span>            # TODO explain this check</span></span>
<span class="line"><span>            # a1, a2 = edges_a[i].geom</span></span>
<span class="line"><span>            # a1 == a2 &amp;&amp; continue</span></span>
<span class="line"><span>            point_a = GI.getpoint(poly_a, i)</span></span>
<span class="line"><span>            f_on_each_a(point_a, i)</span></span>
<span class="line"><span>            f_after_each_a(point_a, i)</span></span>
<span class="line"><span>        end</span></span>
<span class="line"><span>    end</span></span>
<span class="line"><span>    last_a_idx_orig = a_idx_orig</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    # Get the edge in \`poly_a\` that we&#39;re currently processing.</span></span>
<span class="line"><span>    aedge = edges_a[a_idx]</span></span>
<span class="line"><span>    aedge_extent = GI.extent(aedge)</span></span>
<span class="line"><span>    a1t, a2t = aedge.geom</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    # Call \`f_a\` for the start point of the edge.</span></span>
<span class="line"><span>    f_on_each_a(a1t, a_idx_orig)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    for b_idx in b_idxs</span></span>
<span class="line"><span>        b_idx_orig = indices_b[b_idx]</span></span>
<span class="line"><span>        bedge = edges_b[b_idx]</span></span>
<span class="line"><span>        b1t, b2t = bedge.geom</span></span>
<span class="line"><span>        b1t == b2t &amp;&amp; continue</span></span>
<span class="line"><span></span></span>
<span class="line"><span>        if GI.Extents.intersects(aedge_extent, GI.extent(bedge))</span></span>
<span class="line"><span>            LoopStateMachine.@controlflow f_on_each_maybe_intersect(</span></span>
<span class="line"><span>                ((a1t, a2t), a_idx_orig), ((b1t, b2t), b_idx_orig)</span></span>
<span class="line"><span>            )</span></span>
<span class="line"><span>        end</span></span>
<span class="line"><span>    end</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    # Call \`f_after_each_a\` for the start point of the edge - postprocess once</span></span>
<span class="line"><span>    # we&#39;re done with tracing.</span></span>
<span class="line"><span>    f_after_each_a(a1t, a_idx_orig)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    last_a_idx_orig = a_idx_orig</span></span>
<span class="line"><span></span></span>
<span class="line"><span>end</span></span>
<span class="line"><span></span></span>
<span class="line"><span># println(&quot;DOUBLE TREE TRAVERSAL&quot;)</span></span></code></pre></div><p>end</p><p>using Test using GeometryOps using GeometryOps.STRDualQuery using SortTileRecursiveTree using StaticArrays using GI.Extents</p><p>@testset &quot;STRDualQuery&quot; begin @testset &quot;Basic overlapping rectangles&quot; begin # Create two sets of rectangles represented by their corner points # Tree A: Tree B: # [0,0]–[1,1] [0.5,0.5]–<a href="./overlaps with A1">1.5,1.5</a> # [2,2]–[3,3] [2.5,2.5]–<a href="./overlaps with A2">3.5,3.5</a></p><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>    edges_a = [</span></span>
<span class="line"><span>        ((0.0, 0.0), (1.0, 1.0)),  # A1</span></span>
<span class="line"><span>        ((2.0, 2.0), (3.0, 3.0))   # A2</span></span>
<span class="line"><span>    ]</span></span>
<span class="line"><span>    edges_b = [</span></span>
<span class="line"><span>        ((0.5, 0.5), (1.5, 1.5)),  # B1</span></span>
<span class="line"><span>        ((2.5, 2.5), (3.5, 3.5))   # B2</span></span>
<span class="line"><span>    ]</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    # Convert edges to STRtree format</span></span>
<span class="line"><span>    tree_a = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_a])</span></span>
<span class="line"><span>    tree_b = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_b])</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    result = maybe_overlapping_geoms_and_query_lists_in_order(tree_a, tree_b)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    # Check results</span></span>
<span class="line"><span>    @test length(result) == 2</span></span>
<span class="line"><span>    @test result[1][1] == 1  # First edge from tree_a</span></span>
<span class="line"><span>    @test result[1][2] == [1]  # Overlaps with first edge from tree_b</span></span>
<span class="line"><span>    @test result[2][1] == 2  # Second edge from tree_a</span></span>
<span class="line"><span>    @test result[2][2] == [2]  # Overlaps with second edge from tree_b</span></span>
<span class="line"><span>end</span></span>
<span class="line"><span></span></span>
<span class="line"><span>@testset &quot;Non-overlapping geometries&quot; begin</span></span>
<span class="line"><span>    edges_a = [((0.0, 0.0), (1.0, 1.0))]</span></span>
<span class="line"><span>    edges_b = [((10.0, 10.0), (11.0, 11.0))]</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    tree_a = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_a])</span></span>
<span class="line"><span>    tree_b = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_b])</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    result = maybe_overlapping_geoms_and_query_lists_in_order(tree_a, tree_b)</span></span>
<span class="line"><span>    @test isempty(result)</span></span>
<span class="line"><span>end</span></span>
<span class="line"><span></span></span>
<span class="line"><span>@testset &quot;Multiple overlaps&quot; begin</span></span>
<span class="line"><span>    # Create a scenario where one edge overlaps with multiple others</span></span>
<span class="line"><span>    edges_a = [((0.0, 0.0), (2.0, 2.0))]  # One long diagonal line</span></span>
<span class="line"><span>    edges_b = [</span></span>
<span class="line"><span>        ((0.5, 0.5), (1.0, 1.0)),  # B1 overlaps</span></span>
<span class="line"><span>        ((1.0, 1.0), (1.5, 1.5)),  # B2 overlaps</span></span>
<span class="line"><span>        ((1.5, 1.5), (2.0, 2.0)),  # B3 overlaps</span></span>
<span class="line"><span>        ((3.0, 3.0), (4.0, 4.0))   # B4 doesn&#39;t overlap</span></span>
<span class="line"><span>    ]</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    tree_a = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_a])</span></span>
<span class="line"><span>    tree_b = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_b])</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    result = maybe_overlapping_geoms_and_query_lists_in_order(tree_a, tree_b)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    @test length(result) == 1</span></span>
<span class="line"><span>    @test result[1][1] == 1</span></span>
<span class="line"><span>    @test result[1][2] == [1, 2, 3]  # Should find first three edges from tree_b</span></span>
<span class="line"><span>end</span></span>
<span class="line"><span></span></span>
<span class="line"><span>@testset &quot;Empty trees&quot; begin</span></span>
<span class="line"><span>    empty_tree = STRtree(GI.Line{2, Float64}[])</span></span>
<span class="line"><span>    edges_a = [((0.0, 0.0), (1.0, 1.0))]</span></span>
<span class="line"><span>    non_empty_tree = STRtree([GI.Line(SVector{2}(p1, p2)) for (p1, p2) in edges_a])</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    # Test empty tree with non-empty tree</span></span>
<span class="line"><span>    result1 = maybe_overlapping_geoms_and_query_lists_in_order(empty_tree, non_empty_tree)</span></span>
<span class="line"><span>    @test isempty(result1)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    # Test non-empty tree with empty tree</span></span>
<span class="line"><span>    result2 = maybe_overlapping_geoms_and_query_lists_in_order(non_empty_tree, empty_tree)</span></span>
<span class="line"><span>    @test isempty(result2)</span></span>
<span class="line"><span></span></span>
<span class="line"><span>    # Test empty tree with empty tree</span></span>
<span class="line"><span>    result3 = maybe_overlapping_geoms_and_query_lists_in_order(empty_tree, empty_tree)</span></span>
<span class="line"><span>    @test isempty(result3)</span></span>
<span class="line"><span>end</span></span></code></pre></div><p>end</p><hr><p><em>This page was generated using <a href="https://github.com/fredrikekre/Literate.jl" target="_blank" rel="noreferrer">Literate.jl</a>.</em></p>`,28)]))}const g=a(p,[["render",l]]);export{c as __pageData,g as default};
