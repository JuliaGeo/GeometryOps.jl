import{_ as t,c as a,m as s,a as i,a7 as h,o as n}from"./chunks/framework.DHnBKsCQ.js";const k="/GeometryOps.jl/v0.1.2/assets/piujvlk.CNHhYFWU.png",R=JSON.parse('{"title":"Barycentric coordinates","description":"","frontmatter":{},"headers":[],"relativePath":"source/methods/barycentric.md","filePath":"source/methods/barycentric.md","lastUpdated":null}'),l={name:"source/methods/barycentric.md"},p=h("",4),e={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},E={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"10.692ex",height:"2.262ex",role:"img",focusable:"false",viewBox:"0 -750 4726 1000","aria-hidden":"true"},r=h("",1),d=[r],g=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mo",{stretchy:"false"},"("),s("msub",null,[s("mi",null,"λ"),s("mn",null,"1")]),s("mo",null,","),s("msub",null,[s("mi",null,"λ"),s("mn",null,"2")]),s("mo",null,","),s("msub",null,[s("mi",null,"λ"),s("mn",null,"3")]),s("mo",{stretchy:"false"},")")])],-1),y={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},F={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.025ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.357ex",height:"1.025ex",role:"img",focusable:"false",viewBox:"0 -442 600 453","aria-hidden":"true"},o=s("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[s("g",{"data-mml-node":"math"},[s("g",{"data-mml-node":"mi"},[s("path",{"data-c":"1D45B",d:"M21 287Q22 293 24 303T36 341T56 388T89 425T135 442Q171 442 195 424T225 390T231 369Q231 367 232 367L243 378Q304 442 382 442Q436 442 469 415T503 336T465 179T427 52Q427 26 444 26Q450 26 453 27Q482 32 505 65T540 145Q542 153 560 153Q580 153 580 145Q580 144 576 130Q568 101 554 73T508 17T439 -10Q392 -10 371 17T350 73Q350 92 386 193T423 345Q423 404 379 404H374Q288 404 229 303L222 291L189 157Q156 26 151 16Q138 -11 108 -11Q95 -11 87 -5T76 7T74 17Q74 30 112 180T152 343Q153 348 153 366Q153 405 129 405Q91 405 66 305Q60 285 60 284Q58 278 41 278H27Q21 284 21 287Z",style:{"stroke-width":"3"}})])])],-1),c=[o],C=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"n")])],-1),B={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},A={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.025ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.357ex",height:"1.025ex",role:"img",focusable:"false",viewBox:"0 -442 600 453","aria-hidden":"true"},D=s("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[s("g",{"data-mml-node":"math"},[s("g",{"data-mml-node":"mi"},[s("path",{"data-c":"1D45B",d:"M21 287Q22 293 24 303T36 341T56 388T89 425T135 442Q171 442 195 424T225 390T231 369Q231 367 232 367L243 378Q304 442 382 442Q436 442 469 415T503 336T465 179T427 52Q427 26 444 26Q450 26 453 27Q482 32 505 65T540 145Q542 153 560 153Q580 153 580 145Q580 144 576 130Q568 101 554 73T508 17T439 -10Q392 -10 371 17T350 73Q350 92 386 193T423 345Q423 404 379 404H374Q288 404 229 303L222 291L189 157Q156 26 151 16Q138 -11 108 -11Q95 -11 87 -5T76 7T74 17Q74 30 112 180T152 343Q153 348 153 366Q153 405 129 405Q91 405 66 305Q60 285 60 284Q58 278 41 278H27Q21 284 21 287Z",style:{"stroke-width":"3"}})])])],-1),u=[D],T=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"n")])],-1),m={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},b={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.566ex"},xmlns:"http://www.w3.org/2000/svg",width:"14.876ex",height:"2.262ex",role:"img",focusable:"false",viewBox:"0 -750 6575.4 1000","aria-hidden":"true"},Q=h("",1),_=[Q],v=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mo",{stretchy:"false"},"("),s("msub",null,[s("mi",null,"λ"),s("mn",null,"1")]),s("mo",null,","),s("msub",null,[s("mi",null,"λ"),s("mn",null,"2")]),s("mo",null,","),s("mo",null,"."),s("mo",null,"."),s("mo",null,"."),s("mo",null,","),s("msub",null,[s("mi",null,"λ"),s("mi",null,"n")]),s("mo",{stretchy:"false"},")")])],-1),w=h("",35);function f(x,V,P,M,q,N){return n(),a("div",null,[p,s("p",null,[i("In the case of a triangle, barycentric coordinates are a set of three numbers "),s("mjx-container",e,[(n(),a("svg",E,d)),g]),i(", each associated with a vertex of the triangle. Any point within the triangle can be expressed as a weighted average of the vertices, where the weights are the barycentric coordinates. The weights sum to 1, and each is non-negative.")]),s("p",null,[i("For a polygon with "),s("mjx-container",y,[(n(),a("svg",F,c)),C]),i(" vertices, generalized barycentric coordinates are a set of "),s("mjx-container",B,[(n(),a("svg",A,u)),T]),i(" numbers "),s("mjx-container",m,[(n(),a("svg",b,_)),v]),i(", each associated with a vertex of the polygon. Any point within the polygon can be expressed as a weighted average of the vertices, where the weights are the generalized barycentric coordinates.")]),w])}const G=t(l,[["render",f]]);export{R as __pageData,G as default};
