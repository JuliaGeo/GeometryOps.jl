using GeometryOps
using Documenter

DocMeta.setdocmeta!(GeometryOps, :DocTestSetup, :(using GeometryOps); recursive=true)

makedocs(;
    modules=[GeometryOps],
    authors="Anshul Singhvi <anshulsinghvi@gmail.com> and contributors",
    repo="https://github.com/asinghvi17/GeometryOps.jl/blob/{commit}{path}#{line}",
    sitename="GeometryOps.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://asinghvi17.github.io/GeometryOps.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/asinghvi17/GeometryOps.jl",
    devbranch="main",
)
