using GeometryOps
using Documenter
using Literate
using Makie, CairoMakie
CairoMakie.activate!(px_per_unit = 2, type = "png") # TODO: make this svg

DocMeta.setdocmeta!(GeometryOps, :DocTestSetup, :(using GeometryOps); recursive=true)

source_path = joinpath(dirname(@__DIR__), "src")
output_path = joinpath(@__DIR__, "src", "source")
mkpath(output_path)

literate_pages = String[]

withenv("JULIA_DEBUG" => "Literate") do # allow Literate debug output to escape to the terminal!
    for (root_path, dirs, files) in walkdir(source_path)
        output_dir = joinpath(output_path, relpath(root_path, source_path))
        for file in files
            # convert the source file to Markdown
            Literate.markdown(joinpath(root_path, file), output_dir; documenter = false)
            # TODO: make this respect nesting somehow!
            push!(literate_pages, joinpath("source", relpath(root_path, source_path), splitext(file)[1] * ".md"))
        end
    end
end

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
        "Source code" => literate_pages
    ]
)

deploydocs(;
    repo="github.com/asinghvi17/GeometryOps.jl",
    devbranch="main",
    push_preview = true,
)
