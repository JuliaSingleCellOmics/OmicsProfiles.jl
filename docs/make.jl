using OmicsProfiles
using Documenter

DocMeta.setdocmeta!(OmicsProfiles, :DocTestSetup, :(using OmicsProfiles); recursive=true)

makedocs(;
    modules=[OmicsProfiles],
    authors="Yueh-Hua Tu",
    repo="https://github.com/JuliaSingleCellOmics/OmicsProfiles.jl/blob/{commit}{path}#{line}",
    sitename="OmicsProfiles.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliasinglecellomics.github.io/OmicsProfiles.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSingleCellOmics/OmicsProfiles.jl",
    devbranch="main",
)
