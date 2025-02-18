using Documenter, FLOWFarm

DocMeta.setdocmeta!(FLOWFarm, :DocTestSetup, :(using FLOWFarm); recursive=true)

makedocs(;
    modules=[FLOWFarm],
    format=Documenter.HTML(repolink="https://github.com/byuflowlab/FLOWFarm.jl", size_threshold=nothing, size_threshold_warn=nothing),
    pages=[
        "Intro" => "index.md",
        "Quick Start" => "Tutorial.md",
        "Wind Farm Struct" => "WindFarmStruct.md",
        "How-to Guide" => "How_to.md",
        "Theory" => "Explanation.md",
        "Reference" => "Reference.md"
    ],
    repo="https://github.com/byuflowlab/FLOWFarm.jl/blob/{commit}{path}#L{line}",
    sitename="FLOWFarm.jl",
    authors="Jared J. Thomas <jaredthomas68@gmail.com>, PJ Stanley <stanley_andrewpj@yahoo.com>",
    doctest=true
)

deploydocs(
    repo = "github.com/byuflowlab/FLOWFarm.jl.git",
    versions = nothing,
    devbranch = "develop"
)
