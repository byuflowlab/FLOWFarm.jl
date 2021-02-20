using Documenter, FlowFarm

makedocs(;
    modules=[FlowFarm],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "Tutorial.md",
        "How-To" => "How_to.md",
        "Explanation" => "Explanation.md",
        "Reference" => "Reference.md"
    ],
    repo="https://github.com/byuflowlab/FlowFarm.jl/blob/{commit}{path}#L{line}",
    sitename="FlowFarm.jl",
    authors="PJ Stanley <stanley_andrewpj@yahoo.com>",
    assets=String[],
)

deploydocs(
    repo = "github.com/byuflowlab/FlowFarm.jl.git",
)