using Documenter
using CTDirectShooting

makedocs(
    sitename = "CTDirectShooting.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction" => "index.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/CTDirectShooting.jl.git",
    devbranch = "main"
)
