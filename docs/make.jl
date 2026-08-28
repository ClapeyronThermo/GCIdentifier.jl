using Documenter, DocumenterVitepress, GCIdentifier

makedocs(
    sitename = "GCIdentifier.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/ClapeyronThermo/GCIdentifier.jl",
    ),
    warnonly = Documenter.except(),
    authors = "Pierre J. Walker and Andrés Riedemann.",
    pages = [
        "Home" => "index.md",
        "Group Assignment" => "group_search.md",
        "Finding Missing Groups" => "missing_groups.md",
        "Custom Groups" => "custom_groups.md",
        "API" => "api.md",
    ],
)

DocumenterVitepress.deploydocs(
    repo = "github.com/ClapeyronThermo/GCIdentifier.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
