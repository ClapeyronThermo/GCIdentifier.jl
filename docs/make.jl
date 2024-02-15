using Documenter,GCIdentifier

makedocs(sitename = "GCIdentifier.jl",
format = Documenter.HTML(
    # Use clean URLs, unless built as a "local" build
    canonical = "https://ClapeyronThermo.github.io/GCIdentifier.jl/",
),
warnonly = Documenter.except(),
    authors = "Pierre J. Walker and Andrés Riedemann.",
    pages = [
        "Home" => "index.md",
        "Background" => "background.md",
        "Basic Usage" => "basic_usage.md",        
        "API" => "api.md"]

        deploydocs(;
    repo="github.com/ClapeyronThermo/GCIdentifier.jl.git",
)
