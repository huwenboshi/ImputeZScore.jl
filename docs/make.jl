using Documenter, TraitSimulation
 
makedocs(modules=[TraitSimulation], doctest=true)
 
deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/huwenboshi/TraitSimulation.jl.git",
    julia  = "0.5", osname = "linux")
