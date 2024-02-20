# Adding the package src to the load path
push!(LOAD_PATH,"../src/")

using Documenter, LinearResponse

makedocs(sitename = "LinearResponse.jl",
         pages=[
                "Home" => "index.md",
                "Installation" => "installation.md",
                "Quickstart" => "quickstart.md",
                "Functions" => "functions.md",
                "IO" => "io.md"
               ],
         format = Documenter.HTML(prettyurls=false))

deploydocs(repo="github.com/JuliaStellarDynamics/LinearResponse.jl",devbranch="documentation")