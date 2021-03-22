push!(LOAD_PATH,"../src/")
using Documenter, SchoolChoice

makedocs(modules = [SchoolChoice], sitename = "SchoolChoice.jl")

deploydocs(repo = "github.com/pedrovergaramerino/SchoolChoice.jl.git", devbranch = "main")
