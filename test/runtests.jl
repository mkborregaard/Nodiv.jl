using Aqua
using Nodiv
using Phylo
using Random
using RecipesBase
using SpatialEcology
using Statistics
using Test

include("toydata.jl")

@testset "Nodiv.jl" begin
    include("metrics.jl")
    include("analysis.jl")
    include("distances.jl")
    include("grouping.jl")
    include("recipes.jl")
    include("deprecated.jl")
    include("aqua.jl")
end
