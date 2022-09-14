using OmicsProfiles
using SparseArrays
using SHA
using Test

tests = [
    "mtx",
]

@testset "OmicsProfiles.jl" begin
    for t in tests
        include("$(t).jl")
    end
end
