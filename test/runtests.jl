using OmicsProfiles
using SparseArrays
using SHA
using CodecZlib
using DataFrames
using Test

tests = [
    "mtx",
    "omicsprofile",
    "annotatedprofile",
]

@testset "OmicsProfiles.jl" begin
    for t in tests
        include("$(t).jl")
    end
end
