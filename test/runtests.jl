using OmicsProfiles
using SparseArrays
using SHA
using CodecZlib
using DataFrames
using Test

const TEST_PATH = @__DIR__

tests = [
    "mtx",
    "io",
    "omicsprofile",
    "annotatedprofile",
    "filter"
]

@testset "OmicsProfiles.jl" begin
    for t in tests
        include("$(t).jl")
    end
end
