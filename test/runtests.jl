using OmicsProfiles
using SnowyOwl
using SparseArrays
using SHA
using CodecZlib
using DataFrames
using CUDA
using Test

const TEST_PATH = @__DIR__

ENV["DATADEPS_ALWAYS_ACCEPT"] = true

cuda_tests = [
    "device",
]

tests = [
    "io",
    "omicsprofile",
    "annotatedprofile",
    "filter",
]

if CUDA.functional()
    CUDA.allowscalar(false)
    append!(tests, cuda_tests)
else
    @warn "CUDA unavailable, not testing GPU support"
end

@testset "OmicsProfiles.jl" begin
    for t in tests
        include("$(t).jl")
    end
end
