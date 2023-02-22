module OmicsProfiles

using LinearAlgebra
using SparseArrays
using Mmap

using AxisArrays
using CUDA, Adapt
using DataFrames, CSV
using DataStructures
using Graphs, SimpleWeightedGraphs
using MatrixMarket
using TranscodingStreams, CodecZlib

import Base: ==
import DataFrames: nrow, ncol

const PROJECT_PATH = dirname(@__DIR__)
const DenseOrSparse = Union{Matrix,SparseMatrixCSC,CuMatrix}

export
    # io
    read_10x,

    # device
    tocpu,
    togpu,

    # omicsprofile
    AbstractProfile,
    OmicsProfile,
    nrow,
    ncol,
    nvar,
    varnames,
    getvarindex,
    setvarindex!,
    layernames,
    getlayer,
    setlayer!,
    getpipeline,
    setpipeline!,

    # annotatedProfile
    AnnotatedProfile,
    Profile,
    obsnames,
    omicsnames,
    nobs,
    getobsindex,
    setobsindex!


include("io.jl")
include("omicsprofile.jl")
include("annotatedprofile.jl")
include("device.jl")
include("filter.jl")

end
