module OmicsProfiles

using LinearAlgebra
using SparseArrays
using Mmap

using AxisArrays
using DataFrames, CSV
using DataStructures
using Graphs, SimpleWeightedGraphs
using MatrixMarket
using TranscodingStreams, CodecZlib

import Base: ==
import DataFrames: nrow, ncol

const PROJECT_PATH = dirname(@__DIR__)
const DenseOrSparse = Union{Matrix,SparseMatrixCSC}
const FEATURE_COLS = [:featurekey, :featurename, :featuretype, :chromosome, :featurestart, :featureend, :isgene, :genus_species]

export
    # io
    read_mtx,
    read_csv,
    read_features,
    read_barcodes,
    read_genes,
    read_cells,

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
include("filter.jl")

end
