module OmicsProfiles

using LinearAlgebra
using SparseArrays

using DataFrames
using DataStructures
using Graphs, SimpleWeightedGraphs
using TranscodingStreams

import Base: ==
import DataFrames: nrow, ncol

const PROJECT_PATH = dirname(@__DIR__)
const DenseOrSparse = Union{Matrix,SparseMatrixCSC}
const FEATURE_COLS = [:featurekey, :featurename, :featuretype, :chromosome, :featurestart, :featureend, :isgene, :genus_species]

export
    # io
    read_mtx,
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
    geneexpr,

    # annotatedProfile
    AnnotatedProfile,
    Profile,
    obsnames,
    omicsnames,
    nobs,
    getobsindex,
    setobsindex!


include("mtx.jl")
include("io.jl")
include("omicsprofile.jl")
include("annotatedprofile.jl")

end
