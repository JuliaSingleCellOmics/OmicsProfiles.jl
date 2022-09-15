module OmicsProfiles

using LinearAlgebra
using SparseArrays

using DataFrames
using DataStructures
using Graphs, SimpleWeightedGraphs
using TranscodingStreams

import DataFrames: nrow, ncol

const PROJECT_PATH = dirname(@__DIR__)
const DenseOrSparse = Union{Matrix,SparseMatrixCSC}
const DEFAULT_FEATURE_COLS = [:ensembleid, :genesymbol, :type]
const DEFAULT_BARCODE_COLS = [:barcode]
const FEATURE_COLS = [:featurekey, :featurename, :featuretype, :chromosome, :featurestart, :featureend, :isgene, :genus_species]

export
    # omicsprofile
    AbstractProfile,
    OmicsProfile,
    countmatrix,
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
    get_gene_expr,

    # annotatedProfile
    AnnotatedProfile,
    Profile,
    obsnames,
    omicsnames,
    nobs,
    getobsindex,
    setobsindex!


include("mtx.jl")
include("omicsprofile.jl")
include("annotatedprofile.jl")

end
