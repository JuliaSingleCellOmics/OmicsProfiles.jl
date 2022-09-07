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
    OmicsProfile,
    obsnames,
    varnames,
    layernames,
    nrow,
    ncol,
    nvar,
    nobs,
    get_gene_expr


include("mtx.jl")
include("omicsprofile.jl")
include("annotatedprofile.jl")

end
