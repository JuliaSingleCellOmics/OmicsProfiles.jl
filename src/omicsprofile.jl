abstract type AbstractProfile end

"""
    OmicsProfile(countmat, var, varindex; T=float(eltype(countmat)))

# Arguments

- `countmat`: The count matrix for given omics and it will be set to key `:count`.
- `var::DataFrame`: The dataframe contains meta information for features or variables.
- `varindex::Symbol`: The index of dataframe `var`.
- `T`: Element type of `countmat`.

# Examples

```jldoctest
julia> using OmicsProfiles, DataFrames

julia> r, c = (100, 500)
(100, 500)

julia> data = rand(0:100, r, c);

julia> var = DataFrame(index=1:r, C=rand(r), D=rand(r));

julia> prof = OmicsProfile(data, var, :index)
OmicsProfile (nvar = 100):
    var: index, C, D
    layers: count
```
"""
struct OmicsProfile <: AbstractProfile
    var::DataFrame
    varindex::Ref{Symbol}
    varm::Dict{Symbol,AbstractMatrix}
    vargraphs::Dict{Symbol,AbstractGraph}
    layers::Dict{Symbol,DenseOrSparse}
    pipeline::OrderedDict{Symbol,Dict}
end

function OmicsProfile(countmat::M, var::DataFrame, varindex::Symbol; T=float(eltype(countmat))) where {M<:AbstractMatrix}
    @assert nrow(var) == size(countmat, 1) "count matrix should have the same number of variables with `var`."
    not_unique = nonunique(var, varindex)
    any(not_unique) && throw(ArgumentError("not unique var index at rows: $(findall(not_unique))."))
    varm = Dict{Symbol,AbstractMatrix}()
    vargraphs = Dict{Symbol,AbstractGraph}()
    S = Union{Matrix{T},SparseMatrixCSC{T,UInt32}}
    layers = Dict{Symbol,S}(:count => T.(countmat))
    pipeline = OrderedDict{Symbol,Dict}()
    return OmicsProfile(copy(var), Ref(varindex), varm, vargraphs, layers, pipeline)
end

varnames(p::OmicsProfile) = names(p.var)
layernames(p::OmicsProfile) = keys(p.layers)
countmatrix(p::OmicsProfile) = p.layers[:count]

nrow(p::OmicsProfile) = size(countmatrix(p), 1)
ncol(p::OmicsProfile) = size(countmatrix(p), 2)
nvar(p::OmicsProfile) = nrow(p.var)

function Base.show(io::IO, p::OmicsProfile)
    print(io, "OmicsProfile (nvar = ", nvar(p), "):")
    isempty(p.var) || print(io, "\n    var: ", join(varnames(p), ", "))
    isempty(p.varm) || print(io, "\n    varm: ", join(keys(p.varm), ", "))
    isempty(p.vargraphs) || print(io, "\n    vargraphs: ", join(keys(p.vargraphs), ", "))
    isempty(p.layers) || print(io, "\n    layers: ", join(layernames(p), ", "))
    isempty(p.pipeline) || print(io, "\n    pipeline: ", join(keys(p.pipeline), " => "))
end

function Base.copy(p::OmicsProfile)
    return OmicsProfile(copy(p.var), p.varindex, copy(p.varm), copy(p.vargraphs),
        copy(p.layers), copy(p.pipeline))
end

Base.maximum(p::OmicsProfile) = maximum(countmatrix(p))
Base.minimum(p::OmicsProfile) = minimum(countmatrix(p))

Base.size(p::OmicsProfile) = size(countmatrix(p))
Base.axes(p::OmicsProfile) = axes(countmatrix(p))

getvarindex(p::OmicsProfile) = p.varindex[]

function setvarindex!(p::OmicsProfile, index::Symbol)
    not_unique = nonunique(p.var, index)
    any(not_unique) && throw(ArgumentError("not unique var index at rows: $(findall(not_unique))."))
    p.varindex[] = index
    return p
end

getlayer(p::OmicsProfile, i::Symbol) = p.layers[i]

function setlayer!(p::OmicsProfile, x, i::Symbol)
    @assert size(x, 1) == nvar(p)
    p.layers[i] = x
    return p
end

getpipeline(p::OmicsProfile, i::Symbol) = p.pipeline[i]

function setpipeline!(p::OmicsProfile, x, i::Symbol)
    p.pipeline[i] = x
    return p
end

function geneexpr(p::OmicsProfile, gene_name, layer::Symbol=:count)
    varindex = getvarindex(p)
    idx = findall(p.var[!, varindex] .== gene_name)
    return view(getlayer(p, layer), idx, :)
end

function Base.getindex(p::OmicsProfile, inds...)
    new_prof = OmicsProfile(getindex(countmatrix(p), inds[1], inds[2]),
                 getindex(p.var, inds[1], :),
                 getindex(p.obs, inds[2], :))
    for (k, v) in p.layers
        new_prof.layers[k] = getindex(v, inds[1], inds[2])
    end
    new_prof.pipeline = copy(p.pipeline)
    return new_prof
end

# Base.merge!

# merge
