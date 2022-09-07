abstract type AbstractProfile end

mutable struct OmicsProfile <: AbstractProfile
    var::DataFrame
    varindex::Symbol
    varm
    vargraphs::Dict{Symbol,AbstractGraph}
    layers::Dict{Symbol,DenseOrSparse}
    pipeline::OrderedDict{Symbol,Dict}
end

function OmicsProfile(countmat::M, var::DataFrame, varindex::Symbol; T=float(eltype(countmat))) where {M<:AbstractMatrix}
    @assert nrow(var) == size(countmat, 1) "count matrix should have the same number of variables with `var`."
    not_unique = nonunique(var, varindex)
    any(not_unique) || throw(ArgumentError("not unique var index at rows: $(findall(not_unique))."))
    varm
    vargraphs = Dict{Symbol,AbstractGraph}()
    S = Union{Matrix{T},SparseMatrixCSC{T,UInt32}}
    layers = Dict{Symbol,S}(:count => T.(countmat))
    pipeline = OrderedDict{Symbol,Dict}()
    return OmicsProfile(copy(var), varindex, varm, vargraphs, layers, pipeline)
end

varnames(p::OmicsProfile) = names(p.var)
layernames(p::OmicsProfile) = keys(p.layers)
countmatrix(p::OmicsProfile) = p.layers[:count]

nrow(p::OmicsProfile) = size(countmatrix(p), 1)
ncol(p::OmicsProfile) = size(countmatrix(p), 2)
nvar(p::OmicsProfile) = nrow(p.var)

function Base.show(io::IO, p::OmicsProfile)
    println(io, "OmicsProfile: n_var = ", nvar(p))
    isempty(p.var) || println(io, "    var: ", join(varnames(p), ", "))
    isempty(p.layers) || println(io, "    layers: ", join(layernames(p), ", "))
    isempty(p.pipeline) || println(io, "    pipeline: ", join(keys(p.pipeline), " => "))
end

function Base.copy(p::OmicsProfile)
    return OmicsProfile(copy(p.var), p.varindex, copy(p.varm), copy(p.vargraphs),
        copy(p.layers), copy(p.pipeline))
end

Base.maximum(p::OmicsProfile) = maximum(countmatrix(p))
Base.minimum(p::OmicsProfile) = minimum(countmatrix(p))

Base.size(p::OmicsProfile) = size(countmatrix(p))
Base.axes(p::OmicsProfile) = axes(countmatrix(p))

function Base.getindex(p::OmicsProfile, inds...)
    p_ = OmicsProfile(getindex(countmatrix(p), inds[1], inds[2]),
                 getindex(p.var, inds[1], :),
                 getindex(p.obs, inds[2], :))
    for (k, v) in p.layers
        if size(v) == size(p_.data)
            p_.layers[k] = getindex(v, inds[1], inds[2])
        end
    end
    p_.pipeline = copy(p.pipeline)
    return p_
end

Base.setproperty!(p::OmicsProfile, name::Symbol, x) = setproperty!(p, Val(name), x)
Base.setproperty!(p::OmicsProfile, ::Val{S}, x) where S = setfield!(p, S, x)

function Base.setproperty!(p::OmicsProfile, ::Val{:obs}, x)
    @assert nrow(x) == size(p.data, 2)
    setfield!(p, :obs, x)
end

function Base.setproperty!(p::OmicsProfile, ::Val{:var}, x)
    @assert nrow(x) == size(p.data, 1)
    setfield!(p, :var, x)
end

Base.filter(x::Pair{Symbol,T}, p::OmicsProfile) where {T} = filter!(x, copy(p))

function Base.filter!(x::Pair{Symbol,T}, p::OmicsProfile) where {T}
    col, f = x
    sel = f.(p.var[:,col])
    filter!(x, p.var)
    filter_layers!(p, var_idx=sel)
    p.data = p.data[sel, :]
    return p
end

function filter_layers!(p::OmicsProfile; var_idx=(:), obs_idx=(:))
    for k in keys(p.layers)
        if size(p.layers[k]) == size(p.data)
            p.layers[k] = p.layers[k][var_idx, obs_idx]
        end
    end
    return p
end

function get_gene_expr(p::OmicsProfile, gene_name, ::Nothing=nothing)
    idx = collect(1:nrow(p))[p.var.index .== gene_name]
    return view(p.data, idx, :)
end

function get_gene_expr(p::OmicsProfile, gene_name, layer::Symbol)
    idx = collect(1:nrow(p))[p.var.index .== gene_name]
    return view(p.layers[layer], idx, :)
end

function setpipeline!(p::OmicsProfile, x, i::Symbol)
    p.pipeline[i] = x
    return p
end

getpipeline(p::OmicsProfile, i::Symbol) = p.pipeline[i]

# function set_index!(p::OmicsProfile; obs="", var="")
#     if obs != ""
#         # make obs unique
#     end

#     if var != ""
#         # make var unique
#     end
# end
