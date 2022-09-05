mutable struct OmicsProfile{T<:AbstractMatrix,S}
    data::T
    var::DataFrame
    obs::DataFrame
    layers::Dict{Symbol,S}
    pipeline::OrderedDict
end

function OmicsProfile(data::M, var::DataFrame, obs::DataFrame; T=float(eltype(data))) where {M<:AbstractMatrix}
    @assert (nrow(var), nrow(obs)) == size(data)
    S = Union{Matrix{T},SparseMatrixCSC{T,UInt32}}
    data = T.(data)
    layers = Dict{Symbol,S}()
    pipeline = OrderedDict{Symbol,Dict}()
    OmicsProfile{typeof(data),S}(data, copy(var), copy(obs), layers, pipeline)
end

obsnames(p::OmicsProfile) = names(p.obs)
varnames(p::OmicsProfile) = names(p.var)
layernames(p::OmicsProfile) = keys(p.layers)
nrow(p::OmicsProfile) = size(p.data, 1)
ncol(p::OmicsProfile) = size(p.data, 2)
nvar(p::OmicsProfile) = size(p.data, 1)
nobs(p::OmicsProfile) = size(p.data, 2)

function Base.show(io::IO, p::OmicsProfile)
    println(io, "OmicsProfile(n_var × n_obs = ", nrow(p), " × ", ncol(p), ")")
    isempty(p.obs) || println(io, "    obs: ", join(obsnames(p), ", "))
    isempty(p.var) || println(io, "    var: ", join(varnames(p), ", "))
    isempty(p.layers) || println(io, "    layers: ", join(layernames(p), ", "))
    isempty(p.pipeline) || println(io, "    pipeline: ", join(keys(p.pipeline), ", "))
end

Base.copy(p::OmicsProfile) = OmicsProfile(copy(p.data), copy(p.var), copy(p.obs), copy(p.layers), copy(p.pipeline))

Base.maximum(p::OmicsProfile) = maximum(p.data)
Base.minimum(p::OmicsProfile) = minimum(p.data)

Base.size(p::OmicsProfile) = size(p.data)
Base.axes(p::OmicsProfile) = axes(p.data)

function Base.getindex(p::OmicsProfile, inds...)
    p_ = OmicsProfile(getindex(p.data, inds[1], inds[2]),
                 getindex(p.var, inds[1], :),
                 getindex(p.obs, inds[2], :))
    for (k, v) in p.layers
        if size(v) == size(p_.data)
            p_.layers[k] = getindex(v, inds[1], inds[2])
        end
    end
    p_.pipeline = copy(p.pipeline)
    p_
end

Base.setproperty!(prof::OmicsProfile, name::Symbol, x) = setproperty!(prof, Val(name), x)
Base.setproperty!(prof::OmicsProfile, ::Val{S}, x) where S = setfield!(prof, S, x)

function Base.setproperty!(prof::OmicsProfile, ::Val{:obs}, x)
    @assert nrow(x) == size(prof.data, 2)
    setfield!(prof, :obs, x)
end

function Base.setproperty!(prof::OmicsProfile, ::Val{:var}, x)
    @assert nrow(x) == size(prof.data, 1)
    setfield!(prof, :var, x)
end

Base.filter(x::Pair{Symbol,T}, prof::OmicsProfile) where {T} = filter!(x, copy(prof))

function Base.filter!(x::Pair{Symbol,T}, prof::OmicsProfile) where {T}
    col, f = x
    sel = f.(prof.var[:,col])
    filter!(x, prof.var)
    filter_layers!(prof, var_idx=sel)
    prof.data = prof.data[sel, :]
    prof
end

function filter_layers!(prof::OmicsProfile; var_idx=(:), obs_idx=(:))
    for k in keys(prof.layers)
        if size(prof.layers[k]) == size(prof.data)
            prof.layers[k] = prof.layers[k][var_idx, obs_idx]
        end
    end
    prof
end

function get_gene_expr(prof::OmicsProfile, gene_name, ::Nothing=nothing)
    idx = collect(1:nrow(prof))[prof.var.index .== gene_name]
    return view(prof.data, idx, :)
end

function get_gene_expr(prof::OmicsProfile, gene_name, layer::Symbol)
    idx = collect(1:nrow(prof))[prof.var.index .== gene_name]
    return view(prof.layers[layer], idx, :)
end

# function set_index!(p::OmicsProfile; obs="", var="")
#     if obs != ""
#         # make obs unique
#     end

#     if var != ""
#         # make var unique
#     end
# end
