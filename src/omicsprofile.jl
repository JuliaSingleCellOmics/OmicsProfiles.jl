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

varnames(p::OmicsProfile) = names(getfield(p, :var))
layernames(p::OmicsProfile) = keys(getfield(p, :layers))
pipelinenames(p::OmicsProfile) = keys(getfield(p, :pipeline))

nrow(p::OmicsProfile) = size(p.count, 1)
ncol(p::OmicsProfile) = size(p.count, 2)
nvar(p::OmicsProfile) = nrow(p.var)

function Base.show(io::IO, p::OmicsProfile)
    print(io, "OmicsProfile (nvar = ", nvar(p), "):")
    isempty(p.var) || print(io, "\n    var: ", join(varnames(p), ", "))
    isempty(p.varm) || print(io, "\n    varm: ", join(keys(p.varm), ", "))
    isempty(p.vargraphs) || print(io, "\n    vargraphs: ", join(keys(p.vargraphs), ", "))
    isempty(p.layers) || print(io, "\n    layers: ", join(layernames(p), ", "))
    isempty(p.pipeline) || print(io, "\n    pipeline: ", join(keys(p.pipeline), " => "))
end

==(p1::OmicsProfile, p2::OmicsProfile) = p1.var == p2.var && getvarindex(p1) == getvarindex(p2) &&
    p1.varm == p2.varm && p1.vargraphs == p2.vargraphs && p1.layers == p2.layers &&
    p1.pipeline == p2.pipeline

Base.copy(p::OmicsProfile) = OmicsProfile(copy(p.var), p.varindex, copy(p.varm),
    copy(p.vargraphs), copy(p.layers), copy(p.pipeline))

Base.maximum(p::OmicsProfile) = maximum(p.count)
Base.minimum(p::OmicsProfile) = minimum(p.count)

Base.size(p::OmicsProfile) = size(p.count)
Base.axes(p::OmicsProfile) = Base.axes(p.count)

getvarindex(p::OmicsProfile) = getfield(p, :varindex)[]

function setvarindex!(p::OmicsProfile, index::Symbol)
    not_unique = nonunique(p.var, index)
    any(not_unique) && throw(ArgumentError("not unique var index at rows: $(findall(not_unique))."))
    getfield(p, :varindex)[] = index
    return p
end

getlayer(p::OmicsProfile, i::Symbol) = AxisArray(getfield(p, :layers)[i]; row=p.var[!, p.varindex])

function setlayer!(p::OmicsProfile, i::Symbol, x)
    @assert size(x, 1) == nvar(p)
    getfield(p, :layers)[i] = x
    return p
end

getpipeline(p::OmicsProfile, i::Symbol) = getfield(p, :pipeline)[i]

function setpipeline!(p::OmicsProfile, i::Symbol, x)
    getfield(p, :pipeline)[i] = x
    return p
end

function Base.getproperty(p::OmicsProfile, name::Symbol)
    if name == :varindex
        return getvarindex(p)
    elseif name in keys(getfield(p, :varm))
        return getfield(p, :varm)[name]
    elseif name in keys(getfield(p, :vargraphs))
        return getfield(p, :vargraphs)[name]
    elseif name in layernames(p)
        return getlayer(p, name)
    elseif name in pipelinenames(p)
        return getpipeline(p, name)
    else
        return getfield(p, name)
    end
end

function Base.setproperty!(p::OmicsProfile, name::Symbol, x)
    if name == :varindex
        return setvarindex!(p, x)
    elseif name in keys(getfield(p, :varm))
        getfield(p, :varm)[name] = x
        return p
    elseif name in keys(getfield(p, :vargraphs))
        getfield(p, :vargraphs)[name] = x
        return p
    elseif name in layernames(p)
        return setlayer!(p, name, x)
    elseif name in pipelinenames(p)
        return setpipeline!(p, name, x)
    else
        return setfield!(p, name, x)
    end
end

function Base.propertynames(p::OmicsProfile)
    props = keys(getfield(p, :varm)) ∪ keys(getfield(p, :vargraphs)) ∪ layernames(p) ∪
        pipelinenames(p)
    return (:var, :varindex, :varm, :vargraphs, :layers, :pipeline, props...)
end

function Base.getindex(p::OmicsProfile, inds...)
    new_prof = OmicsProfile(getindex(p.count, inds[1], :),
                getindex(p.var, inds[1], :), getvarindex(p))
    for name in layernames(p)
        l = getindex(getlayer(p, name), inds[1], :)
        setlayer!(new_prof, name, parent(l))
    end
    for name in pipelinenames(p)
        pl = copy(getpipeline(p, name))
        setpipeline!(new_prof, name, pl)
    end
    return new_prof
end

# Base.view(p::OmicsProfile, inds...)

# Base.merge!

# merge
