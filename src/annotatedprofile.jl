"""
    AnnotatedProfile(omics, obs, obsindex)
    AnnotatedProfile(op, omicsname, obs, obsindex)

A data container preserves multi-omics data for single cell sequencing analysis. Multi-omics
`omics` in dictionary or, in separation, single omics `op` with its name `omicsname` and
observations metadata `obs` are retained in the container. It records data and progression
mainly about observations.

# Arguments

- `omics::Dict{Symbol,OmicsProfile}`: A collection of omics profile with their names in keys.
- `op::OmicsProfile`: Single omics profile.
- `omicsname::Symbol`: The name of given `OmicsProfile`.
- `obs::DataFrame`: The dataframe contains meta information for observations.
- `obsindex::Symbol`: The index of dataframe `obs`.

# Examples

```jldoctest
julia> ngenes, ncells = (100, 500)
(100, 500)

julia> data = rand(0:10, ngenes, ncells);

julia> var = DataFrame(genesymbol=1:ngenes);

julia> obs = DataFrame(barcode=1:ncells);

julia> omic = OmicsProfile(data, var, :genesymbol);

julia> ap = AnnotatedProfile(omic, :RNA, obs, :barcode)
AnnotatedProfile (nobs = 500):
    obs: barcode
RNA => OmicsProfile (nvar = 100):
    var: genesymbol
    layers: count
```

See also [`Profile`](@ref) for routine use and [`OmicsProfile`](@ref) for single omics data container.
"""
struct AnnotatedProfile <: AbstractProfile
    omics::Dict{Symbol,OmicsProfile}
    obs::DataFrame
    obsindex::Ref{Symbol}
    obsm::Dict{Symbol,AbstractMatrix}
    obsgraphs::Dict{Symbol,AbstractGraph}
    pipeline::OrderedDict{Symbol,Dict}
end

function AnnotatedProfile(omics::Dict{Symbol,OmicsProfile}, obs::DataFrame, obsindex::Symbol)
    @assert all(nrow(obs) == ncol(p) for p in values(omics)) "all omics data should have the same number of observations with `obs`."
    obsm = Dict{Symbol,AbstractMatrix}()
    obsgraphs = Dict{Symbol,AbstractGraph}()
    pipeline = OrderedDict{Symbol,Dict}()
    return AnnotatedProfile(omics, obs, Ref(obsindex), obsm, obsgraphs, pipeline)
end

function AnnotatedProfile(p::OmicsProfile, omicsname, obs::DataFrame, obsindex::Symbol)
    omics = Dict{Symbol,OmicsProfile}(Symbol(omicsname) => p)
    return AnnotatedProfile(omics, obs, obsindex)
end

"""
    Profile(countmat, omicsname, var, obs; varindex=:genesymbol, obsindex=:barcode,
            T=float(eltype(countmat)))

Constructor for establishing `AnnotatedProfile` with a `OmicsProfile` inside.

# Arguments

- `countmat`: The count matrix for given omics and it will be set to key `:count`.
- `omicsname::Symbol`: The name of given `OmicsProfile`.
- `var::DataFrame`: The dataframe contains meta information for features or variables.
- `obs::DataFrame`: The dataframe contains meta information for observations.
- `varindex::Symbol`: The index of dataframe `var`.
- `obsindex::Symbol`: The index of dataframe `obs`.
- `T`: Element type of `countmat`.

# Examples

```jldoctest
julia> ngenes, ncells = (100, 500)
(100, 500)

julia> data = rand(0:10, ngenes, ncells);

julia> var = DataFrame(genesymbol=1:ngenes);

julia> obs = DataFrame(barcode=1:ncells);

julia> prof = Profile(data, :RNA, var, obs, varindex=:genesymbol, obsindex=:barcode)
AnnotatedProfile (nobs = 500):
    obs: barcode
RNA => OmicsProfile (nvar = 100):
    var: genesymbol
    layers: count
```

See also [`AnnotatedProfile`](@ref) for multi-omics data container and [`OmicsProfile`](@ref)
for single omics data container.
"""
function Profile(countmat::AbstractMatrix, omicsname, var::DataFrame, obs::DataFrame;
        varindex::Symbol=:genesymbol, obsindex::Symbol=:barcode, T=float(eltype(countmat)))
    p = OmicsProfile(countmat, var, varindex; T=T)
    return AnnotatedProfile(p, omicsname, obs, obsindex)
end

omicsnames(p::AnnotatedProfile) = keys(getfield(p, :omics))
obsnames(p::AnnotatedProfile) = names(p.obs)
varnames(p::AnnotatedProfile, omicsname::Symbol) = varnames(p.omics[omicsname])
layernames(p::AnnotatedProfile, omicsname::Symbol) = layernames(p.omics[omicsname])
pipelinenames(p::AnnotatedProfile) = keys(getfield(p, :pipeline))

nvar(p::AnnotatedProfile, omicsname::Symbol) = nvar(p.omics[omicsname])
nobs(p::AnnotatedProfile) = nrow(p.obs)

function Base.show(io::IO, p::AnnotatedProfile)
    print(io, "AnnotatedProfile (nobs = ", nobs(p), "):")
    isempty(p.obs) || print(io, "\n    obs: ", join(obsnames(p), ", "))
    isempty(p.obsm) || print(io, "\n    obsm: ", join(keys(p.obsm), ", "))
    isempty(p.obsgraphs) || print(io, "\n    obsgraphs: ", join(keys(p.obsgraphs), ", "))
    isempty(p.pipeline) || print(io, "\n    pipeline: ", join(keys(p.pipeline), " => "))
    for (kind, omic) in p.omics
        print(io, "\n$kind => ", omic)
    end
end

==(p1::AnnotatedProfile, p2::AnnotatedProfile) = p1.omics == p2.omics &&
    p1.obs == p2.obs && getobsindex(p1) == getobsindex(p2) && p1.obsm == p2.obsm &&
    p1.obsgraphs == p2.obsgraphs && p1.pipeline == p2.pipeline

Base.copy(p::AnnotatedProfile) = AnnotatedProfile(copy(p.omics), copy(p.obs),
    Ref(getobsindex(p)), copy(p.obsm), copy(p.obsgraphs), copy(p.pipeline))

getobsindex(p::AnnotatedProfile) = getfield(p, :obsindex)[]

function setobsindex!(p::AnnotatedProfile, index::Symbol)
    not_unique = nonunique(p.obs, index)
    any(not_unique) && throw(ArgumentError("not unique obs index at rows: $(findall(not_unique))."))
    getfield(p, :obsindex)[] = index
    return p
end

getpipeline(p::AnnotatedProfile, i::Symbol) = getfield(p, :pipeline)[i]

function setpipeline!(p::AnnotatedProfile, i::Symbol, x)
    getfield(p, :pipeline)[i] = x
    return p
end

function Base.getproperty(p::AnnotatedProfile, name::Symbol)
    if name == :obsindex
        return getobsindex(p)
    elseif name in omicsnames(p)
        return getfield(p, :omics)[name]
    elseif name in keys(getfield(p, :obsm))
        return getfield(p, :obsm)[name]
    elseif name in keys(getfield(p, :obsgraphs))
        return getfield(p, :obsgraphs)[name]
    elseif name in keys(getfield(p, :pipeline))
        return getpipeline(p, name)
    else
        return getfield(p, name)
    end
end

function Base.setproperty!(p::AnnotatedProfile, name::Symbol, x)
    if name == :obsindex
        return setobsindex!(p, x)
    elseif name in keys(getfield(p, :obsm))
        getfield(p, :obsm)[name] = x
        return p
    elseif name in keys(getfield(p, :obsgraphs))
        getfield(p, :obsgraphs)[name] = x
        return p
    elseif name in pipelinenames(p)
        return setpipeline!(p, name, x)
    else
        return setfield!(p, name, x)
    end
end

function Base.propertynames(p::AnnotatedProfile)
    props = keys(getfield(p, :obsm)) ∪ keys(getfield(p, :obsgraphs)) ∪ pipelinenames(p)
    return (:omics, :obs, :obsindex, :obsm, :obsgraphs, :pipeline, props...)
end
