"""
    AnnotatedProfile(omics, obs, obsindex)
    AnnotatedProfile(op, omicsname, obs, obsindex)

# Arguments

- `omics::Dict{Symbol,OmicsProfile}`: A collection of omics profile with their names in keys.
- `op::OmicsProfile`: Single omics profile.
- `omicsname::Symbol`: The name of given `OmicsProfile`.
- `obs::DataFrame`: The dataframe contains meta information for observations.
- `obsindex::Symbol`: The index of dataframe `obs`.
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
    Profile(countmat, omicsname, var, obs; varindex=:genesymbol, obsindex::Symbol=:barcode,
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
"""
function Profile(countmat::AbstractMatrix, omicsname, var::DataFrame, obs::DataFrame;
        varindex::Symbol=:genesymbol, obsindex::Symbol=:barcode, T=float(eltype(countmat)))
    p = OmicsProfile(countmat, var, varindex; T=T)
    return AnnotatedProfile(p, omicsname, obs, obsindex)
end

obsnames(p::AnnotatedProfile) = names(p.obs)
varnames(p::AnnotatedProfile, omicsname::Symbol) = varnames(p.omics[omicsname])
layernames(p::AnnotatedProfile, omicsname::Symbol) = layernames(p.omics[omicsname])
omicsnames(p::AnnotatedProfile) = keys(p.omics)
countmatrix(p::AnnotatedProfile, omicsname::Symbol) = countmatrix(p.omics[omicsname])

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

# Base.copy(p::AnnotatedProfile)

getobsindex(p::AnnotatedProfile) = p.obsindex[]

function setobsindex!(p::AnnotatedProfile, index::Symbol)
    not_unique = nonunique(p.obs, index)
    any(not_unique) && throw(ArgumentError("not unique obs index at rows: $(findall(not_unique))."))
    p.obsindex[] = index
    return p
end

getpipeline(p::AnnotatedProfile, i::Symbol) = p.pipeline[i]

function setpipeline!(p::AnnotatedProfile, x, i::Symbol)
    p.pipeline[i] = x
    return p
end
