mutable struct AnnotatedProfile <: AbstractProfile
    omics::Dict{Symbol,OmicsProfile}
    obs::DataFrame
    obsindex::Symbol
    obsm
    obsgraphs::Dict{Symbol,AbstractGraph}
    pipeline::OrderedDict{Symbol,Dict}
end

function AnnotatedProfile(omics::Dict{Symbol,OmicsProfile}, obs::DataFrame, obsindex::Symbol,
        )
    @assert all(nrow(obs) == ncol(p) for p in values(omics)) "all omics data should have the same number of observations with `obs`."
    obsm
    obsgraphs = Dict{Symbol,AbstractGraph}()
    pipeline = OrderedDict{Symbol,Dict}()
    return AnnotatedProfile(omics, obs, obsindex, obsm, obsgraphs, pipeline)
end

function AnnotatedProfile(p::OmicsProfile, omicsname, obs::DataFrame, obsindex::Symbol)
    omics = Dict{Symbol,OmicsProfile}(Symbol(omicsname) => p)
    return AnnotatedProfile(omics, obs, obsindex)
end

function Profile(countmat, omicsname, var::DataFrame, obs::DataFrame;
        varindex::Symbol=:genesymbol, obsindex::Symbol=:barcode)
    p = OmicsProfile(countmat, var, varindex)
    return AnnotatedProfile(p, omicsname, obs, obsindex)
end

obsnames(p::AnnotatedProfile) = names(p.obs)
varnames(p::AnnotatedProfile, kind::Symbol) = varnames(p.omics[kind])
layernames(p::AnnotatedProfile, kind::Symbol) = layernames(p.omics[kind])
omicsnames(p::AnnotatedProfile) = keys(p.omics)
countmatrix(p::AnnotatedProfile, kind::Symbol) = countmatrix(p.omics[kind])

nvar(p::AnnotatedProfile, kind::Symbol) = nvar(p.omics[kind])
nobs(p::AnnotatedProfile) = nrow(p.obs)

function Base.show(io::IO, p::AnnotatedProfile)
    println(io, "AnnotatedProfile: n_obs = ", nobs(p))
    isempty(p.obs) || println(io, "    obs: ", join(obsnames(p), ", "))
    isempty(p.pipeline) || println(io, "    pipeline: ", join(keys(p.pipeline), " => "))
    for (kind, omic) in p.omics
        print(io, "$kind => ", omic)
    end
end

function setpipeline!(p::AnnotatedProfile, x, i::Symbol)
    p.pipeline[i] = x
    return p
end

getpipeline(p::AnnotatedProfile, i::Symbol) = p.pipeline[i]
