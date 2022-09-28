Base.filter(x::Pair{Symbol,<:Any}, p::AbstractProfile) = filter!(x, copy(p))

function Base.filter!(x::Pair{Symbol,T}, p::OmicsProfile) where {T}
    col, f = x
    @assert string(col) in names(p.var)
    sel = f.(p.var[:, col])
    filter!(x, p.var)
    isempty(p.varm) || filter_mappings!(p, sel)
    isempty(p.layers) || filter_layers!(p, sel)
    return p
end

function Base.filter!(x::Pair{Symbol,T}, p::AnnotatedProfile) where {T}
    col, f = x
    if col in names(p.obs)
        sel = f.(p.obs[:, col])
        filter!(x, p.obs)
        isempty(p.obsm) || filter_mappings!(p, sel)
        isempty(p.omics) || filter_layers!(p; obsidx=sel)
    else
        for omicsname in omicsnames(p)
            op = p.omics[omicsname]
            col in names(op.var) && filter!(x, op)
        end
    end
    return p
end

function filter_layers!(p::OmicsProfile, varidx=(:), obsidx=(:))
    for k in layernames(p)
        p.layers[k] = p.layers[k][varidx, obsidx]
    end
    return p
end

function filter_layers!(p::AnnotatedProfile; varidx=(:), obsidx=(:))
    for omicsname in omicsnames(p)
        filter_layers!(p.omics[omicsname], varidx, obsidx)
    end
    return p
end

function filter_mappings!(p::OmicsProfile, varidx=(:))
    for k in keys(p.varm)
        p.varm[k] = p.varm[k][varidx, :]
    end
    return p
end

function filter_mappings!(p::AnnotatedProfile, obsidx=(:))
    for k in keys(p.obsm)
        p.obsm[k] = p.obsm[k][:, obsidx]
    end
    return p
end

# function filter_graphs!(p::OmicsProfile, varidx=(:))

# end

# function filter_graphs!(p::AnnotatedProfile, obsidx=(:))

# end
