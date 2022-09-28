"""
    filter(col => func, p)

Returns a copy of profile `p` containing data from `p` for which `func` returns `true`.

# Arguments

- `col`: Column name to filter.
- `func`: A predicate function which returns `true`/`false`.
- `p::AbstractProfile`: A profile to filter on.

# Examples

```jldoctest
julia> using OmicsProfiles, DataFrames

julia> ngenes, ncells = (100, 500)
(100, 500)

julia> X = rand(0:10, ngenes, ncells);

julia> var = DataFrame(genesymbol=1:ngenes, A=rand(ngenes));

julia> prof = OmicsProfile(X, var, :genesymbol)
OmicsProfile (nvar = 100):
    var: genesymbol, A
    layers: count

julia> prof.varm[:a] = rand(ngenes, 50);

julia> prof2 = filter(:A => x -> x > 0.5, prof)
OmicsProfile (nvar = 46):
    var: genesymbol, A
    varm: a
    layers: count

julia> obs = DataFrame(barcode=1:ncells, B=rand(ncells));

julia> ap = AnnotatedProfile(prof, :RNA, obs, :barcode)
AnnotatedProfile (nobs = 500):
    obs: barcode, B
RNA => OmicsProfile (nvar = 100):
    var: genesymbol, A
    varm: a
    layers: count

julia> ap2 = filter(:B => x -> x > 0.5, ap)
AnnotatedProfile (nobs = 253):
    obs: barcode, B
RNA => OmicsProfile (nvar = 100):
    var: genesymbol, A
    varm: a
    layers: count
```
"""
Base.filter(x::Pair{Symbol,<:Any}, p::AbstractProfile) = filter!(x, copy(p))

function Base.filter!(x::Pair{Symbol,T}, p::OmicsProfile) where {T}
    col, f = x
    @assert string(col) in names(p.var)
    sel = f.(p.var[!, col])
    filter!(x, p.var)
    isempty(p.varm) || filter_mappings!(p, sel)
    isempty(p.layers) || filter_layers!(p, sel)
    return p
end

function Base.filter!(x::Pair{Symbol,T}, p::AnnotatedProfile) where {T}
    col, f = x
    if string(col) in names(p.obs)
        sel = f.(p.obs[!, col])
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

function filter_graphs!(p::OmicsProfile, varidx=(:))
    return p
end

function filter_graphs!(p::AnnotatedProfile, obsidx=(:))
    return p
end
