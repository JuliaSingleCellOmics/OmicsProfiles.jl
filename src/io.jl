const HCA_GENES_HEADER = [:featurekey, :featurename, :featuretype, :chromosome,
    :featurestart, :featureend, :isgene, :genus_species]

read_mtx(filename::AbstractString) = MatrixMarket.mmread(filename)

read_csv(file, delim::Char, header) = CSV.File(file; delim=delim, header=header) |> DataFrame

function read_csv(filename::AbstractString; delim=nothing, header::Vector{Symbol}=Symbol[])
    header = isempty(header) ? 1 : header
    delim = predict_delim(delim, filename)

    if endswith(filename, ".gz")
        file = transcode(GzipDecompressor, Mmap.mmap(filename))
    else
        file = filename
    end
    return read_csv(file, delim, header)
end

predict_delim(delim::Char, filename::AbstractString) = delim

function predict_delim(::Nothing, filename::AbstractString)
    if occursin(".csv", filename)
        return ','
    elseif occursin(".tsv", filename)
        return '\t'
    else
        @warn "extensions of filename $filename is not identified, assign delim = ','"
        return ','
    end
end

find_file(path::AbstractString, infix::AbstractString, ext::AbstractString) =
    find_file(path, infix, [ext])

find_file(path::AbstractString, infix::AbstractString, exts::AbstractVector) =
    find_file(path, [infix], exts)

find_file(path::AbstractString, infixes::AbstractVector, ext::AbstractString) =
    find_file(path, infixes, [ext])

function find_file(path::AbstractString, infixes::AbstractVector, exts::AbstractVector)
    for filename in readdir(path)
        for infix in infixes
            for ext in exts
                if occursin("$infix.$ext", filename)
                    return joinpath(path, filename)
                end
            end
        end
    end

    error("attempt to read [$(join(infixes, '/'))].[$(join(exts, '/'))].[gz] file in $path, but file not found.")
end

"""
    read_10x(path; make_unique=true, omicsname=:RNA, varnames=[:ensembleid, :genesymbol],
        obsnames=[:barcode], varindex=:genesymbol, obsindex=:barcode)

Reads 10x dataset from `path` and returns an `AnnotatedProfile`.

# Arguments

- `path::AbstractString`: The path to 10x dataset folder.
- `make_unique::Bool`: To make `varindex` unique or not if given `varindex` is not unique.
- `omicsname::Symbol`: The name of given sequencing data.
- `varnames::AbstractVector`: The header or column names of feature file.
- `obsnames::AbstractVector`: The header or column names of barcode file.
- `varindex::Symbol`: The index of dataframe `var`.
- `obsindex::Symbol`: The index of dataframe `obs`.

# Examples

```jldoctest
julia> using SnowyOwl

julia> prof = read_10x(SnowyOwl.Dataset.pbmc3k_folder(), varindex=:ensembleid)
AnnotatedProfile (nobs = 2700):
    obs: barcode
RNA => OmicsProfile (nvar = 32738):
    var: ensembleid, genesymbol
    layers: count
```
"""
function read_10x(path::AbstractString; make_unique::Bool=true, omicsname::Symbol=:RNA,
    varnames::AbstractVector=[:ensembleid, :genesymbol], obsnames::AbstractVector=[:barcode],
    varindex::Symbol=:genesymbol, obsindex::Symbol=:barcode)
    varfile = find_file(path, ["genes", "features"], ["csv", "tsv"])
    obsfile = find_file(path, "barcodes", ["csv", "tsv"])
    exprfile = find_file(path, "matrix", "mtx")
    var = read_csv(varfile, header=varnames)
    obs = read_csv(obsfile, header=obsnames)
    X = SparseMatrixCSC{Float64,UInt32}(read_mtx(exprfile))
    make_unique && (var = make_index_unique(var, varindex))
    return Profile(X, omicsname, var, obs, varindex=varindex, obsindex=obsindex)
end

function make_index_unique(var, varindex::Symbol)
    # @show nrow(unique(var, [varindex]))
    var
end

# function read_10x_h5()

# end

# function read_h5ad()

# end

# read_hdf
# read_loom
