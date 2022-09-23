read_mtx(filename::String) = mmread(filename)

read_csv(file, delim::Char, header::Vector{Symbol}) =
    CSV.File(file; delim=delim, header=header) |> DataFrame

function read_csv(filename::String; delim=nothing, header::Vector{Symbol}=Symbol[])
    header = isempty(header) ? 1 : header

    if isnothing(delim)
        delim = ','
    elseif endswith(filename, ".csv") || endswith(filename, r"\.csv\.\w+")
        delim = ','
    elseif endswith(filename, ".tsv") || endswith(filename, r"\.tsv\.\w+")
        delim = '\t'
    end

    if endswith(filename, ".gz")
        file = transcode(GzipDecompressor, Mmap.mmap(filename))
    else
        file = filename
    end
    return read_csv(file, delim, header)
end

read_features(filename::String, header::Vector{Symbol}=[:ensembleid, :genesymbol, :type]) =
    read_csv(filename, header=header)
read_barcodes(filename::String, header::Vector{Symbol}=[:barcode]) =
    read_csv(filename, header=header)

## Human Cell Atlas (HCA)
read_genes(filename::String) = read_csv(filename)
read_cells(filename::String) = read_csv(filename)
