abstract type AbstractDevice end

struct CPUDevice <: AbstractDevice end

passingfunc(::CPUDevice) = cpu_sparsify

cpu_sparsify(x) = adapt(Array, x) |> sparse

struct CUDADevice <: AbstractDevice end

passingfunc(::CUDADevice) = densify_gpu

densify_gpu(x) = x |> Matrix |> cu

function todevice(p::OmicsProfile, device::AbstractDevice)
    f = passingfunc(device)
    for k in keys(p.varm)
        p.varm[k] = p.varm[k] |> f
    end

    for k in keys(p.layers)
        p.layers[k] = p.layers[k] |> f
    end

    return p
end

function todevice(p::AnnotatedProfile, device::AbstractDevice)
    f = passingfunc(device)
    for k in keys(p.omics)
        p.omics[k] = todevice(p.omics[k], device)
    end

    for k in keys(p.obsm)
        p.obsm[k] = p.obsm[k] |> f
    end

    return p
end

"""
    tocpu(x)

Returns `x` which has array contents moved to main memory.

# Arguments

- `x`: The object to move to CPU memory. Currently supports `OmicsProfile` and `AnnotatedProfile`.

# Examples

```jldoctest
julia> using CUDA, DataFrames

julia> ngenes, ncells = (100, 500)
(100, 500)

julia> data = rand(0:10, ngenes, ncells);

julia> var = DataFrame(genesymbol=1:ngenes);

julia> prof = Profile(data, :RNA, var, obs) |> togpu
AnnotatedProfile (nobs = 500):
    obs: barcode
RNA => OmicsProfile (nvar = 100):
    var: genesymbol
    layers: count

julia> typeof(prof.RNA.layers[:count])
CuArray{Float32, 2, CUDA.Mem.DeviceBuffer}

julia> prof = prof |> tocpu
AnnotatedProfile (nobs = 500):
    obs: barcode
RNA => OmicsProfile (nvar = 100):
    var: genesymbol
    layers: count

julia> typeof(prof.RNA.layers[:count])
SparseArrays.SparseMatrixCSC{Float32, Int64}
```

See also [`togpu`](@ref) for moving `x` to GPU.
"""
tocpu(x) = todevice(x, CPUDevice())

"""
    togpu(x)

Returns `x` which has array contents moved to GPU memory.

# Arguments

- `x`: The object to move to GPU memory. Currently supports `OmicsProfile` and `AnnotatedProfile`.

# Examples

```jldoctest
julia> using CUDA, DataFrames

julia> ngenes, ncells = (100, 500)
(100, 500)

julia> data = rand(0:10, ngenes, ncells);

julia> var = DataFrame(genesymbol=1:ngenes);

julia> prof = Profile(data, :RNA, var, obs) |> togpu
AnnotatedProfile (nobs = 500):
    obs: barcode
RNA => OmicsProfile (nvar = 100):
    var: genesymbol
    layers: count

julia> typeof(prof.RNA.layers[:count])
CuArray{Float32, 2, CUDA.Mem.DeviceBuffer}
```

See also [`togpu`](@ref) for moving `x` to GPU.
"""
togpu(x) = todevice(x, CUDADevice())
