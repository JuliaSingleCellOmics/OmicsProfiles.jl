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

tocpu(x) = todevice(x, CPUDevice())
togpu(x) = todevice(x, CUDADevice())
