@testset "omicsprofile" begin
    ngenes, ncells = (100, 500)
    data = rand(0:10, ngenes, ncells)
    var = DataFrame(index=string.(1:ngenes), A=rand(ngenes), B=repeat([1], ngenes))
    @test_throws AssertionError OmicsProfile(data[1:50, :], var, :index)

    prof = OmicsProfile(data, var, :index)
    @test varnames(prof) == ["index", "A", "B"]
    @test prof.count == data
    @test prof.varindex == :index
    @test nrow(prof) == ngenes
    @test ncol(prof) == ncells
    @test nvar(prof) == ngenes
    @test maximum(prof) == 10
    @test minimum(prof) == 0
    @test size(prof) == (ngenes, ncells)
    @test Base.axes(prof) == (Base.OneTo(ngenes), Base.OneTo(ncells))

    prof2 = copy(prof)
    @test prof2 !== prof
    @test prof == prof2

    prof.var[!, :newindex] = collect(ngenes:-1:1)
    prof.varindex = :newindex
    @test prof.varindex == :newindex
    @test_throws ArgumentError prof.varindex = :D
    prof.varindex = :index
    select!(prof.var, Not(:newindex))

    prof.varm[:pcs] = rand(ngenes, 50)
    @test size(prof.pcs) == (ngenes, 50)
    prof.pcs = rand(ngenes, 80)
    @test size(prof.pcs) == (ngenes, 80)

    prof.layers[:a] = rand(ngenes, ncells)
    prof.layers[:b] = rand(ngenes, ngenes)
    @test collect(layernames(prof)) == [:a, :b, :count]
    @test size(prof.a) == (ngenes, ncells)
    @test size(prof.b) == (ngenes, ngenes)
    prof.a = collect(prof.count)
    @test prof.a == prof.count

    prof.pipeline[:qc_metrics] = Dict(:a => 1)
    prof.pipeline[:normalize] = Dict(:b => 2)
    prof.pipeline[:log_transform] = Dict(:c => 3)
    @test collect(keys(prof.pipeline)) == [:qc_metrics, :normalize, :log_transform]
    @test prof.qc_metrics[:a] == 1
    @test prof.normalize[:b] == 2
    @test prof.log_transform[:c] == 3
    prof.qc_metrics = Dict(:a => 5)
    @test prof.qc_metrics[:a] == 5

    @test repr(prof) == "OmicsProfile (nvar = 100):\n    var: index, A, B\n    varm: pcs\n    layers: a, b, count\n    pipeline: qc_metrics => normalize => log_transform"

    @test vec(prof.count["50", :]) == data[50, :]
    @test vec(prof.a["50", :]) == prof.layers[:a][50, :]
    @test prof.a[["50", "100"], :] == prof.layers[:a][[50, 100], :]

    # indexing
    prof2 = prof[1:50, :]
    @test prof2.count == prof.count[1:50, :]

    idx = rand([false,true], ngenes)
    prof3 = prof[idx, :]
    @test prof3.count == prof.count[idx, :]
    @test prof3.var == prof.var[idx, :]
end
