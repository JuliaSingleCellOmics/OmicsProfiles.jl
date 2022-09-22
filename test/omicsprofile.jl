@testset "omicsprofile" begin
    ngenes, ncells = (100, 500)
    data = rand(0:10, ngenes, ncells)
    var = DataFrame(index=1:ngenes, A=rand(ngenes), B=repeat([1], ngenes))
    @test_throws AssertionError OmicsProfile(data[1:50, :], var, :index)

    prof = OmicsProfile(data, var, :index)
    @test varnames(prof) == ["index", "A", "B"]
    @test countmatrix(prof) == data
    @test getvarindex(prof) == :index
    @test nrow(prof) == ngenes
    @test ncol(prof) == ncells
    @test nvar(prof) == ngenes
    @test maximum(prof) == 10
    @test minimum(prof) == 0
    @test size(prof) == (ngenes, ncells)
    @test axes(prof) == (Base.OneTo(ngenes), Base.OneTo(ncells))

    prof2 = copy(prof)
    @test prof2 !== prof
    @test prof == prof2

    prof.var[!, :newindex] = collect(ngenes:-1:1)
    setvarindex!(prof, :newindex)
    @test getvarindex(prof) == :newindex
    @test_throws ArgumentError setvarindex!(prof, :D)
    setvarindex!(prof, :index)
    select!(prof.var, Not(:newindex))

    setlayer!(prof, rand(ngenes, ncells), :a)
    setlayer!(prof, rand(ngenes, ngenes), :b)
    @test collect(layernames(prof)) == [:a, :b, :count]
    @test size(getlayer(prof, :a)) == (ngenes, ncells)
    @test size(getlayer(prof, :b)) == (ngenes, ngenes)

    setpipeline!(prof, Dict(:a => 1), :qc_metrics)
    setpipeline!(prof, Dict(:b => 2), :normalize)
    setpipeline!(prof, Dict(:c => 3), :log_transform)
    @test collect(keys(prof.pipeline)) == [:qc_metrics, :normalize, :log_transform]
    @test getpipeline(prof, :qc_metrics)[:a] == 1
    @test getpipeline(prof, :normalize)[:b] == 2
    @test getpipeline(prof, :log_transform)[:c] == 3

    @test repr(prof) == "OmicsProfile (nvar = 100):\n    var: index, A, B\n    layers: a, b, count\n    pipeline: qc_metrics => normalize => log_transform"

    @test vec(geneexpr(prof, 50)) == data[50, :]
    @test vec(geneexpr(prof, 50, :a)) == prof.layers[:a][50, :]
end
