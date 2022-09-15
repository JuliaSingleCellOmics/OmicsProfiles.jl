@testset "omicsprofile" begin
    r, c = (100, 500)
    data = rand(0:10, r, c)
    var = DataFrame(index=1:r, A=rand(r), B=repeat([1], r))
    @test_throws AssertionError OmicsProfile(data[1:50, :], var, :index)

    prof = OmicsProfile(data, var, :index)
    @test varnames(prof) == ["index", "A", "B"]
    @test countmatrix(prof) == data
    @test getvarindex(prof) == :index
    @test nrow(prof) == r
    @test ncol(prof) == c
    @test nvar(prof) == r
    @test maximum(prof) == 10
    @test minimum(prof) == 0
    @test size(prof) == (r, c)
    @test axes(prof) == (Base.OneTo(r), Base.OneTo(c))

    prof2 = copy(prof)
    @test prof2 !== prof
    @test countmatrix(prof) == countmatrix(prof2)
    @test prof.var == prof2.var

    prof.var[!, :newindex] = collect(r:-1:1)
    setvarindex!(prof, :newindex)
    @test getvarindex(prof) == :newindex
    @test_throws ArgumentError setvarindex!(prof, :D)
    setvarindex!(prof, :index)
    select!(prof.var, Not(:newindex))

    setlayer!(prof, rand(r, c), :a)
    setlayer!(prof, rand(r, r), :b)
    @test collect(layernames(prof)) == [:a, :b, :count]
    @test size(getlayer(prof, :a)) == (r, c)
    @test size(getlayer(prof, :b)) == (r, r)

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
