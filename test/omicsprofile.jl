@testset "omicsprofile" begin
    r, c = (100, 500)
    data = rand(0:10, r, c)
    obs = DataFrame(index=1:c, A=rand(c), B=rand(c))
    var = DataFrame(index=1:r, C=rand(r), D=rand(r))
    @test_throws AssertionError OmicsProfile(data, obs, var)

    prof = OmicsProfile(data, var, obs)
    prof.layers[:a] = rand(r, c)
    prof.layers[:b] = rand(c, r)
    @test obsnames(prof) == ["index", "A", "B"]
    @test varnames(prof) == ["index", "C", "D"]
    @test nrow(prof) == r
    @test ncol(prof) == c
    @test nvar(prof) == r
    @test nobs(prof) == c
    @test maximum(prof) == 10
    @test minimum(prof) == 0
    @test size(prof) == (r, c)
    @test axes(prof) == (Base.OneTo(r), Base.OneTo(c))
    @test_throws AssertionError prof.obs = var
    @test_throws AssertionError prof.var = obs
    @test size(prof.layers[:a]) == (r, c)
    @test size(prof.layers[:b]) == (c, r)
    @test vec(get_gene_expr(prof, 50)) == data[50, :]
    @test vec(get_gene_expr(prof, 50, :a)) == prof.layers[:a][50, :]
    @test isempty(get_gene_expr(prof, r+11, :a))

    prof2 = copy(prof)
    @test prof2 !== prof
    @test prof.data == prof2.data
    @test prof.var == prof2.var
    @test prof.obs == prof2.obs
    
    prof2 = filter(:C => x -> x > 0, prof)
    @test prof2 !== prof
    @test prof2.var == prof.var[prof.var.C .> 0, :]
    @test prof2.data == prof.data[prof.var.C .> 0, :]
    @test prof2.obs == prof.obs
    @test prof2.layers[:a] == prof.layers[:a][prof.var.C .> 0, :]
    @test prof2.layers[:b] == prof.layers[:b]

    filter!(:C => x -> x > 0, prof2)
    @test prof2.var == prof.var[prof.var.C .> 0, :]
    @test prof2.data == prof.data[prof.var.C .> 0, :]
    @test prof2.obs == prof.obs
    @test prof2.layers[:a] == prof.layers[:a][prof.var.C .> 0, :]
    @test prof2.layers[:b] == prof.layers[:b]

    prof2 = prof[1:50, :]
    @test prof2.data == prof.data[1:50, :]

    idx = rand([false,true], r)
    prof3 = prof[idx, :]
    @test prof3.data == prof.data[idx, :]
    @test prof3.var == prof.var[idx, :]
    @test prof3.obs == prof.obs

    idx2 = rand([false,true], c)
    prof4 = prof[:, idx2]
    @test prof4.data == prof.data[:, idx2]
    @test prof4.var == prof.var
    @test prof4.obs == prof.obs[idx2, :]
end
