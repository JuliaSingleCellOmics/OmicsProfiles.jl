@testset "annotatedprofile" begin
    ngenes, ncells = (100, 500)
    data = rand(0:10, ngenes, ncells)
    var = DataFrame(genesymbol=1:ngenes, A=rand(ngenes), B=repeat([1], ngenes))
    obs = DataFrame(barcode=1:ncells, C=rand(ncells), D=rand(ncells))

    omic = OmicsProfile(data, var, :genesymbol)
    ap = AnnotatedProfile(omic, :RNA, obs, :barcode)
    @test collect(omicsnames(ap)) == [:RNA]
    @test varnames(ap, :RNA) == ["genesymbol", "A", "B"]
    @test obsnames(ap) == ["barcode", "C", "D"]
    @test ap.RNA.count == data
    @test ap.obsindex == :barcode

    @test nvar(ap, :RNA) == ngenes
    @test nobs(ap) == ncells

    ap.obs[!, :newindex] = collect(ncells:-1:1)
    ap.obsindex = :newindex
    @test ap.obsindex == :newindex
    @test_throws ArgumentError ap.obsindex = :A
    ap.obsindex = :barcode
    select!(ap.obs, Not(:newindex))

    ap.obsm[:pca] = rand(50, ncells)
    @test size(ap.pca) == (50, ncells)
    ap.pca = rand(30, ncells)
    @test size(ap.pca) == (30, ncells)

    ap.pipeline[:merge_omics] = Dict(:a => 1)
    ap.pipeline[:integration] = Dict(:b => 2)
    @test collect(keys(ap.pipeline)) == [:merge_omics, :integration]
    @test ap.merge_omics[:a] == 1
    @test ap.integration[:b] == 2
    ap.merge_omics = Dict(:a => 5)
    @test ap.merge_omics[:a] == 5
    @test propertynames(ap) == (:omics, :obs, :obsindex, :obsm, :obsgraphs, :pipeline,
        :integration, :pca, :merge_omics)

    @test repr(ap) == "AnnotatedProfile (nobs = 500):\n    obs: barcode, C, D\n    obsm: pca\n    pipeline: merge_omics => integration" *
        "\nRNA => OmicsProfile (nvar = 100):\n    var: genesymbol, A, B\n    layers: count"

    ap2 = copy(ap)
    @test ap2 !== ap
    @test ap2 == ap

    p = Profile(data, :RNA, var, obs)
    @test p isa AnnotatedProfile

    idx = rand([false,true], ncells)
    # ap2 = ap[:, idx]
    # @test ap2.RNA.count == ap.RNA.count[:, idx]
    # @test ap2.RNA.var == ap.RNA.var
    # @test ap2.obs == ap.obs[idx, :]
end
