@testset "annotatedprofile" begin
    r, c = (100, 500)
    data = rand(0:10, r, c)
    var = DataFrame(genesymbol=1:r, A=rand(r), B=repeat([1], r))
    obs = DataFrame(barcode=1:c, C=rand(c), D=rand(c))

    omic = OmicsProfile(data, var, :genesymbol)
    ap = AnnotatedProfile(omic, :RNA, obs, :barcode)
    @test collect(omicsnames(ap)) == [:RNA]
    @test varnames(ap, :RNA) == ["genesymbol", "A", "B"]
    @test obsnames(ap) == ["barcode", "C", "D"]
    @test countmatrix(ap, :RNA) == data
    @test getobsindex(ap) == :barcode

    @test nvar(ap, :RNA) == r
    @test nobs(ap) == c

    setpipeline!(ap, Dict(:a => 1), :merge_omics)
    setpipeline!(ap, Dict(:b => 2), :integration)
    @test collect(keys(ap.pipeline)) == [:merge_omics, :integration]
    @test getpipeline(ap, :merge_omics)[:a] == 1
    @test getpipeline(ap, :integration)[:b] == 2

    @test repr(ap) == "AnnotatedProfile (nobs = 500):\n    obs: barcode, C, D\n    pipeline: merge_omics => integration" *
        "\nRNA => OmicsProfile (nvar = 100):\n    var: genesymbol, A, B\n    layers: count"
    p = Profile(data, :RNA, var, obs)
end
