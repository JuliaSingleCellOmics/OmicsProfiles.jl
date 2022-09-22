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
    @test countmatrix(ap, :RNA) == data
    @test getobsindex(ap) == :barcode

    @test nvar(ap, :RNA) == ngenes
    @test nobs(ap) == ncells

    setpipeline!(ap, Dict(:a => 1), :merge_omics)
    setpipeline!(ap, Dict(:b => 2), :integration)
    @test collect(keys(ap.pipeline)) == [:merge_omics, :integration]
    @test getpipeline(ap, :merge_omics)[:a] == 1
    @test getpipeline(ap, :integration)[:b] == 2

    @test repr(ap) == "AnnotatedProfile (nobs = 500):\n    obs: barcode, C, D\n    pipeline: merge_omics => integration" *
        "\nRNA => OmicsProfile (nvar = 100):\n    var: genesymbol, A, B\n    layers: count"

    ap2 = copy(ap)
    @test ap2 !== ap
    @test ap2 == ap

    p = Profile(data, :RNA, var, obs)
    @test p isa AnnotatedProfile
end
