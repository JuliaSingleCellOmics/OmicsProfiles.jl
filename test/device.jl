@testset "device" begin
    ngenes, ncells = (100, 500)
    data = rand(0:10, ngenes, ncells)
    var = DataFrame(genesymbol=1:ngenes, A=rand(ngenes))
    obs = DataFrame(barcode=1:ncells, C=rand(ncells))

    omic = OmicsProfile(data, var, :genesymbol)
    omic.varm[:pcs] = rand(ngenes, 50)

    ap = AnnotatedProfile(omic, :RNA, obs, :barcode)
    ap.obsm[:pca] = rand(50, ncells)

    ap = ap |> togpu

    # OmicsProfile
    @test ap.RNA.varm[:pcs] isa CuMatrix
    @test ap.RNA.count isa CuMatrix

    # AnnotatedProfile
    @test ap.obsm[:pca] isa CuMatrix
end
