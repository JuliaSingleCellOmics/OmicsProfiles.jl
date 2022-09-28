@testset "filter" begin
    ngenes, ncells = (100, 500)
    data = rand(0:10, ngenes, ncells)
    var = DataFrame(genesymbol=1:ngenes, A=rand(ngenes), B=repeat([1], ngenes))
    prof = OmicsProfile(data, var, :genesymbol)
    prof.varm[:a] = rand(ngenes, 50)

    @testset "OmicsProfile" begin
        prof2 = filter(:A => x -> x > 0, prof)
        @test prof2 !== prof
        @test prof2.var == prof.var[prof.var.A .> 0, :]
        @test prof2.count == prof.count[prof.var.A .> 0, :]
        @test prof2.varm[:a] == prof.varm[:a][prof.var.A .> 0, :]

        filter!(:A => x -> x > 0, prof2)
        @test prof2.var == prof.var[prof.var.A .> 0, :]
        @test prof2.count == prof.count[prof.var.A .> 0, :]
        @test prof2.varm[:a] == prof.varm[:a][prof.var.A .> 0, :]
    end

    @testset "AnnotatedProfile" begin
        obs = DataFrame(barcode=1:ncells, C=rand(ncells), D=rand(ncells))
        ap = AnnotatedProfile(prof, :RNA, obs, :barcode)
        ap.obsm[:a] = rand(10, ncells)

        ap2 = filter(:A => x -> x > 0, ap)
        @test ap2 !== ap
        @test ap2.RNA.var == ap.RNA.var[prof.var.A .> 0, :]
        @test ap2.RNA.count == ap.RNA.count[prof.var.A .> 0, :]
        @test ap2.RNA.varm[:a] == ap.RNA.varm[:a][prof.var.A .> 0, :]

        filter!(:A => x -> x > 0, ap2)
        @test ap2.RNA.var == ap.RNA.var[prof.var.A .> 0, :]
        @test ap2.RNA.count == ap.RNA.count[prof.var.A .> 0, :]
        @test ap2.RNA.varm[:a] == ap.RNA.varm[:a][prof.var.A .> 0, :]

        ap2 = filter(:C => x -> x > 0, ap)
        @test ap2 !== ap
        @test ap2.obs == ap.obs[ap.obs.C .> 0, :]
        @test ap2.RNA.count == ap.RNA.count[:, ap.obs.C .> 0]
        @test ap2.obsm[:a] == ap.obsm[:a][:, ap.obs.C .> 0]

        filter!(:C => x -> x > 0, ap2)
        @test ap2.obs == ap.obs[ap.obs.C .> 0, :]
        @test ap2.RNA.count == ap.RNA.count[:, ap.obs.C .> 0]
        @test ap2.obsm[:a] == ap.obsm[:a][:, ap.obs.C .> 0]
    end
end
