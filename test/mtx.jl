@testset "read/write mtx" begin
    mtx_filename = joinpath(OmicsProfiles.PROJECT_PATH, "test", "test.mtx")
    res = sparse(
        [5, 4, 1, 2, 6],
        [1, 5, 1, 4, 7],
        [1, 1, 1, 1, 1],
        11, 12
    )

    rows, cols, entries, rep, field, symm = OmicsProfiles.readinfo(mtx_filename)
    @test rows == 11
    @test cols == 12
    @test entries == 5
    @test rep == "coordinate"
    @test field == "integer"
    @test symm == "general"

    A = OmicsProfiles.mmread(mtx_filename)
    @test A isa SparseMatrixCSC
    @test A == res

    newfilename = replace(mtx_filename, "test.mtx" => "test_write.mtx")
    OmicsProfiles.mmwrite(newfilename, res)

    f = open(mtx_filename)
    sha_test = bytes2hex(sha256(read(f)))
    close(f)

    f = open(newfilename)
    sha_new = bytes2hex(sha256(read(f)))
    close(f)

    @test sha_test == sha_new
    rm(newfilename)
end