var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = OmicsProfiles","category":"page"},{"location":"#OmicsProfiles","page":"Home","title":"OmicsProfiles","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for OmicsProfiles.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [OmicsProfiles]","category":"page"},{"location":"#OmicsProfiles.AnnotatedProfile","page":"Home","title":"OmicsProfiles.AnnotatedProfile","text":"AnnotatedProfile(omics, obs, obsindex)\nAnnotatedProfile(op, omicsname, obs, obsindex)\n\nArguments\n\nomics::Dict{Symbol,OmicsProfile}: A collection of omics profile with their names in keys.\nop::OmicsProfile: Single omics profile.\nomicsname::Symbol: The name of given OmicsProfile.\nobs::DataFrame: The dataframe contains meta information for observations.\nobsindex::Symbol: The index of dataframe obs.\n\n\n\n\n\n","category":"type"},{"location":"#OmicsProfiles.OmicsProfile","page":"Home","title":"OmicsProfiles.OmicsProfile","text":"OmicsProfile(countmat, var, varindex; T=float(eltype(countmat)))\n\nArguments\n\ncountmat: The count matrix for given omics and it will be set to key :count.\nvar::DataFrame: The dataframe contains meta information for features or variables.\nvarindex::Symbol: The index of dataframe var.\nT: Element type of countmat.\n\nExamples\n\njulia> using OmicsProfiles, DataFrames\n\njulia> r, c = (100, 500)\n(100, 500)\n\njulia> data = rand(0:100, r, c);\n\njulia> var = DataFrame(index=1:r, C=rand(r), D=rand(r));\n\njulia> prof = OmicsProfile(data, var, :index)\nOmicsProfile (nvar = 100):\n    var: index, C, D\n    layers: count\n\n\n\n\n\n","category":"type"},{"location":"#OmicsProfiles.Profile-Tuple{AbstractMatrix, Any, DataFrames.DataFrame, DataFrames.DataFrame}","page":"Home","title":"OmicsProfiles.Profile","text":"Profile(countmat, omicsname, var, obs; varindex=:genesymbol, obsindex::Symbol=:barcode,\n        T=float(eltype(countmat)))\n\nConstructor for establishing AnnotatedProfile with a OmicsProfile inside.\n\nArguments\n\ncountmat: The count matrix for given omics and it will be set to key :count.\nomicsname::Symbol: The name of given OmicsProfile.\nvar::DataFrame: The dataframe contains meta information for features or variables.\nobs::DataFrame: The dataframe contains meta information for observations.\nvarindex::Symbol: The index of dataframe var.\nobsindex::Symbol: The index of dataframe obs.\nT: Element type of countmat.\n\n\n\n\n\n","category":"method"},{"location":"#OmicsProfiles.mmread","page":"Home","title":"OmicsProfiles.mmread","text":"mmread(filename, retcoord=false)\n\nRead the contents of the Matrix Market file filename into a matrix, which will be either sparse or dense, depending on the Matrix Market format indicated by coordinate (coordinate sparse storage), or array (dense array storage).\n\nArguments\n\nfilename::String: The file to read.\nretcoord::Bool: If it is true, the rows, column and value vectors are returned along   with the header information.\n\n\n\n\n\n","category":"function"},{"location":"#OmicsProfiles.mmwrite-Tuple{String, SparseArrays.SparseMatrixCSC}","page":"Home","title":"OmicsProfiles.mmwrite","text":"mmwrite(filename, matrix)\n\nWrite a sparse matrix to .mtx file format.\n\nArguments\n\nfilename::String: The file to write.\nmatrix::SparseMatrixCSC: The sparse matrix to write.\n\n\n\n\n\n","category":"method"},{"location":"#OmicsProfiles.readinfo-Tuple{String}","page":"Home","title":"OmicsProfiles.readinfo","text":"readinfo(file)\n\nRead header information on the size and structure from file. The actual data matrix is not parsed.\n\nArguments\n\nfile: The filename or io stream.\n\n\n\n\n\n","category":"method"}]
}
