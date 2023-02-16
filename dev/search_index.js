var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = OmicsProfiles","category":"page"},{"location":"#OmicsProfiles","page":"Home","title":"OmicsProfiles","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for OmicsProfiles.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [OmicsProfiles]","category":"page"},{"location":"#OmicsProfiles.AnnotatedProfile","page":"Home","title":"OmicsProfiles.AnnotatedProfile","text":"AnnotatedProfile(omics, obs, obsindex)\nAnnotatedProfile(op, omicsname, obs, obsindex)\n\nArguments\n\nomics::Dict{Symbol,OmicsProfile}: A collection of omics profile with their names in keys.\nop::OmicsProfile: Single omics profile.\nomicsname::Symbol: The name of given OmicsProfile.\nobs::DataFrame: The dataframe contains meta information for observations.\nobsindex::Symbol: The index of dataframe obs.\n\n\n\n\n\n","category":"type"},{"location":"#OmicsProfiles.OmicsProfile","page":"Home","title":"OmicsProfiles.OmicsProfile","text":"OmicsProfile(countmat, var, varindex; T=float(eltype(countmat)))\n\nArguments\n\ncountmat: The count matrix for given omics and it will be set to key :count.\nvar::DataFrame: The dataframe contains meta information for features or variables.\nvarindex::Symbol: The index of dataframe var.\nT: Element type of countmat.\n\nExamples\n\njulia> using OmicsProfiles, DataFrames\n\njulia> r, c = (100, 500)\n(100, 500)\n\njulia> data = rand(0:100, r, c);\n\njulia> var = DataFrame(index=1:r, C=rand(r), D=rand(r));\n\njulia> prof = OmicsProfile(data, var, :index)\nOmicsProfile (nvar = 100):\n    var: index, C, D\n    layers: count\n\n\n\n\n\n","category":"type"},{"location":"#Base.filter-Tuple{Pair{Symbol}, AbstractProfile}","page":"Home","title":"Base.filter","text":"filter(col => func, p)\n\nReturns a copy of profile p containing data from p for which func returns true.\n\nArguments\n\ncol: Column name to filter.\nfunc: A predicate function which returns true/false.\np::AbstractProfile: A profile to filter on.\n\nExamples\n\njulia> using OmicsProfiles, DataFrames\n\njulia> ngenes, ncells = (100, 500)\n(100, 500)\n\njulia> X = rand(0:10, ngenes, ncells);\n\njulia> var = DataFrame(genesymbol=1:ngenes, A=rand(ngenes));\n\njulia> prof = OmicsProfile(X, var, :genesymbol)\nOmicsProfile (nvar = 100):\n    var: genesymbol, A\n    layers: count\n\njulia> prof.varm[:a] = rand(ngenes, 50);\n\njulia> prof2 = filter(:A => x -> x > 0.5, prof)\nOmicsProfile (nvar = 46):\n    var: genesymbol, A\n    varm: a\n    layers: count\n\njulia> obs = DataFrame(barcode=1:ncells, B=rand(ncells));\n\njulia> ap = AnnotatedProfile(prof, :RNA, obs, :barcode)\nAnnotatedProfile (nobs = 500):\n    obs: barcode, B\nRNA => OmicsProfile (nvar = 100):\n    var: genesymbol, A\n    varm: a\n    layers: count\n\njulia> ap2 = filter(:B => x -> x > 0.5, ap)\nAnnotatedProfile (nobs = 253):\n    obs: barcode, B\nRNA => OmicsProfile (nvar = 100):\n    var: genesymbol, A\n    varm: a\n    layers: count\n\n\n\n\n\n","category":"method"},{"location":"#OmicsProfiles.Profile-Tuple{AbstractMatrix, Any, DataFrames.DataFrame, DataFrames.DataFrame}","page":"Home","title":"OmicsProfiles.Profile","text":"Profile(countmat, omicsname, var, obs; varindex=:genesymbol, obsindex::Symbol=:barcode,\n        T=float(eltype(countmat)))\n\nConstructor for establishing AnnotatedProfile with a OmicsProfile inside.\n\nArguments\n\ncountmat: The count matrix for given omics and it will be set to key :count.\nomicsname::Symbol: The name of given OmicsProfile.\nvar::DataFrame: The dataframe contains meta information for features or variables.\nobs::DataFrame: The dataframe contains meta information for observations.\nvarindex::Symbol: The index of dataframe var.\nobsindex::Symbol: The index of dataframe obs.\nT: Element type of countmat.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"<figure>\n    <img src=\"../assets/profile.png\" width=\"80%\" alt=\"profile.png\" /><br>\n    <figcaption><em>AnnotatedProfile supports multi-omics data with respect to the same set of cells, while OmicsProfile holds single omic data with layers which contains data with the same set of features.</em></figcaption>\n</figure>","category":"page"}]
}