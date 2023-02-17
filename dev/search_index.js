var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = OmicsProfiles","category":"page"},{"location":"#[OmicsProfiles](https://github.com/yuehhua/OmicsProfiles.jl)","page":"Home","title":"OmicsProfiles","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<figure>\n    <img src=\"assets/profile.png\" width=\"80%\" alt=\"profile.png\" /><br>\n    <figcaption><em>AnnotatedProfile supports multi-omics data with respect to the same set of cells, while OmicsProfile holds single omic data with layers which contains data with the same set of features.</em></figcaption>\n</figure>","category":"page"},{"location":"","page":"Home","title":"Home","text":"Profile\nOmicsProfile\nAnnotatedProfile\nread_10x","category":"page"},{"location":"#OmicsProfiles.Profile","page":"Home","title":"OmicsProfiles.Profile","text":"Profile(countmat, omicsname, var, obs; varindex=:genesymbol, obsindex::Symbol=:barcode,\n        T=float(eltype(countmat)))\n\nConstructor for establishing AnnotatedProfile with a OmicsProfile inside.\n\nArguments\n\ncountmat: The count matrix for given omics and it will be set to key :count.\nomicsname::Symbol: The name of given OmicsProfile.\nvar::DataFrame: The dataframe contains meta information for features or variables.\nobs::DataFrame: The dataframe contains meta information for observations.\nvarindex::Symbol: The index of dataframe var.\nobsindex::Symbol: The index of dataframe obs.\nT: Element type of countmat.\n\n\n\n\n\n","category":"function"},{"location":"#OmicsProfiles.OmicsProfile","page":"Home","title":"OmicsProfiles.OmicsProfile","text":"OmicsProfile(countmat, var, varindex; T=float(eltype(countmat)))\n\nArguments\n\ncountmat: The count matrix for given omics and it will be set to key :count.\nvar::DataFrame: The dataframe contains meta information for features or variables.\nvarindex::Symbol: The index of dataframe var.\nT: Element type of countmat.\n\nExamples\n\njulia> r, c = (100, 500)\n(100, 500)\n\njulia> data = rand(0:100, r, c);\n\njulia> var = DataFrame(index=1:r, C=rand(r), D=rand(r));\n\njulia> prof = OmicsProfile(data, var, :index)\nOmicsProfile (nvar = 100):\n    var: index, C, D\n    layers: count\n\n\n\n\n\n","category":"type"},{"location":"#OmicsProfiles.AnnotatedProfile","page":"Home","title":"OmicsProfiles.AnnotatedProfile","text":"AnnotatedProfile(omics, obs, obsindex)\nAnnotatedProfile(op, omicsname, obs, obsindex)\n\nArguments\n\nomics::Dict{Symbol,OmicsProfile}: A collection of omics profile with their names in keys.\nop::OmicsProfile: Single omics profile.\nomicsname::Symbol: The name of given OmicsProfile.\nobs::DataFrame: The dataframe contains meta information for observations.\nobsindex::Symbol: The index of dataframe obs.\n\n\n\n\n\n","category":"type"},{"location":"#OmicsProfiles.read_10x","page":"Home","title":"OmicsProfiles.read_10x","text":"read_10x(path; make_unique=true, omicsname=:RNA, varnames=[:ensembleid, :genesymbol],\n    obsnames=[:barcode], varindex=:genesymbol, obsindex=:barcode)\n\nReads 10x dataset from path and returns an AnnotatedProfile.\n\nArguments\n\npath::AbstractString: The path to 10x dataset folder.\nmake_unique::Bool: To make varindex unique or not if given varindex is not unique.\nomicsname::Symbol: The name of given sequencing data.\nvarnames::AbstractVector: The header or column names of feature file.\nobsnames::AbstractVector: The header or column names of barcode file.\nvarindex::Symbol: The index of dataframe var.\nobsindex::Symbol: The index of dataframe obs.\n\nExamples\n\njulia> using SnowyOwl\n\njulia> prof = read_10x(SnowyOwl.Dataset.pbmc3k_folder(), varindex=:ensembleid)\nAnnotatedProfile (nobs = 2700):\n    obs: barcode\nRNA => OmicsProfile (nvar = 32738):\n    var: ensembleid, genesymbol\n    layers: count\n\n\n\n\n\n","category":"function"}]
}