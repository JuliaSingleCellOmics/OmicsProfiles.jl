# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0]

- Support `OmicsProfile` and `AnnotatedProfile` moving to GPU by `tocpu` and `togpu`
- Support reading distinct name of feature files, e.g. `genes.csv` and `features.csv`

## [0.1.7]

- improve documentation

## [0.1.6]

- `OmicsProfile` and `AnnotatedProfile` support `propertynames`
- support `read_10x` to read 10x folder

## [0.1.5]

- add `filter` for `OmicsProfile` and `AnnotatedProfile`

## [0.1.4]

- support `p.RNA.count[:gene1, :]` syntax to get gene expression

## [0.1.3]

- support `read_mtx` and `read_csv`

## [0.1.2]

- override `getproperty` and `setproperty!`

## [0.1.1]

- `OmicsProfile` and `AnnotatedProfile` supports `==` and `copy`
