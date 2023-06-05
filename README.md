[![Build Status](https://github.com/ClapeyronThermo/GCIdentifier.jl/workflows/CI/badge.svg)](https://github.com/ClapeyronThermo/GCIdentifier.jl/actions)
[![codecov](https://codecov.io/gh/ClapeyronThermo/GCIdentifier.jl/branch/master/graph/badge.svg?token=ZVGGR4AAFB)](https://codecov.io/gh/ClapeyronThermo/GCIdentifier.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://clapeyronthermo.github.io/GCIdentifier.jl/dev)
[![project chat](https://img.shields.io/badge/zulip-join_chat-brightgreen.svg)](https://julialang.zulipchat.com/#narrow/stream/265161-Clapeyron.2Ejl)

# GCIdentifier.jl

Identifies subgroups, given a SMILES.

## Introduction

Most Group contribution-based (GC) property prediction models require, at least, a list of substructure fragments ("groups") that conform a molecule, and their amount. While it is possible to build group lists for simple molecules rather easily (e.g water, carbon dioxide, ethanol), it becames harder with increasing molecule size. Furthermore different GC models (UNIFAC Activity model, SAFT-γ-mie equation of State, Joback ideal contributions), have different definitions on what groups are supported and how are named.

Luckily, the process can be automated. By using SMILES as imput, 

## Usage

this library exports:
- `get_groups_from_smiles`: the main function, given a SMILES string, it will return the list of groups and their respective amounts. 
- `RDKitLib`, `MolecularGraphJL` : structs to select the underlying molecular library used to perform substructure search.
- Group lists. the current groups are as defined by the corresponding equation of state in [`Clapeyron.jl`](https://github.com/ClapeyronThermo/Clapeyron.jl):
    - `UNIFACGroups`
    - `SAFTGammaMieGroups`
    - `JobackGroups`
    - `gcPCSAFTGroups`

you can query the functions in the following way:

```
julia> get_groups_from_smiles("CCO",UNIFACGroups)
("CCO", ["CH3" => 1, "CH2" => 1, "OH(P)" => 1])
```