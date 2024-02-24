# Group Assignment
The primary function of GCIdentifier is to automatically assign groups to a molecule given a specific group-contribution method. This assumes that the groups in the method being used are extensive enough to cover most molecules that might be of interest.

Let's consider the case where we want to get the groups for ibuprofen from the UNIFAC group-contribution method. The SMILES representation for ibuprofen is CC(Cc1ccc(cc1)C(C(=O)O)C)C. To get the corresponding groups, simply use:
```julia
julia> (smiles,groups) = get_groups_from_smiles("CC(Cc1ccc(cc1)C(C(=O)O)C)C", UNIFACGroups)
("CC(Cc1ccc(cc1)C(C(=O)O)C)C", ["COOH" => 1, "CH3" => 3, "CH" => 1, "ACH" => 4, "ACCH2" => 1, "ACCH" => 1])
```
where `smiles` will output the molecular SMILES and `groups` is a vector of pairs where the first element is the group name and the second is the number of times the group occurs within the molecule. If this function fails, it is usually because, either, the SMILES is unphysical, or the method used doesn't cover all atoms present within a molecule. In the case of the latter, one can still obtain the groups that have been identified by specifying `check=false` in the optional arguments. For example, SAFT-$\gamma$ Mie does not have the functional group for ketones:
```julia
julia> (smiles,groups) = get_groups_from_smiles("CCC(=O)CC", SAFTgammaMieGroups; check=false)
("CCC(=O)CC", ["CH3" => 2, "CH2" => 2])
```
To propose new groups that cover the missing atoms, take a look at our [`find_missing_groups` function](./missing_groups.md).

## Connectivity
There are certain group-contribution approaches, such as gcPCP-SAFT and s-SAFT-$\gamma$ Mie where information about how groups are linked to each other is required. It is possible to obtain information about the connectivity between groups from GCIdentifier by simply specifying `connectivity=true` within the `get_groups_from_smiles` function. For example, in the case of acetone:
```julia
julia> (smiles,groups,connectivity) = get_groups_from_smiles("CC(=O)C", gcPPCSAFTGroups; connectivity=true)
("CC(=O)C", ["C=O" => 1, "CH3" => 2], [("C=O", "CH3") => 2])
```
where `connectivity` is a vector of pairs where the first element is the groups involved in the link and the second element is the number of times this link occurs.

## Extensions
Unfortunately, the process of obtaining the SMILES representation of a molecule can be itself a challenge. As such, we have included an extension of GCIdentifier where, if called in conjunction with [ChemicalIdentifiers](https://github.com/longemen3000/ChemicalIdentifiers.jl), one can simply specify the molecule name:
```julia
julia> using ChemicalIdentifiers

julia> (component,groups) = get_groups_from_name("ibuprofen",UNIFACGroups)
("ibuprofen", ["COOH" => 1, "CH3" => 3, "CH" => 1, "ACH" => 4, "ACCH2" => 1, "ACCH" => 1])
```
This should greatly simplify the use of GCIdentifier. These groups can then be used in packages such as [Clapeyron](https://github.com/ClapeyronThermo/Clapeyron.jl) to be used to obtain our desired properties, such as the solubility of ibuprofen in water:
```julia
julia> using Clapeyron

julia> liquid = UNIFAC(["water",(component,groups)])

julia> solid = SolidHfus(["water","ibuprofen"])

julia> model = CompositeModel(["water","ibuprofen"]; solid=solid, liquid=liquid)

julia> sle_solubility(model,1e5,298.15, [1.,0.]; solute=["ibuprofen"])[2]
5.547514144547524e-7
```