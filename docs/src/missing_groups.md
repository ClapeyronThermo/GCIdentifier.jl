# Find Missing Groups
In some cases, our group contribution method will not have every group needed to cover every atom within a molecule. The example we gave previously was how SAFT-$\gamma$ Mie lacked a group for ketones:
```julia
julia> (smiles,groups) = get_groups_from_smiles("CCC(=O)CC", SAFTgammaMieGroups)
ERROR: Could not find all groups for CCC(=O)CC
```
In this case, it could be possible that GCIdentifier simply hasn't included the missing group within our database or perhaps that group needs to be parameterised. To find out which groups might help fill-in the missing atoms, one can use `find_missing_groups_from_smiles`:
```julia
julia> groups = find_missing_groups_from_smiles("CCC(=O)CC", SAFTgammaMieGroups)
3-element Vector{GCPair}:
 GCPair("[CX3;H0;!R]", "C=")
 GCPair("[OX1;H0;!R]", "O=")
 GCPair("[CX3;H0;!R](=[OX1;H0;!R])", "C=O=")
```
where `groups` is a vector of `GCPair`s with proposed names of the groups. As we can see, GCIdentifier has proposed three potential groups, where the last is a combination of the first two. From this list, the users can decide which group might be best to parameterise. However, we also have our own [internal heuristics](./api.md) for proposing a _minimal_ group representation of molecules within GCIdentifier which will reduce this list to what we recommend:
```julia
julia> groups = find_missing_groups_from_smiles("CCC(=O)CC", SAFTgammaMieGroups; reduced=true)
1-element Vector{GCPair}:
 GCPair("[CX3;H0;!R](=[OX1;H0;!R])", "C=O=")
```
In the case of the ketone, we would only really want to parameterise this group.

## Automatically fragment a molecule
Let us now consider an extreme case where we are trying to fragment a molecule into groups with no reference group-contribution approach. In the case of our ketone, the code will fragment the molecule into atomic groups, along with combinations of those atomic groups:
```julia
julia> groups = find_missing_groups_from_smiles("CCCC(=O)CCC")
5-element Vector{GCPair}:
 GCPair("[CX4;H3;!R]", "CH3")
 GCPair("[CX4;H2;!R]", "CH2")
 GCPair("[CX3;H0;!R]", "C=")
 GCPair("[OX1;H0;!R]", "O=")
 GCPair("[CX3;H0;!R](=[OX1;H0;!R])", "C=O=")

julia> groups = find_missing_groups_from_smiles("CCCC(=O)CCC"; reduced=true)
3-element Vector{GCPair}:
 GCPair("[CX4;H3;!R]", "CH3")
 GCPair("[CX4;H2;!R]", "CH2")
 GCPair("[CX3;H0;!R](=[OX1;H0;!R])", "C=O=")
```
As we can see, GCIdentifier proposes all the groups that one requires to represent this ketone. However, one aspect of this fragmentation that users may care about is that the methylene (CH$_2$) groups bonded to the ketone group versus those bonded to the methyl group could technically be treated different, due to the difference in environment within which they exist. This distinction can be made by adding the `environment=true` optional argument:
```julia
julia> groups = find_missing_groups_from_smiles("CCCC(=O)CCC"; reduced=true, environment=true)
6-element Vector{GCPair}:
 GCPair("[CX4;H3;!R;$([CX4;H3;!R]([CX4;H2;!R]))]", "CH3(CH2)")
 GCPair("[CX4;H2;!R;$([CX4;H2;!R]([CX4;H3;!R])([CX4;H2;!R]))]", "CH2(CH3CH2)")
 GCPair("[CX4;H2;!R;$([CX4;H2;!R]([CX4;H2;!R])([CX3;H0;!R]))]", "CH2(CH2C=)")
 GCPair("[CX3;H0;!R;$([CX3;H0;!R]([CX4;H2;!R])(=[OX1;H0;!R])([CX4;H2;!R]))](=[OX1;H0;!R;$([OX1;H0;!R](=[CX3;H0;!R]))])", "C=(CH2O=CH2)O=(C=)")
 GCPair("[CX4;H2;!R;$([CX4;H2;!R]([CX3;H0;!R])([CX4;H2;!R]))]", "CH2(C=CH2)")
 GCPair("[CX4;H2;!R;$([CX4;H2;!R]([CX4;H2;!R])([CX4;H3;!R]))]", "CH2(CH2CH3)")
```
As we can see in the above, GCIdentifier now proposes many more groups as we now care about the environment within which each group exists. 

One last flexible element of the `find_missing_groups_from_smiles` function is related to the size of the groups. In all of the above cases, the groups proposed only involve one or two heavy atoms. This is fine for most small molecules. However, for larger ones, such as adenylic acid, the groups proposed may not necessarily be the best representation:
```julia
julia> groups = find_missing_groups_from_smiles("C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)O)O)O)N"; reduced=true)
12-element Vector{GCPair}:
 GCPair("[cX3;H1;R][nX2;H0;R]", "aCHaN")
 GCPair("[cX3;H0;R][nX2;H0;R]", "aCaN")
 GCPair("[cX3;H0;R][NX3;H2;!R]", "aCNH2")
 GCPair("[cX3;H0;R][nX3;H0;R]", "aCaN")
 GCPair("[cX3;H1;R][nX3;H0;R]", "aCHaN")
 GCPair("[CX4;H1;R][nX3;H0;R]", "cCHaN")
 GCPair("[CX4;H1;R][OX2;H0;R]", "cCHcO")
 GCPair("[CX4;H1;R][OX2;H1;!R]", "cCHOH")
 GCPair("[CX4;H1;R][CX4;H2;!R]", "cCHCH2")
 GCPair("[OX2;H0;!R]([PX4;H0;!R])", "OP")
 GCPair("[OX1;H0;!R]", "O=")
 GCPair("[OX2;H1;!R]", "OH")
```
Most groups here are quite reasonable, with the exception of the "O=" and "OH" groups as those will be bonded directly to the phosphorous atom, which we would expect results in very different "O=" and "OH" groups than you'd find on, for example, alcohols. While we could again use the `environment` capability of GCIdentifier to simply specify the environment in which these groups exist, this will result in far more groups being proposed. Ideally, we would like to combine the phosphate group into one. This can be done by specifying a larger `max_group_size` in the optional arguments:
```julia
julia> groups = find_missing_groups_from_smiles("C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)O)O)O)N"; reduced=true, max_group_size=5)
10-element Vector{GCPair}:
 GCPair("[cX3;H1;R][nX2;H0;R]", "aCHaN")
 GCPair("[cX3;H0;R][nX2;H0;R]", "aCaN")
 GCPair("[cX3;H0;R][NX3;H2;!R]", "aCNH2")
 GCPair("[cX3;H0;R][nX3;H0;R]", "aCaN")
 GCPair("[cX3;H1;R][nX3;H0;R]", "aCHaN")
 GCPair("[CX4;H1;R][nX3;H0;R]", "cCHaN")
 GCPair("[CX4;H1;R][OX2;H0;R]", "cCHcO")
 GCPair("[CX4;H1;R][OX2;H1;!R]", "cCHOH")
 GCPair("[CX4;H1;R][CX4;H2;!R]", "cCHCH2")
 GCPair("[PX4;H0;!R]([OX2;H0;!R])(=[OX1;H0;!R])([OX2;H1;!R])([OX2;H1;!R])", "POO=OHOH")
```
This is now a very reasonable representation of adenylic acid.