# Defining Custom Groups
Within GCIdentifier, we support the following group-contribution methods:
* Joback's Method (`JobackGroups`)
* Original UNIFAC (`ogUNIFACGroups`)
* (Dortmund) UNIFAC (`UNIFACGroups`)
* gcPC-SAFT (`gcPCSAFTGroups`)
* gcPCP-SAFT (`gcPCPSAFTGroups`)
* SAFT-$\gamma$ Mie (`SAFTgammaMieGroups`)
* SAFT-$\gamma$ Mie, generated from the Chemeo database (`SAFTgammaMieChemeoGroups`)

There are many more available that we have yet to implement. If you wish do so yourself, then all that needs to be done is define a vector of `GCPair`s. `GCPair` is a struct contain the group SMARTS and name:
```julia
julia> group = GCPair("[CX4H3]","CH3")

julia> group.name
"CH3"

julia> group.smarts
"[CX4H3]"
```
While the group name is entirely arbitrary, one must be very careful when defining the SMARTS as this is the information used by GCIdentifier (and MolecularGraph) to identify the groups. To learn more about SMARTS, one can consult:
* [Guide to understanding the SMARTS nomenclature](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
* [Tool to visualise SMARTS](https://smarts.plus)

Once the user is certain of their SMARTS representation and made the vector of `GCPair`s, then one can simply feed in this vector (`GroupList`) into `get_groups_from_smiles`:
```julia
julia> GroupList = [GCPair("[CX4H3]","CH3"),GCPair("[CX4H2]","CH2")]

julia> get_groups_from_smiles("CCCC", GroupList)
("CCCC", ["CH3" => 2, "CH2" => 2])
```
If you've defined your own `GroupList` for an existing (or new) group-contribution method, feel free to open a pull request!

