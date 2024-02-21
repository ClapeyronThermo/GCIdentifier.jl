"""
    find_missing_groups_from_smiles(smiles::String, groups, lib = DEFAULTLIB;max_group_size = nothing, environment=false, reduced=false)

Given a SMILES string and a group list (`groups::Vector{GCPair}`), returns a list of potential groups (`new_groups::Vector{GCPair}`) which could cover those atoms not covered within `groups`. If no `groups` vector is provided, it will simply generate all possible groups for the molecule.

A set of heuristics are built into the code when it comes to combining heavy atoms into large groups:
1. If a carbon atom is bonded to another carbon atom, unless only one of the carbons is on a ring, they will not be combined into a group.
2. All other combinations of atoms are allowed.
The logic behind the first heuristic is due to the fact that neighbouring atoms with similar electronegativities won't have a great impact on each other's properties. As such, they are not combined into a group. In the future, this approach could be extended to use HNMR data to determine which atoms can be combined into the same group.

Optional arguments:
- `max_group_size::Int`: The maximum number of atoms within a group to be generated. If `nothing`, the maximum size is however many atoms a central atom is bonded to.
- `environment::Bool`: If true, the groups SMARTS will include information about the environment of the group is in. For example, in pentane, if environment is false, there will only be one CH2 group, whereas, if environment is true, there will be two CH2 groups, one bonded to CH3 and one bonded to another CH2.
- `reduced::Bool`: If true, the groups will be generated such that the minimum number of groups required to represent the molecule, based on `max_group_size`, will be generated. If false, all possible groups will be generated.

## Example

```julia
julia> find_missing_groups_from_smiles("CC(=O)O")
7-element Vector{GCIdentifier.GCPair}:
 GCIdentifier.GCPair("[CX4;H3;!R]", "CH3")
 GCIdentifier.GCPair("[CX3;H0;!R]", "C=")
 GCIdentifier.GCPair("[OX1;H0;!R]", "O=")
 GCIdentifier.GCPair("[OX2;H1;!R]", "OH")
 GCIdentifier.GCPair("[CX3;H0;!R](=[OX1;H0;!R])", "C=O=")
 GCIdentifier.GCPair("[CX3;H0;!R]([OX2;H1;!R])", "C=OH")
 GCIdentifier.GCPair("[CX3;H0;!R](=[OX1;H0;!R])([OX2;H1;!R])", "C=O=OH")
```
"""
function find_missing_groups_from_smiles(smiles, groups=nothing, lib=MolecularGraphJL(); max_group_size=nothing, environment=false, reduced=false)
    mol = get_mol(lib, smiles)
    __bonds = __getbondlist(lib,mol)
    atoms = get_atoms(lib,mol)

    if isnothing(groups)
        missing_atoms = ones(Bool, length(atoms))
    else
        smatches_idx_expanded, atom_coverage = find_covered_atoms(mol, groups, lib, atoms, __bonds, false)
        missing_atoms = (sum(atom_coverage, dims=1) .== 0)[:]
    end

    atom_type = string.(MolecularGraph.atom_symbol(mol))
    graph = MolecularGraph.to_dict(mol)["graph"]
    bond_mat = zeros(Int, length(atom_type), length(atom_type))
    for i in 1:length(graph)
        bond_mat[graph[i][1], graph[i][2]] = MolecularGraph.props(mol, graph[i][1], graph[i][2]).order
        bond_mat[graph[i][2], graph[i][1]] = MolecularGraph.props(mol, graph[i][1], graph[i][2]).order
    end
    atom_type = atom_type[missing_atoms]
    bond_mat = bond_mat[missing_atoms, missing_atoms]

    is_bonded = bond_mat .> 0

    ring = MolecularGraph.is_in_ring(mol)[missing_atoms]
    ring_string = [if i "" else "!" end for i in ring]

    hydrogens = MolecularGraph.total_hydrogens(mol)[missing_atoms]
    bonds = MolecularGraph.connectivity(mol)[missing_atoms]
    aromatic = MolecularGraph.is_aromatic(mol)[missing_atoms]
    hybrid = MolecularGraph.hybridization(mol)[missing_atoms]
    atom_type = [if aromatic[i] lowercase(atom_type[i]) else atom_type[i] end for i in 1:length(atom_type)]

    natoms = length(atom_type)
    smarts = @. "["*atom_type*"X"*string(bonds)*";H"*string(hydrogens)*";"*ring_string*"R"*"]"
    names = @. generate_group_name(atom_type, bonds, hydrogens, ring, aromatic, hybrid)
    if environment
        new_smarts = deepcopy(smarts)
        new_names = deepcopy(names)
        for i in 1:natoms
            bond_orders = bond_mat[i, is_bonded[i,:]]
            bonded_smarts = smarts[is_bonded[i,:]]
            aromatic_bondable = aromatic[is_bonded[i,:]]
            
            new_smarts[i] = new_smarts[i][1:end-1]*raw";$("*new_smarts[i]*prod([if aromatic_bondable[j] "("*bonded_smarts[j]*")"
                                                                                elseif bond_orders[j]==2 "(="*bonded_smarts[j]*")"
                                                                                elseif bond_orders[j]==3 "(#"*bonded_smarts[j]*")" 
                                                                                else "("*bonded_smarts[j]*")" end for j in 1:length(bonded_smarts)])*")]"
            new_names[i] = new_names[i]*"("*prod(names[is_bonded[i,:]])*")"
        end
        smarts = new_smarts
        names = new_names
    end

    if max_group_size == 1
        unique_smarts = unique(smarts)
        unique_names = []
        for i in 1:length(unique_smarts)
            push!(unique_names, names[findall(x->x==unique_smarts[i], smarts)[1]])
        end
    
        new_groups = [GCPair(unique_smarts[i], unique_names[i]) for i in 1:length(unique_smarts)]

        return new_groups
    end

    occurrence = ones(Int, length(smarts))

    for i in 1:natoms
        # println(occurrence)
        # println(smarts[i])
        # If carbon atom on ring
        if (atom_type[i] == "C" || atom_type[i] == "c") && ring[i] && (occurrence[i] > 0 || reduced)
            for j in 1:natoms
                if is_bonded[i,j]
                    # Bond it to any other atom that is non-carbon atom on a ring or carbon atom not on a ring
                    if (atom_type[j] == "C" && !ring[j] || (atom_type[j] != "C" && atom_type[j] != "c"))
                        push!(smarts, smarts[i]*smarts[j])
                        push!(names, names[i]*names[j])
                        occurrence[i] -= 1
                        occurrence[j] -= 1
                        append!(occurrence, [1])
                    end
                end
            end
        elseif ((atom_type[i] != "C" && atom_type[i] != "c") || (atom_type[i] == "C" && !ring[i])) && (occurrence[i] > 0 || !reduced)
            # Bond it to any other atom that is not a carbon atom
            nbonds = sum(is_bonded[i,:].&occurrence[1:natoms] .> 0)
            bondable_atom_types = atom_type[1:natoms][is_bonded[i,:].&occurrence[1:natoms] .> 0]

            is_carbon = bondable_atom_types .== "C" .|| bondable_atom_types .== "c"
            ncarbons = sum(is_carbon)
            nbonds = nbonds - ncarbons
            idx_bondable = is_bonded[i,:].&occurrence[1:natoms] .> 0

            bondable_smarts = smarts[1:natoms][is_bonded[i,:].&occurrence[1:natoms] .> 0][.!(is_carbon)]
            bondable_names = names[1:natoms][is_bonded[i,:].&occurrence[1:natoms] .> 0][.!(is_carbon)]
            bondable_atom_types = bondable_atom_types[.!(is_carbon)]
            bond_orders = bond_mat[i, is_bonded[i,:].&occurrence[1:natoms] .> 0][.!(is_carbon)]

            if isnothing(max_group_size)
                max_size = nbonds+1
            else
                max_size = max_group_size
            end

            if reduced
                min_size = max_size-1
            else  
                min_size = 1
            end
            # println(smarts[i])

            # println(max_size-1)
            # println(min_size)
            # println(nbonds)

            for k in min_size:max_size-1
                # println(k)
                combs = Combinatorics.combinations(1:nbonds, k)
                for comb in combs
                    # println(bondable_atom_types[comb])
                    if any(bondable_atom_types[comb] .== "C" .|| bondable_atom_types[comb] .== "c")
                        # println(comb)
                        continue
                    else
                        new_smarts = deepcopy(smarts[i])
                        new_names = deepcopy(names[i])
                        for l in 1:length(comb)
                            if bond_orders[comb[l]] == 2
                                new_smarts *= "(="*bondable_smarts[comb[l]]*")"
                            elseif bond_orders[comb[l]] == 3
                                new_smarts *= "(#"*bondable_smarts[comb[l]]*")"
                            else
                                new_smarts *= "("*bondable_smarts[comb[l]]*")"
                            end
                            new_names *= bondable_names[comb[l]]
                        end
                        push!(smarts, new_smarts)
                        push!(names, new_names)
                        occurrence[i] -= 1
                        idx = 1:natoms
                        idx = idx[idx_bondable][.!(is_carbon)][comb]
                        occurrence[idx] .-= 1
                        append!(occurrence, [1])
                    end
                end
            end
        end
    end

    if reduced 
        unique_smarts = unique(smarts[occurrence .> 0])
    else
        unique_smarts = unique(smarts)
    end
    # find the names of the unique smarts
    unique_names = []
    occurrence = zeros(Int, length(unique_smarts))
    for i in 1:length(unique_smarts)
        push!(unique_names, names[findall(x->x==unique_smarts[i], smarts)[1]])
        query_i = get_qmol(lib,unique_smarts[i])

        occurrence[i] = length(get_substruct_matches(lib,mol,query_i,__bonds))

        # println(unique_smarts[i], " ", unique_names[i], " ", occurrence[i])
    end

    new_groups = [GCPair(unique_smarts[i], unique_names[i]) for i in 1:length(unique_smarts)]

    return new_groups
end

function generate_group_name(atom_type::String, bond::Int, hydrogens::Int, ring::Bool, aromatic::Bool, hybrid::Symbol)
    name = ""
    if  !ring # If is not on ring
        if hydrogens == 0
            name = atom_type
        elseif hydrogens == 1
            name = atom_type*"H"
        else
            name = atom_type*"H"*string(hydrogens)
        end
    elseif ring && lowercase(atom_type) != atom_type # If is on ring and not aromatic
        if hydrogens == 0
            name = "c"*atom_type
        elseif hydrogens == 1
            name = "c"*atom_type*"H"
        else
            name = "c"*atom_type*"H"*string(hydrogens)
        end
    elseif lowercase(atom_type) == atom_type # If is aromatic
        if hydrogens == 0
            name = "a"*uppercase(atom_type)
        elseif hydrogens == 1
            name = "a"*uppercase(atom_type)*"H"
        else
            name = "a"*uppercase(atom_type)*"H"*string(hydrogens)
        end
    end

    if hybrid == :sp2 && !aromatic
        name *= "="
    elseif hybrid == :sp
        name *= "#"
    end
    return name
end

export find_missing_groups_from_smiles