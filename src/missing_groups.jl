function find_missing_groups_from_smiles(smiles,groups, lib=MolecularGraphJL(); max_group_size=5, environment=false)
    mol = get_mol(lib, smiles)

    __bonds = __getbondlist(lib,mol)

    atom_type = string.(MolecularGraph.atom_symbol(mol))
    graph = MolecularGraph.to_dict(mol)["graph"]
    bond_mat = zeros(Int, length(atom_type), length(atom_type))
    for i in 1:length(graph)
        bond_mat[graph[i][1], graph[i][2]] = MolecularGraph.props(mol, graph[i][1], graph[i][2]).order
        bond_mat[graph[i][2], graph[i][1]] = MolecularGraph.props(mol, graph[i][1], graph[i][2]).order
    end

    is_bonded = bond_mat .> 0

    ring = MolecularGraph.is_in_ring(mol)
    ring_string = [if i "" else "!" end for i in ring]

    hydrogens = MolecularGraph.total_hydrogens(mol)
    bonds = MolecularGraph.connectivity(mol)
    aromatic = MolecularGraph.is_aromatic(mol)
    hybrid = MolecularGraph.hybridization(mol)
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

    for i in 1:natoms
        for j in 1:natoms
            if is_bonded[i,j]
                # If carbon not on ring
                if atom_type[i] == "C" && !ring[i] 
                    # Bond it to any other atom that is non-carbon atom not on ring or carbon atom on a ring
                    if (atom_type[j] == "C" && ring[j] || atom_type[j] != "C")
                        push!(smarts, smarts[i]*smarts[j])
                        push!(names, names[i]*names[j])
                    end
                # If carbon atom on ring
                elseif (atom_type[i] == "C" || atom_type[i] == "c") && ring[i]
                    # Bond it to any other atom that is non-carbon atom on a ring or carbon atom not on a ring
                    if (atom_type[j] == "C" && !ring[j] || (atom_type[j] != "C" && atom_type[j] != "c"))
                        push!(smarts, smarts[i]*smarts[j])
                        push!(names, names[i]*names[j])
                    end
                elseif (atom_type[i] != "C" && atom_type[i] != "c")
                    # Bond it to any other atom that is not a carbon atom
                    nbonds = sum(is_bonded[i,:])
                    bondable_smarts = smarts[1:natoms][is_bonded[i,:]]
                    bondable_names = names[1:natoms][is_bonded[i,:]]
                    bondable_atom_types = atom_type[1:natoms][is_bonded[i,:]]
                    bond_orders = bond_mat[i, is_bonded[i,:]]
                    for k in 1:max_group_size-1
                        combs = Combinatorics.combinations(1:nbonds, k)
                        for comb in combs
                            if any(bondable_atom_types[comb] .== "C" .|| bondable_atom_types[comb] .== "c")
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
                            end
                        end
                    end
                end
            end
        end
    end

    unique_smarts = unique(smarts)
    # find the names of the unique smarts
    unique_names = []
    occurence = zeros(Int, length(unique_smarts))
    for i in 1:length(unique_smarts)
        push!(unique_names, names[findall(x->x==unique_smarts[i], smarts)[1]])
        query_i = get_qmol(lib,unique_smarts[i])

        occurence[i] = length(get_substruct_matches(lib,mol,query_i,__bonds))

        # println(unique_smarts[i], " ", unique_names[i], " ", occurence[i])
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