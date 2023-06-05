
"""
    GCPair

Struct used to hold a description of a group. contains the SMARTS string necessary to match the group within a SMILES query, and the assigned name.

"""
struct GCPair
    smarts::String
    name::String
end

smarts(x::GCPair) = x.smarts
name(x::GCPair) = x.name

#sorting comparison between 2 smatches
function _isless_smatch(smatch1,smatch2)
    #fallback, if one is not matched, throw to the end
    len_smatch1 = length(smatch1)
    len_smatch2 = length(smatch1)
    if len_smatch1 == 0 || len_smatch2 == 0
        return len_smatch1 < len_smatch2
    end

    #first criteria: bigger groups go first,
    #the group with the bigger group
    atom_size1 = length(smatch1[1]["atoms"])
    atom_size2 = length(smatch2[1]["atoms"])
    if atom_size1 != atom_size2
        return atom_size1 < atom_size2
    end

    #second criteria
    #if the size of the match is the same,
    #return the one with the least amount of matches first
    if len_smatch1 != len_smatch2
        len_smatch1 > len_smatch2
    end

    #third criteria
    #return the match with more bonds
    bond_count1 = length(smatch1[1]["bonds"])
    bond_count2 = length(smatch2[1]["bonds"])
    if bond_count1 != bond_count2
        return bond_count1 < bond_count2
    end

    #no more comparizons?
    return false
end

"""
    get_grouplist(x)

Should return a `Vector{GCPair}` containing the available groups for SMILES matching.

"""
function get_grouplist end
get_grouplist(x::Vector{GCPair}) = x

"""
    get_groups_from_smiles(smiles::String,groups;connectivity = false)

Given a SMILES string and a group list (`groups::Vector{GCPair}`), returns a list of groups and their corresponding amount.

If `connectivity` is true, then it will additionally return a vector containing the amount of bonds between each pair.

## Example

```julia
julia> get_groups_from_smiles("CCO",UNIFACGroups)
("CCO", ["CH3" => 1, "CH2" => 1, "OH(P)" => 1])

julia> get_groups_from_smiles("CCO",JobackGroups,connectivity = true)
("CCO", ["-CH3" => 1, "-CH2-" => 1, "-OH (alcohol)" => 1], [("-CH3", "-CH2-") => 1, ("-CH2-", "-OH (alcohol)") => 1])
```
"""
function get_groups_from_smiles(smiles::String,groups;connectivity = false)
    groups = get_grouplist(groups)
    return get_groups_from_smiles(smiles,groups;connectivity = connectivity)
end

function get_groups_from_smiles(smiles::String,groups::Vector{GCPair},lib =DEFAULTLIB;connectivity=false,check = true)
    mol = get_mol(lib,smiles)
    __bonds = __getbondlist(lib,mol)
    group_list = Int[]
    group_id = Int[]
    group_occ_list = Int[]
    atoms_list = Int[]
    coverage_atoms = Vector{Int}[]
    coverage_bonds = Vector{Int}[]

    smatches = []
    smatches_idx = Int[]

    #step 0.a, find all groups that could get a match
    for i in 1:length(groups)
        query_i = get_qmol(lib,smarts(groups[i]))
        if has_substruct_match(lib,mol,query_i)
            push!(smatches,get_substruct_matches(lib,mol,query_i,__bonds))
            push!(smatches_idx,i)
        end
    end

    #step 0.b, sort the matches by the amount of matched atoms. biggest groups come first.
    perm = sortperm(smatches,lt = _isless_smatch,rev = true)
    smatches = smatches[perm]
    smatches_idx = smatches_idx[perm]
    #TODO: identify overlapped groups here?

    for (idx,smatch) in pairs(smatches)
        i = smatches_idx[idx]
        for j in 1:length(smatch)
            #add first match.
            if isempty(atoms_list)
                push!(group_list,i)
                push!(group_id,i)
                push!(group_occ_list,1)
                append!(atoms_list,smatch[j]["atoms"])
                push!(coverage_atoms,copy(smatch[j]["atoms"]))
                push!(coverage_bonds,copy(smatch[j]["bonds"]))
                continue #go to next iteration
            end

            # If no atoms covered by this group are already covered by other groups
            atoms_covered_j = count(in(atoms_list),smatch[j]["atoms"])
            if atoms_covered_j == 0
                append!(atoms_list,smatch[j]["atoms"])
                if !(i in group_list)
                    push!(group_list,i)
                    push!(group_id,i)
                    push!(group_occ_list,1)
                    push!(coverage_atoms,copy(smatch[j]["atoms"]))
                    push!(coverage_bonds,copy(smatch[j]["bonds"]))
                else
                    group_occ_list[end] += 1
                    append!(coverage_atoms[end],smatch[j]["atoms"])
                    append!(coverage_bonds[end],smatch[j]["bonds"])
                end
                continue
            end

            # Check which groups group i has an overlap with
            id = 0
            ng_rm = 0
            for k in 1:length(group_id)
                id += 1

                # Does group 1 cover any atoms of group id
                in_coverage_j = count(in(coverage_atoms[id]), smatch[j]["atoms"])
                in_coverage_j == 0 && continue
                
                # We only care if group i covers _more_ atoms than group k
                length(smatch[j]["atoms"]) < length(coverage_atoms[id]) && continue

                if ((length(smatch[j]["atoms"])>length(coverage_atoms[id])) &
                    # Also make sure that group i covers all the atoms of group k
                    (in_coverage_j == length(coverage_atoms[id]))) |
                    (length(smatch[j]["bonds"])>length(coverage_bonds[id]))
                    # find out which atoms are covered
                    overlap_atoms = coverage_atoms[id][coverage_atoms[id] .∈ [smatch[j]["atoms"]]]
                    id_rm = group_id[id]
                    name_rm = group_list[id]
                    #bond_rm =  coverage_bonds[id][coverage_bonds[id] .∈ [smatch[j]["bonds"]]]
                    filter!(e -> e ∉ overlap_atoms,atoms_list)
                    filter!(e -> e ∉ overlap_atoms,coverage_atoms[id])
                    group_occ_list[id] -= 1

                    # If group k no longer covers any atoms, remove it
                    if group_occ_list[id] == 0
                        filter!(e-> e ≠ 0,group_occ_list)
                        filter!(e-> !isempty(e),coverage_atoms)
                        deleteat!(coverage_bonds,id)
                        filter!(e-> e ≠ id_rm,group_id)
                        filter!(e-> e ≠ name_rm,group_list)
                        id -= 1
                    end
                    ng_rm += 1
                end
            end

            ng_rm == 0 && continue

            if !(i in group_list)
                push!(group_list,i)
                push!(group_id,i)
                push!(group_occ_list,1)
                append!(coverage_atoms,[smatch[j]["atoms"]])
                append!(coverage_bonds,[smatch[j]["bonds"]])
                append!(atoms_list,smatch[j]["atoms"])
            else
                group_occ_list[end] += 1
                append!(atoms_list,smatch[j]["atoms"])
                append!(coverage_atoms[end],smatch[j]["atoms"])
                append!(coverage_bonds[end],smatch[j]["bonds"])
            end
        end
    end

    if check
        atoms = get_atoms(lib,mol)
        if !(count(in(atoms),atoms_list)==length(atoms))
            error("Could not find all groups for "*smiles)
        end
    end

    gcpairs = [name(groups[group_id[i]]) => group_occ_list[i] for i in 1:length(group_id)]

    if connectivity
        return (smiles,gcpairs,get_connectivity(mol,group_id,groups,lib))
    else
        return (smiles,gcpairs)
    end
end

function get_connectivity(mol,group_id,groups,lib = DEFAULTLIB)

    ngroups = length(group_id)
    A = zeros(ngroups,ngroups)
    connectivity = Pair{NTuple{2,String},Int}[]
    for i in 1:ngroups
        gci = groups[group_id[i]]
        smart1 = smarts(gci)
        smart2 = smarts(gci)
        querie = get_qmol(lib,smart1*smart2)
        smatch = get_substruct_matches(lib,mol,querie)
        name_i = name(gci)
        A[i,i] = length(smatch)
        if A[i,i]!=0
            append!(connectivity,[(name_i,name_i)=>Int(A[i,i])])
        end

        for j in i+1:ngroups
            gcj = groups[group_id[j]]
            smart2 = smarts(gcj)
            querie = get_qmol(lib,smart1*smart2)
            smatch = get_substruct_matches(lib,mol,querie)
            A[i,j] = length(smatch)
            name_j = name(gcj)
            if A[i,j]!=0
                append!(connectivity,[(name_i,name_j)=>Int(A[i,j])])
            end
        end
    end
    return connectivity
end



#TODO: move this to Clapeyron?
"""
    @gcstring_str(str)

    given a string of the form "Group1:n1;Group2:2", returns ["Group1" => n1,"Group2 => n2]

"""
macro gcstring_str(str)
    gcpairs = split(str,';')
    res = Pair{String,Int}[]
    for gci in gcpairs
        gc,_ni = split_2(gci,':')
        ni = parse(Int,_ni)
        push!(res,gc => ni)
    end
    res
end

export get_groups_from_name, get_groups_from_smiles