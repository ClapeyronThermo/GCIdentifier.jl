#we use MolecularGraph for windows, but we need the same api.

function get_mol(smiles)
    mol = MolecularGraph.smilestomol(smiles)
    return mol
end

function get_qmol(smarts)
    return MolecularGraph.smartstomol(smarts)
end

function get_atoms(mol)
    0:(length(mol.graph.fadjlist) - 1)
end

function has_substruct_match(mol,query)
    MolecularGraph.has_substruct_match(mol,query)
end


function __getbondlist(mol)
    res = Set{NTuple{2,Int}}()
    adjlist = mol.graph.fadjlist
    for (i,xi) in pairs(adjlist)
        for j in xi
            push!(res,minmax(i,j))
        end
    end
    rvec = sort!(collect(res))
end

function get_substruct_matches(mol,query,bonds = __getbondlist(mol))
    matches = MolecularGraph.substruct_matches(mol,query)
    res = Dict{String,Vector{Int}}[]
    #convert to RDKitLib expected struct.
    for match in matches
        dictᵢ = Dict{String,Vector{Int}}()
        atomsᵢ = collect(keys(match))
        if length(atomsᵢ) == 1
            bondsᵢ = Int[]
        else
            bondsᵢ = __getbonds_mg(bonds,atomsᵢ)
        end
        atomsᵢ .-= 1
        dictᵢ["atoms"] = atomsᵢ
        dictᵢ["bonds"] = bondsᵢ
        push!(res,dictᵢ)
    end
    return res
end

function __getbonds_mg(list,atoms)
    res_set = Set{Int}()
    n = length(atoms)
    for i in 1:n
        for j in (i+1):n
            ij = minmax(atoms[i],atoms[j])
            bond_idx = findfirst(isequal(ij),list)
            if bond_idx !== nothing
                push!(res_set,bond_idx)
            end
        end
    end
    res = sort!(collect(res_set))
    res .-= 1
    return res
end
