#we use MolecularGraph for windows, but we need the same api.

"""
    RDKitLib

Struct used to select the RDKit library (via `RDKitMinimalLib.jl` package) to perform the group search. Default in Linux and Mac.
"""
struct RDKitLib end

"""
    RDKitLib

Struct used to select the `MolecularGraph.jl` library to perform the group search. Default in Windows.
"""
struct MolecularGraphJL end

#get useful structure from SMILES
function get_mol(::RDKitLib,smiles)
    mol = RDKitMinimalLib.get_mol(smiles)
    return mol
end

function get_mol(::MolecularGraphJL,smiles)
    mol = MolecularGraph.smilestomol(smiles)
    return mol
end

#get useful structure from SMARTS
function get_qmol(::RDKitLib,smarts)
    return RDKitMinimalLib.get_qmol(smarts)
end

function get_qmol(::MolecularGraphJL,smarts)
    return MolecularGraph.smartstomol(smarts)
end


#used in the check function
function get_atoms(::RDKitLib,mol)
    0:(length(RDKitMinimalLib.get_substruct_matches(mol,mol)[1]["atoms"]) - 1)
end

function get_atoms(::MolecularGraphJL,mol)
    0:(length(mol.graph.fadjlist) - 1)
end


#check if query has an substruct match
function has_substruct_match(::RDKitLib,mol,query)
    !isempty(RDKitMinimalLib.get_substruct_match(mol,query))
end

function has_substruct_match(::MolecularGraphJL,mol,query)
    MolecularGraph.has_substruct_match(mol,query)
end

#function used for MolecularGraph.jl. the RDKitLib version returns nothing

__getbondlist(::RDKitLib,mol) = nothing

function __getbondlist(::MolecularGraphJL,mol)
    res = Set{NTuple{2,Int}}()
    adjlist = mol.graph.fadjlist
    for (i,xi) in pairs(adjlist)
        for j in xi
            push!(res,minmax(i,j))
        end
    end
    rvec = sort!(collect(res))
end
#parse substruct matches, use RDKitLib format.
function get_substruct_matches(::RDKitLib,mol,query,__bonds = nothing)
    matches = RDKitMinimalLib.get_substruct_matches(mol,query)
    #stabilize type of result
    res = Dict{String,Vector{Int}}[]
    for match in matches
        dictᵢ = Dict{String,Vector{Int}}()
        dictᵢ["atoms"] = convert(Vector{Int},match["atoms"])
        dictᵢ["bonds"] = convert(Vector{Int},match["bonds"])
        push!(res,dictᵢ)
    end
    return res
end

function get_substruct_matches(lib::MolecularGraphJL,mol,query,bonds = __getbondlist(lib,mol))
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
#select default molecule engine

@static if Sys.iswindows()
    const DEFAULTLIB = MolecularGraphJL()
else
    const DEFAULTLIB = RDKitLib()
end

export RDKitLib, MolecularGraphJL