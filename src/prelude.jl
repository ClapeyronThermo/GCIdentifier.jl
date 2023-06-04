#we use MolecularGraph for windows, but we need the same api.


struct RDKit end
struct MolecularGraphJL end

#get useful structure from SMILES
function get_mol(::RDKit,smiles)
    mol = RDKitMinimalLib.get_mol(smiles)
    return mol
end

function get_mol(::MolecularGraphJL,smiles)
    mol = MolecularGraph.smilestomol(smiles)
    return mol
end

#get useful structure from SMARTS
function get_qmol(::RDKit,smarts)
    return RDKitMinimalLib.get_qmol(smarts)
end

function get_qmol(::MolecularGraphJL,smarts)
    return MolecularGraph.smartstomol(smarts)
end


#used in the check function
function get_atoms(::RDKit,mol)
    RDKitMinimalLib.get_substruct_matches(mol,mol)[1]["atoms"]
end

function get_atoms(::MolecularGraphJL,mol)
    0:(length(mol.graph.fadjlist) - 1)
end


#check if query has an substruct match
function has_substruct_match(::RDKit,mol,query)
    !isempty(RDKitMinimalLib.get_substruct_match(mol,query))
end

function has_substruct_match(::MolecularGraphJL,mol,query)
    MolecularGraph.has_substruct_match(mol,query)
end

#parse substruct matches, use RDKit format.

function get_substruct_matches(::RDKit,mol,query)
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

function get_substruct_matches(::MolecularGraphJL,mol,query)
    matches = MolecularGraph.substruct_matches(mol,query)
    res = Dict{String,Vector{Int}}[]
    @show collect(matches)
    #convert to RDKit expected struct.
    for match in matches
        dictᵢ = Dict{String,Vector{Int}}()
        atomsᵢ = collect(keys(match))
        atomsᵢ .-= 1

        dictᵢ["atoms"] = atomsᵢ
        dictᵢ["bonds"] = atomsᵢ
        push!(res,dictᵢ)
    end
    return res
end

#select default molecule engine

@static if Sys.iswindows()
    const DEFAULTLIB = MolecularGraphJL()
else
    const DEFAULTLIB = RDKit()
end