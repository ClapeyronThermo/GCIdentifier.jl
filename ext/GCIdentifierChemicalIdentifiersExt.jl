module GCIdentifierChemicalIdentifiersExt

using ChemicalIdentifiers
using GCIdentifier

const GC = GCIdentifier

function GC.get_groups_from_name(component::String,groups;connectivity = false)
    groups = GC.get_grouplist(groups)
    return GC.get_groups_from_name(component,groups,connectivity=connectivity)
end

function GC.get_groups_from_name(component::String,groups::Vector{GC.GCPair},lib = GC.DEFAULTLIB;connectivity=false,check = true)
    res = search_chemical(component)
    smiles = res.smiles
    gcpairs = get_groups_from_smiles(smiles,groups,lib;connectivity=connectivity,check = check)
    if connectivity == true
        (smiles,groups_found,connectivity) = gcpairs
        return (component,groups_found,connectivity)
    else
        (smiles,groups_found) = gcpairs
        return (component,groups_found)
    end
end

end #module
