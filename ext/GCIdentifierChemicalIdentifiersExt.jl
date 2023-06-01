module GCIdentifierChemicalIdentifiersExt

if !isdefined(Base,:get_extension)
    using ChemicalIdentifiers
    using GCIdentifier
else
    using ..ChemicalIdentifiers
    using ..GCIdentifier
end

const GC = GCIdentifier

function GC.get_groups_from_name(component::String,groups::Array{String};connectivity=false,check = true)
    res = search_chemical(component)
    smiles = res.smiles
    gcpairs = get_groups_from_smiles(smiles,groups;connectivity=connectivity,check = check)
    if connectivity == true
        (smiles,groups_found,connectivity) = gcpairs
        return (component,groups_found,connectivity)
    else
        (smiles,groups_found) = gcpairs
        return (component,groups_found)
    end
end

end #module
