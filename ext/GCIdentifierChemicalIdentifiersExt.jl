module GCIdentifierChemicalIdentifiersExt

if !isdefined(Base,:get_extension)
    using ChemicalIdentifiers
    using GCIdentifier
else
    using ..ChemicalIdentifiers
    using ..GCIdentifier
end

const GC = GCIdentifier

function GC.get_groups_from_name(component::String,groups::Array{String};connectivity=false)
    res = search_chemical(component)
    smiles = res.smiles
    if connectivity == true
        (smiles,groups_found,connectivity) = get_groups_from_smiles(smiles,groups;connectivity=connectivity)
        return (component,groups_found,connectivity)
    else
        (smiles,groups_found) = get_groups_from_smiles(smiles,groups;connectivity=connectivity)
        return (component,groups_found)
    end
end

end #module
