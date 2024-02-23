module GCIdentifier
import RDKitMinimalLib
using Combinatorics

@static if !isdefined(Base,:eachsplit)
    eachsplit(str::AbstractString, dlm; limit::Integer=0, keepempty::Bool=true) = split(str,dlm;limit,keepempty)
    eachsplit(str::AbstractString; limit::Integer=0, keepempty::Bool=false)  = split(str;limit,keepempty)
end

split_2(str) = NTuple{2}(eachsplit(str, limit=2))
split_2(str,dlm) = NTuple{2}(eachsplit(str,dlm, limit=2))

#TODO: windows support with MolecularGraph
import MolecularGraph

include("prelude.jl")
include("group_search.jl")
include("missing_groups.jl")
include("database/database.jl")

function get_groups_from_name end #overload this if ChemicalIdentifiers is loaded.

if !isdefined(Base,:get_extension)
    using Clapeyron,ChemicalIdentifiers
    include("../ext/GCIdentifierClapeyronExt.jl")
    include("../ext/GCIdentifierChemicalIdentifiersExt.jl")
end

end # module
