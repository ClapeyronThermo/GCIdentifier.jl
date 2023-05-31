module GCIdentifier
using RDKitMinimalLib, ChemicalIdentifiers

@static if !isdefined(Base,:eachsplit)
    eachsplit(str::AbstractString, dlm; limit::Integer=0, keepempty::Bool=true) = split(str,dlm;limit,keepempty)
    eachsplit(str::AbstractString; limit::Integer=0, keepempty::Bool=false)  = split(str;limit,keepempty)
end

split_2(str) = NTuple{2}(eachsplit(str, limit=2))
split_2(str,dlm) = NTuple{2}(eachsplit(str,dlm, limit=2))

include("group_search.jl")
include("database/database.jl")

if !isdefined(Base,:get_extension)
    using Clapeyron
    include("../ext/GCIdentifierClapeyronExt.jl")
end
end # module