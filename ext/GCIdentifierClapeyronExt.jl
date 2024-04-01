module GCIdentifierClapeyronExt

using Clapeyron
using GCIdentifier

const GC = GCIdentifier

GC.get_grouplist(m::Clapeyron.EoSModel) = GC.get_grouplist(typeof(m))
GC.get_grouplist(::Type{T}) where T <: Clapeyron.UNIFAC = GC.UNIFACGroups
GC.get_grouplist(::Type{T}) where T <: Clapeyron.SAFTgammaMie = GC.SAFTgammaMieGroups
GC.get_grouplist(::Type{T}) where T <: Clapeyron.JobackIdeal = GC.JobackGroups
GC.get_grouplist(::Type{T}) where T <: Clapeyron.HomogcPCPSAFT = GC.gcPCSAFTGroups
GC.get_grouplist(::Type{T}) where T <: Clapeyron.gcPPCSAFT = GC.gcPPCSAFTGroups

end #module
