module GCIdentifierClapeyronExt

if isdefined(Base,:get_extension)
    using Clapeyron
    using GCIdentifier
else
    using ..Clapeyron
    using ..GCIdentifier
end

const GC = GCIdentifier

GC.get_grouplist(m::Clapeyron.EoSModel) = GC.get_grouplist(typeof(m))
GC.get_grouplist(::Type{T}) where T <: Clapeyron.UNIFAC = GC.UNIFACGroups
GC.get_grouplist(::Type{T}) where T <: Clapeyron.SAFTgammaMie = GC.SAFTgammaMieGroups
GC.get_grouplist(::Type{T}) where T <: Clapeyron.JobackIdeal = GC.JobackGroups
GC.get_grouplist(::Type{T}) where T <: Clapeyron.gcPCSAFT = GC.gcPCSAFTgroups

end #module
