module GCIdentifierClapeyronExt

if !isdefined(Base,:get_extension)
    using Clapeyron
    using GCIdentifier
else
    using ..Clapeyron
    using ..GCIdentifier
end

const GC = GCIdentifier

GC.get_grouplist(m::Clapeyron.EoSModel) = GC.get_grouplist(typeof(m))
GC.get_grouplist(Type{T}) where T <: UNIFAC = GC.UNIFACGroups
GC.get_grouplist(Type{T}) where T <: SAFTgammaMie = GC.SAFTgammaMieGroups
GC.get_grouplist(Type{T}) where T <: Joback = GC.JobackGroups
GC.get_grouplist(Type{T}) where T <: gcPCSAFT = GC.gcPCSAFTgroups

end #module
