using GCIdentifier
using Test
using GCIdentifier: @gcstring_str

@testset "All tests" begin
    include("test_group_search.jl")
    include("test_missing_groups.jl")
    include("test_ext_clapeyron.jl")
    include("test_ext_chemicalidentifiers.jl")
end