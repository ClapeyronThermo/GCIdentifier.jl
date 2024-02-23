using ChemicalIdentifiers

@testset "Extension - Clapeyron" begin
    name = "ibuprofen"
    (component, groups) = get_groups_from_name(name, UNIFACGroups)
    @test isequal(groups, ["COOH" => 1, 
                           "CH3" => 3, 
                           "CH" => 1, 
                           "ACH" => 4, 
                           "ACCH2" => 1, 
                           "ACCH" => 1])
end