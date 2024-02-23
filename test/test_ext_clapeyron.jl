using Clapeyron

@testset "Extension - Clapeyron" begin
    smiles = "c1ccccc1C(CCCC)(CCCC(C))C"
    (component, groups) = get_groups_from_smiles(smiles, JobackIdeal)
    @test isequal(groups,["-CH3" => 3, 
                          "-CH2-" => 7, 
                          ">C<" => 1, 
                          "ring=CH-" => 5, 
                          "ring=C<" => 1])
end