using GCIdentifier
using Test
using GCIdentifier: @gcstring_str

function test_gcmatch(groups,smiles,result)
    obtained = Set(get_groups_from_smiles(smiles,groups)[2])
    evaluated = Set(result)
    @test isequal(obtained,evaluated)
end

test_gcmatch(groups) = (smiles,result) -> test_gcmatch(groups,smiles,result)

@testset "UNIFAC" begin
    #http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html
    unifac = test_gcmatch(UNIFACGroups)
    unifac("CC",gcstring"CH3:2")
    unifac("CCCC",gcstring"CH3:2;CH2:2")
    unifac("CC(C)C",gcstring"CH3:3;CH:1")
    unifac("CC(C)(C)C",gcstring"CH3:4;C:1")
end