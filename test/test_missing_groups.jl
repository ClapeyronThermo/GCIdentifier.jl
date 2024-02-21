@testset "Missing groups functionality" begin
    @testset "Ketones in SAFT-γ Mie" begin
        smiles = "CCC(=O)CC"
        @test_throws "Could not find all groups for " get_groups_from_smiles(smiles,SAFTgammaMieGroups)
        new_groups = find_missing_groups_from_smiles(smiles, SAFTgammaMieGroups)
        @test isequal(new_groups,[GCIdentifier.GCPair("[CX3;H0;!R]", "C="),
                                  GCIdentifier.GCPair("[OX1;H0;!R]", "O="),
                                  GCIdentifier.GCPair("[CX3;H0;!R](=[OX1;H0;!R])", "C=O=")])
        
        new_groups = find_missing_groups_from_smiles(smiles, SAFTgammaMieGroups; reduced = true)
        @test isequal(new_groups,[GCIdentifier.GCPair("[CX3;H0;!R](=[OX1;H0;!R])", "C=O=")])
    end

    @testset "Generic groups for AMP" begin
        smiles = "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)O)O)O)N"
        groups = find_missing_groups_from_smiles(smiles)
        @test isequal(groups[end-1], GCIdentifier.GCPair("[PX4;H0;!R]([OX2;H0;!R])(=[OX1;H0;!R])([OX2;H1;!R])([OX2;H1;!R])", "POO=OHOH"))
        
        groups = find_missing_groups_from_smiles(smiles; reduced = true)
        @test isequal(groups[end-2], GCIdentifier.GCPair("[OX2;H0;!R]([PX4;H0;!R])", "OP"))

        groups = find_missing_groups_from_smiles(smiles; reduced = true, max_group_size = 5)
        @test isequal(groups[end], GCIdentifier.GCPair("[PX4;H0;!R]([OX2;H0;!R])(=[OX1;H0;!R])([OX2;H1;!R])([OX2;H1;!R])", "POO=OHOH"))
    end
    
end
