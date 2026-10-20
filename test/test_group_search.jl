function test_gcmatch(groups,smiles,result)
    evaluated = Set(get_groups_from_smiles(smiles,groups,check = false)[2])
    reference = Set(result)
    @test isequal(evaluated,reference)
end

test_gcmatch(groups) = (smiles,result) -> test_gcmatch(groups,smiles,result)

@testset "UNIFAC examples" begin
    #http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html
    unifac = test_gcmatch(UNIFACGroups)

    #alkane group
    unifac("CC",gcstring"CH3:2") #ethane
    unifac("CCCC",gcstring"CH3:2;CH2:2") #n-butane
    unifac("CC(C)C",gcstring"CH3:3;CH:1") #isobutane
    unifac("CC(C)(C)C",gcstring"CH3:4;C:1") #neopentane

    #alpha-olefin group
    unifac("CCCCC=C",gcstring"CH3:1;CH2:3;CH2=CH:1") #hexene-1
    unifac("CCCC=CC",gcstring"CH3:2;CH2:2;CH=CH:1") #hexene-2
    unifac("CCC(C)=C",gcstring"CH3:2;CH2:1;CH2=C:1") #2-methyl-1-butene
    unifac("CC=C(C)C",gcstring"CH3:3;CH=C:1") #2-methyl-2-butene
    unifac("CC(=C(C)C)C",gcstring"CH3:4;C=C:1") #2,3-dimethylbutene

    #aromatic carbon
    unifac("c1c2ccccc2ccc1",gcstring"ACH:8;AC:2") #napthaline
    unifac("c1ccccc1C=C",gcstring"ACH:5;AC:1;CH2=CH:1") #styrene

    #aromatic carbon-alkane
    unifac("Cc1ccccc1",gcstring"ACH:5;ACCH3:1") #toluene
    unifac("CCc1ccccc1",gcstring"ACH:5;ACCH2:1;CH3:1") #ethylbenzene
    unifac("CC(C)c1ccccc1",gcstring"ACH:5;ACCH:1;CH3:2") #cumene

    #alcohol
    unifac("CC(O)C",gcstring"CH3:2;CH:1;OH (S):1") #2-propanol
    unifac("CCO",gcstring"CH3:1;CH2:1;OH (P):1") #ethanol

    #methanol
    unifac("CO",gcstring"CH3OH:1") #methanol

    #water
    unifac("O",gcstring"H2O:1") #water

    #aromatic carbon-alcohol
    unifac("Oc1ccccc1",gcstring"ACH:5;ACOH:1") #phenol

    #carbonyl
    unifac("O=C(C)CC",gcstring"CH3:1;CH2:1;CH3CO:1") #butanone
    unifac("O=C(CC)CC",gcstring"CH3:2;CH2:1;CH2CO:1") #pentanone-3

    #aldehyde
    unifac("CCC=O",gcstring"CH3:1;CH2:1;HCO:1") #propionaldehyde , fails at the matching stage

    #acetate group
    unifac("CCCCOC(=O)C",gcstring"CH3:1;CH2:3;CH3COO:1") #Butyl acetate #fails at the matching stage
    unifac("O=C(OC)CC",gcstring"CH3:2;CH2COO:1") #methyl propionate #fails at the matching stage

    #formate group
    unifac("O=COCC",gcstring"CH3:1;CH2:1;HCOO:1") #ethyl formate

    #ether
    unifac("COC",gcstring"CH3:1;CH3O:1") #dimethyl ether
    unifac("CCOCC",gcstring"CH3:2;CH2:1;CH2O:1") #diethyl ether
    unifac("O(C(C)C)C(C)C",gcstring"CH3:4;CH:1;CHO:1") #diisopropyl ether
    unifac("C1CCOC1",gcstring"CY-CH2:2;THF:1") #tetrahydrofuran #check?

    #primary amine
    unifac("CN",gcstring"CH3NH2:1") #methylamine
    unifac("CCN",gcstring"CH3:1;CH2NH2:1") #ethylamine
    unifac("CC(C)N",gcstring"CH3:2;CHNH2:1") #isopropyl amine

    #secondary amine group
    unifac("CNC",gcstring"CH3:1;CH3NH:1") #dimethylamine
    unifac("CCNCC",gcstring"CH3:2;CH2:1;CH2NH:1") #diethylamine
    unifac("CC(C)NC(C)C",gcstring"CH3:4;CH:1;CHNH:1") #diisopropyl amine

    #tertiary amine
    unifac("CN(C)C",gcstring"CH3:2;CH3N:1") #trimethylamine
    unifac("CCN(CC)CC",gcstring"CH3:3;CH2:2;CH2N:1") #triethylamine

    #aromatic amine
    unifac("c1ccc(cc1)N",gcstring"ACH:5;ACNH2:1") #aniline

    #furfural
    unifac("c1cc(oc1)C=O",gcstring"FURFURAL:1") #furfural

    #carboxylic acids
    unifac("CC(=O)O",gcstring"CH3:1;COOH:1") #acetic acid

    #cyclic ethers
    unifac("C1OCOCO1",gcstring"TRIOXAN:3")
    unifac("C1OCOCC1",gcstring"CY-CH2:1;CY-CH2O:2")
    unifac("C1OCCC1",gcstring"CY-CH2:2;THF:1")

    #non-unique group assignment
    unifac("c1ccccc1COCCOCC",gcstring"ACH:5;ACCH2:1;CH2O:2;CH2:1;CH3:1")
end

@testset "gcPPC-SAFT Connectivity" begin
    # Test a highly branched hydrocarbon
    smiles = "c1ccccc1C(CCCC)(CCCC(C))C"
    (component, groups, connectivity) = get_groups_from_smiles(smiles, gcPPCSAFTGroups; connectivity=true)

    @test isequal(["CH2" => 7, "aCH" => 5, "CH3" => 3, "aC" => 1, "C" => 1] |> Set,Set(groups))
    @test isequal([("aCH", "aCH") => 4, ("CH2", "C") => 2, ("C", "aC") => 1, ("aCH", "aC") => 2, ("CH3", "CH2") => 2, ("CH2", "CH2") => 5, ("CH3", "C") => 1] |> Set,Set(connectivity))

    # Test case where order matters
    smiles = "CC(=O)OC"
    (component, groups, connectivity) = get_groups_from_smiles(smiles, gcPPCSAFTGroups; connectivity=true)

    @test isequal(["CH3" => 2, "COO" => 1] |> Set,Set(groups))
    @test isequal([("CH3", "COO") => 2] |> Set,Set(connectivity))

    # Test case where there is overlap between groups
    smiles = "C(COCCOCCOCCOCCO)O"
    (component, groups, connectivity) = get_groups_from_smiles(smiles, gcPPCSAFTGroups; connectivity=true)

    @test isequal(["CH2" => 6, "OH" => 2, "OCH2" => 4] |> Set,Set(groups))
    @test isequal([("CH2", "CH2") => 2, ("CH2", "OH") => 2, ("CH2", "OCH2") => 6, ("OCH2", "OCH2") => 1] |> Set,Set(connectivity))

    @testset "cases" begin
        function test_assign(groups, smiles, result)
            @test sort(get_groups_from_smiles(smiles, groups)[2]) == sort(result)
        end

        test_assign(groups) = (smiles, result) -> test_assign(groups, smiles, result)
        
        # Issues
        ogunifac = test_assign(ogUNIFACGroups)
        modunifac = test_assign(UNIFACGroups)

        # groups with the same name" begin
        dup = test_assign([GCPair(raw"[CX4;H3]","CH3"), GCPair(raw"[CX4;H2]","CH2"), GCPair(raw"[CX4;H1]","CH"), GCPair(raw"[OH1;$(O[CH2])]","OH"), GCPair(raw"[OH1;$(O[CH1])]","OH")])
        dup(raw"OCC(O)CO", gcstring"CH:1;CH2:2;OH:3")
    
        # stereo SMILES" begin
        modunifac(raw"CC[C@@H](C)O", gcstring"CH:1;CH2:1;CH3:2;OH (S):1")
    
        # exact cover of overlapping groups" begin
        ogunifac(raw"CCOCCl", gcstring"CH2CL:1;CH2O:1;CH3:1")
    
        # ogUNIFAC ring carbons" begin
        ogunifac(raw"C1CCCCC1", gcstring"CH2:6")
    
        # aldehydes" begin
        ogunifac(raw"O=Cc1ccccc1", gcstring"AC:1;ACH:5;HCO:1")
        modunifac(raw"O=Cc1ccccc1", gcstring"AC:1;ACH:5;HCO:1")
    
        # carboxylic acids" begin
        ogunifac(raw"CC(C)CC(=O)O", gcstring"CH:1;CH2:1;CH3:2;COOH:1")
        modunifac(raw"CC(C)CC(=O)O", gcstring"CH:1;CH2:1;CH3:2;COOH:1")
    
        # ogUNIFAC without modUNIFAC groups" begin
        ogunifac(raw"C1COCO1", gcstring"CH2:1;THF:2")
    
        # ogUNIFAC pyridines" begin
        ogunifac(raw"Cc1cccnc1", gcstring"C5H4N:1;CH3:1")
    
        # ogUNIFAC thiophenes" begin
        ogunifac(raw"c1ccsc1", gcstring"C4H4S:1")
    
        # ogUNIFAC cyclic ethers" begin
        ogunifac(raw"CC1CCC(C)O1", gcstring"CH:1;CH2:2;CH3:2;CHO:1")
    
        # ogUNIFAC siloxanes" begin
        ogunifac(raw"C[Si](C)(C)O[Si](C)(C)C", gcstring"CH3:6;SI:1;SIO:1")
    
        # ogUNIFAC hydroxy ethers" begin
        ogunifac(raw"C=COCCO", gcstring"C2H5O2:1;CH2=CH:1")
    
        # ogUNIFAC chlorofluorocarbons" begin
        ogunifac(raw"CC(C)C(F)Cl", gcstring"CH:1;CH3:2;HCCLF:1")
    
        # ogUNIFAC morpholines" begin
        ogunifac(raw"CC(C)N1CCOCC1", gcstring"CH:1;CH2:2;CH2N:1;CH3:2;THF:1")
    
        # hydroxyl groups not bonded to sp3 carbon" begin
        ogunifac(raw"C[Si](C)(C)O", gcstring"CH3:3;OH:1;SI:1")
        modunifac(raw"O=C1OC(C(O)CO)C(O)=C1O", gcstring"C=C:1;CH:1;CH2:1;COO:1;CY-CH:1;OH (P):1;OH (S):1;OH (T):2")
    
        # ogUNIFAC imidazolium" begin
        ogunifac(raw"CC[n+]1ccn(C)c1.O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F", gcstring"BTI:1;CH2:1;CH3:2;IMIDAZOL:1")
    
        # sulfones" begin
        ogunifac(raw"CCCS(=O)(=O)CCC", gcstring"(CH2)2SU:1;CH2:2;CH3:2")
        modunifac(raw"CCCS(=O)(=O)CCC", gcstring"(CH2)2SU:1;CH2:2;CH3:2")
    
        # cyclic sulfides" begin
        ogunifac(raw"C1CCSC1", gcstring"CH2:3;CH2S:1")
        modunifac(raw"C1CCSC1", gcstring"CH2S:1;CY-CH2:3")
    
        # modUNIFAC aromatic carbon next to ring" begin
        modunifac(raw"c1ccc2c(c1)CCC2", gcstring"AC:2;ACH:4;CY-CH2:3")
    
        # modUNIFAC cyclic ethers" begin
        modunifac(raw"C[C@H]1O[C@@H](C)O[C@@H](C)O1", gcstring"CH3:3;TRIOXAN:3")
    
        # modUNIFAC epoxides" begin
        modunifac(raw"CC[C@@H]1OC1(C)C", gcstring"CH2:1;CH3:3;COCH:1")
    
        # modUNIFAC lactams" begin
        modunifac(raw"CCN1CCCCCC1=O", gcstring"CH3:1;CY-CH2:5;NEP:1")
    
        # modUNIFAC tertiary amides" begin
        modunifac(raw"CCC(=O)N(CC)CC", gcstring"AM(CH2)2:1;CH2:1;CH3:3")
    
        # modUNIFAC disulfides" begin
        modunifac(raw"CCCSSCCC", gcstring"-S-S-:1;CH2:4;CH3:2")
    
        # modUNIFAC carbonates" begin
        modunifac(raw"CCCCOC(=O)OCCCC", gcstring"(CH2)2CB:1;CH2:4;CH3:2")
    
        # modUNIFAC ionic liquid cations" begin
        modunifac(raw"CCCC[n+]1ccc(C)cc1.F[B-](F)(F)F", gcstring"BF4-:1;C5H4N+:1;CH2:3;CH3:2")
    
        # modUNIFAC alkyl sulfates" begin
        modunifac(raw"COCCOS(=O)(=O)[O-].Cn1cc[n+](C)c1", gcstring"C3H3N2+:1;CH2:2;CH3:2;CH3O:1;SO4:1")
    
        # modUNIFAC group priority" begin
        modunifac(raw"CCCCOCc1ccccc1", gcstring"ACCH2:1;ACH:5;CH2:2;CH2O:1;CH3:1")
    end
end

# @testset "DDB decompositions" begin
#     ogunifac = test_assign(ogUNIFACGroups)
#     modunifac = test_assign(UNIFACGroups)

#     @testset "groups with the same name" begin
#         dup = test_assign([GCPair(raw"[CX4;H3]","CH3"), GCPair(raw"[CX4;H2]","CH2"), GCPair(raw"[CX4;H1]","CH"), GCPair(raw"[OH1;$(O[CH2])]","OH"), GCPair(raw"[OH1;$(O[CH1])]","OH")])
#         dup(raw"OCC(O)CO", gcstring"CH:1;CH2:2;OH:3")
#         dup(raw"CC(O)CO", gcstring"CH:1;CH2:1;CH3:1;OH:2")
#         dup(raw"OCC(O)C(O)CO", gcstring"CH:2;CH2:2;OH:4")
#     end
#     @testset "stereo SMILES" begin
#         modunifac(raw"CC[C@@H](C)O", gcstring"CH:1;CH2:1;CH3:2;OH (S):1")
#         modunifac(raw"C[C@@H]1CCCCC1=O", gcstring"CH2CO:1;CH3:1;CY-CH:1;CY-CH2:3")
#         modunifac(raw"CC[C@@H](C)C(C)C", gcstring"CH:2;CH2:1;CH3:4")
#     end
#     @testset "exact cover of overlapping groups" begin
#         ogunifac(raw"CCOCCl", gcstring"CH2CL:1;CH2O:1;CH3:1")
#         ogunifac(raw"NCC(N)=O", gcstring"CH2NH2:1;CONH2:1")
#         ogunifac(raw"C1OCOCOCO1", gcstring"THF:4")
#     end
#     @testset "ogUNIFAC ring carbons" begin
#         ogunifac(raw"C1CCCCC1", gcstring"CH2:6")
#         ogunifac(raw"CC1CCCCC1", gcstring"CH:1;CH2:5;CH3:1")
#         ogunifac(raw"OC1CCCCC1", gcstring"CH:1;CH2:5;OH:1")
#     end
#     @testset "aldehydes" begin
#         ogunifac(raw"CC=O", gcstring"CH3:1;HCO:1")
#         ogunifac(raw"CCC=O", gcstring"CH2:1;CH3:1;HCO:1")
#         ogunifac(raw"O=Cc1ccccc1", gcstring"AC:1;ACH:5;HCO:1")
#         modunifac(raw"CC=O", gcstring"CH3:1;HCO:1")
#         modunifac(raw"CCC=O", gcstring"CH2:1;CH3:1;HCO:1")
#         modunifac(raw"O=Cc1ccccc1", gcstring"AC:1;ACH:5;HCO:1")
#     end
#     @testset "carboxylic acids" begin
#         ogunifac(raw"CCC(=O)O", gcstring"CH2:1;CH3:1;COOH:1")
#         ogunifac(raw"CC(C)CC(=O)O", gcstring"CH:1;CH2:1;CH3:2;COOH:1")
#         ogunifac(raw"O=CO", gcstring"HCOOH:1")
#         modunifac(raw"CCC(=O)O", gcstring"CH2:1;CH3:1;COOH:1")
#         modunifac(raw"CC(C)CC(=O)O", gcstring"CH:1;CH2:1;CH3:2;COOH:1")
#         modunifac(raw"O=CO", gcstring"HCOOH:1")
#     end
#     @testset "ogUNIFAC without modUNIFAC groups" begin
#         ogunifac(raw"C1COCO1", gcstring"CH2:1;THF:2")
#         ogunifac(raw"C=Cc1ccncc1", gcstring"C5H4N:1;CH2=CH:1")
#         ogunifac(raw"CN1CCCC1=O", gcstring"NMP:1")
#     end
#     @testset "ogUNIFAC pyridines" begin
#         ogunifac(raw"Cc1cccnc1", gcstring"C5H4N:1;CH3:1")
#         ogunifac(raw"Cc1ccnc(C)c1", gcstring"C5H3N:1;CH3:2")
#         ogunifac(raw"Brc1cccnc1", gcstring"BR:1;C5H4N:1")
#     end
#     @testset "ogUNIFAC thiophenes" begin
#         ogunifac(raw"c1ccsc1", gcstring"C4H4S:1")
#         ogunifac(raw"Cc1ccsc1", gcstring"C4H3S:1;CH3:1")
#         ogunifac(raw"c1ccc2sccc2c1", gcstring"ACH:4;C4H2S:1")
#     end
#     @testset "ogUNIFAC cyclic ethers" begin
#         ogunifac(raw"C1CCOC1", gcstring"CH2:3;THF:1")
#         ogunifac(raw"C1COCCO1", gcstring"CH2:2;THF:2")
#         ogunifac(raw"CC1CCC(C)O1", gcstring"CH:1;CH2:2;CH3:2;CHO:1")
#     end
#     @testset "ogUNIFAC siloxanes" begin
#         ogunifac(raw"C[Si](C)(C)O[Si](C)(C)C", gcstring"CH3:6;SI:1;SIO:1")
#         ogunifac(raw"CO[Si](C)(C)C", gcstring"CH3:3;CH3O:1;SI:1")
#         ogunifac(raw"C[Si](C)(C)O[Si](C)(C)O[Si](C)(C)C", gcstring"CH3:8;SI:1;SIO:2")
#     end
#     @testset "ogUNIFAC hydroxy ethers" begin
#         ogunifac(raw"C=COCCO", gcstring"C2H5O2:1;CH2=CH:1")
#         ogunifac(raw"CC(O)COC(C)C", gcstring"C2H4O2:1;CH:1;CH3:3")
#         ogunifac(raw"CCOC(CO)OCC", gcstring"C2H4O2:1;CH2:1;CH2O:1;CH3:2")
#     end
#     @testset "ogUNIFAC chlorofluorocarbons" begin
#         ogunifac(raw"CC(C)C(F)Cl", gcstring"CH:1;CH3:2;HCCLF:1")
#         ogunifac(raw"FC(F)(F)CC(F)(Cl)Cl", gcstring"CCL2F:1;CF3:1;CH2:1")
#         ogunifac(raw"CCCC(F)Cl", gcstring"CH2:2;CH3:1;HCCLF:1")
#     end
#     @testset "ogUNIFAC morpholines" begin
#         ogunifac(raw"CC(C)N1CCOCC1", gcstring"CH:1;CH2:2;CH2N:1;CH3:2;THF:1")
#         ogunifac(raw"OCCN1CCOCC1", gcstring"CH2:4;CH2N:1;OH:1;THF:1")
#         ogunifac(raw"CCN1CCOCC1", gcstring"CH2:3;CH2N:1;CH3:1;THF:1")
#     end
#     @testset "hydroxyl groups not bonded to sp3 carbon" begin
#         ogunifac(raw"C[Si](C)(C)O", gcstring"CH3:3;OH:1;SI:1")
#         ogunifac(raw"OC1=CC=CC1", gcstring"CH2:1;CH=C:1;CH=CH:1;OH:1")
#         ogunifac(raw"CCN(O)CC", gcstring"CH2:1;CH2N:1;CH3:2;OH:1")
#         modunifac(raw"O=C1OC(C(O)CO)C(O)=C1O", gcstring"C=C:1;CH:1;CH2:1;COO:1;CY-CH:1;OH (P):1;OH (S):1;OH (T):2")
#         modunifac(raw"OC1=CC=CC1", gcstring"CH=C:1;CH=CH:1;CY-CH2:1;OH (T):1")
#         modunifac(raw"OC1=CCC=C1", gcstring"CH=C:1;CH=CH:1;CY-CH2:1;OH (T):1")
#     end
#     @testset "ogUNIFAC imidazolium" begin
#         ogunifac(raw"CC[n+]1ccn(C)c1.O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F", gcstring"BTI:1;CH2:1;CH3:2;IMIDAZOL:1")
#         ogunifac(raw"Cn1cc[n+](C)c1.O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F", gcstring"BTI:1;CH3:2;IMIDAZOL:1")
#         ogunifac(raw"CCCC[n+]1ccn(C)c1.O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F", gcstring"BTI:1;CH2:3;CH3:2;IMIDAZOL:1")
#     end
#     @testset "sulfones" begin
#         ogunifac(raw"CCCS(=O)(=O)CCC", gcstring"(CH2)2SU:1;CH2:2;CH3:2")
#         ogunifac(raw"CCCCS(=O)(=O)CCCC", gcstring"(CH2)2SU:1;CH2:4;CH3:2")
#         ogunifac(raw"CCS(=O)(=O)CC", gcstring"(CH2)2SU:1;CH3:2")
#         modunifac(raw"CCCS(=O)(=O)CCC", gcstring"(CH2)2SU:1;CH2:2;CH3:2")
#         modunifac(raw"CCCCS(=O)(=O)CCCC", gcstring"(CH2)2SU:1;CH2:4;CH3:2")
#         modunifac(raw"CCS(=O)(=O)CC", gcstring"(CH2)2SU:1;CH3:2")
#     end
#     @testset "cyclic sulfides" begin
#         ogunifac(raw"C1CCSC1", gcstring"CH2:3;CH2S:1")
#         ogunifac(raw"C1CCSCC1", gcstring"CH2:4;CH2S:1")
#         ogunifac(raw"C1CCC(SC2CCCCC2)CC1", gcstring"CH:1;CH2:10;CHS:1")
#         modunifac(raw"C1CCSC1", gcstring"CH2S:1;CY-CH2:3")
#         modunifac(raw"C1CCSCC1", gcstring"CH2S:1;CY-CH2:4")
#         modunifac(raw"C1CCC(SC2CCCCC2)CC1", gcstring"CHS:1;CY-CH:1;CY-CH2:10")
#     end
#     @testset "modUNIFAC aromatic carbon next to ring" begin
#         modunifac(raw"c1ccc2c(c1)CCC2", gcstring"AC:2;ACH:4;CY-CH2:3")
#         modunifac(raw"c1ccc(C2CCCCC2)cc1", gcstring"AC:1;ACH:5;CY-CH:1;CY-CH2:5")
#         modunifac(raw"CC1CCc2ccccc21", gcstring"AC:2;ACH:4;CH3:1;CY-CH:1;CY-CH2:2")
#     end
#     @testset "modUNIFAC cyclic ethers" begin
#         modunifac(raw"CC1CCCCO1", gcstring"CH3:1;CY-CH2:3;THF:1")
#         modunifac(raw"CC1OCCO1", gcstring"CH3:1;CY-CH2O:2")
#         modunifac(raw"C[C@H]1O[C@@H](C)O[C@@H](C)O1", gcstring"CH3:3;TRIOXAN:3")
#     end
#     @testset "modUNIFAC epoxides" begin
#         modunifac(raw"CC1OC1C", gcstring"CH3:2;HCOCH:1")
#         modunifac(raw"C1CCC2OC2C1", gcstring"CY-CH2:4;HCOCH:1")
#         modunifac(raw"CC[C@@H]1OC1(C)C", gcstring"CH2:1;CH3:3;COCH:1")
#     end
#     @testset "modUNIFAC lactams" begin
#         modunifac(raw"CCN1CCCCCC1=O", gcstring"CH3:1;CY-CH2:5;NEP:1")
#         modunifac(raw"CC(C)N1CCCC1=O", gcstring"CH3:2;CY-CH2:3;NIPP:1")
#         modunifac(raw"CCN1CCCC1=O", gcstring"CH3:1;CY-CH2:3;NEP:1")
#     end
#     @testset "modUNIFAC tertiary amides" begin
#         modunifac(raw"CCC(=O)N(CC)CC", gcstring"AM(CH2)2:1;CH2:1;CH3:3")
#         modunifac(raw"CN(C)C(=O)c1ccccc1", gcstring"AC:1;ACH:5;AM(CH3)2:1")
#         modunifac(raw"CCCCN(CCCC)C(C)=O", gcstring"AM(CH2)2:1;CH2:4;CH3:3")
#     end
#     @testset "modUNIFAC disulfides" begin
#         modunifac(raw"CSSC", gcstring"-S-S-:1;CH3:2")
#         modunifac(raw"CCSSCC", gcstring"-S-S-:1;CH2:2;CH3:2")
#         modunifac(raw"CCCSSCCC", gcstring"-S-S-:1;CH2:4;CH3:2")
#     end
#     @testset "modUNIFAC carbonates" begin
#         modunifac(raw"CCOC(=O)OCC", gcstring"(CH2)2CB:1;CH3:2")
#         modunifac(raw"CCOC(=O)OC", gcstring"CH2CH3CB:1;CH3:1")
#         modunifac(raw"CCCCOC(=O)OCCCC", gcstring"(CH2)2CB:1;CH2:4;CH3:2")
#     end
#     @testset "modUNIFAC ionic liquid cations" begin
#         modunifac(raw"CC[n+]1ccn(C)c1.O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F", gcstring"BTI-:1;C3H3N2+:1;CH2:1;CH3:2")
#         modunifac(raw"CCn1cc[n+](C)c1C.O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F", gcstring"BTI-:1;C3H2N2+:1;CH2:1;CH3:3")
#         modunifac(raw"CCCC[n+]1ccc(C)cc1.F[B-](F)(F)F", gcstring"BF4-:1;C5H4N+:1;CH2:3;CH3:2")
#     end
#     @testset "modUNIFAC alkyl sulfates" begin
#         modunifac(raw"CCOS(=O)(=O)[O-].CC[n+]1ccn(C)c1", gcstring"C3H3N2+:1;CH2:2;CH3:3;SO4:1")
#         modunifac(raw"COS(=O)(=O)[O-].Cn1cc[n+](C)c1", gcstring"C3H3N2+:1;CH3:3;SO4:1")
#         modunifac(raw"COCCOS(=O)(=O)[O-].Cn1cc[n+](C)c1", gcstring"C3H3N2+:1;CH2:2;CH3:2;CH3O:1;SO4:1")
#     end
#     @testset "modUNIFAC group priority" begin
#         modunifac(raw"CCOCc1ccccc1", gcstring"ACCH2:1;ACH:5;CH2O:1;CH3:1")
#         modunifac(raw"CC(C)OC(C)c1ccccc1", gcstring"ACCH:1;ACH:5;CH3:3;CHO:1")
#         modunifac(raw"CCCCOCc1ccccc1", gcstring"ACCH2:1;ACH:5;CH2:2;CH2O:1;CH3:1")
#     end
# end
