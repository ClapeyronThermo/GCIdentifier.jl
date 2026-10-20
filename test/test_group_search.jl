function test_gcmatch(groups,smiles,result)
    evaluated = Set(get_groups_from_smiles(smiles,groups,check = false)[2])
    reference = Set(result)
    @test isequal(evaluated,reference)
end

test_gcmatch(groups) = (smiles,result) -> test_gcmatch(groups,smiles,result)

@testset "UNIFAC examples" begin
    #http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html
    unifac = test_gcmatch(UNIFACGroups)
    ogunifac = test_gcmatch(ogUNIFACGroups)

    #alkane group
    unifac("CC",gcstring"CH3:2") #ethane
    unifac("CCCC",gcstring"CH3:2;CH2:2") #n-butane
    unifac("CC(C)C",gcstring"CH3:3;CH:1") #isobutane
    unifac("CC(C)(C)C",gcstring"CH3:4;C:1") #neopentane
    ogunifac("C1CCCCC1",gcstring"CH2:6") #cyclohexane

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
    unifac("c1ccc2c(c1)CCC2",gcstring"AC:2;ACH:4;CY-CH2:3") #indane

    #alcohol
    unifac("CC(O)C",gcstring"CH3:2;CH:1;OH (S):1") #2-propanol
    unifac("CCO",gcstring"CH3:1;CH2:1;OH (P):1") #ethanol
    unifac("CC[C@@H](C)O",gcstring"CH:1;CH2:1;CH3:2;OH (S):1") #(R)-2-butanol
    unifac("O=C1OC(C(O)CO)C(O)=C1O",gcstring"C=C:1;CH:1;CH2:1;COO:1;CY-CH:1;OH (P):1;OH (S):1;OH (T):2") #ascorbic acid
    ogunifac("C[Si](C)(C)O",gcstring"CH3:3;OH:1;SI:1") #trimethylsilanol
    test_gcmatch([GCPair(raw"[CX4;H3]","CH3"), GCPair(raw"[CX4;H2]","CH2"), GCPair(raw"[CX4;H1]","CH"), GCPair(raw"[OH1;$(O[CH2])]","OH"), GCPair(raw"[OH1;$(O[CH1])]","OH")])("OCC(O)CO",gcstring"CH:1;CH2:2;OH:3") #glycerol

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
    unifac("O=Cc1ccccc1",gcstring"AC:1;ACH:5;HCO:1") #benzaldehyde
    ogunifac("O=Cc1ccccc1",gcstring"AC:1;ACH:5;HCO:1") #benzaldehyde

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
    ogunifac("CCOCCl",gcstring"CH2CL:1;CH2O:1;CH3:1") #chloromethyl ethyl ether

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
    unifac("CC(C)CC(=O)O",gcstring"CH:1;CH2:1;CH3:2;COOH:1") #isovaleric acid
    ogunifac("CC(C)CC(=O)O",gcstring"CH:1;CH2:1;CH3:2;COOH:1") #isovaleric acid

    #cyclic ethers
    unifac("C1OCOCO1",gcstring"TRIOXAN:3")
    unifac("C1OCOCC1",gcstring"CY-CH2:1;CY-CH2O:2")
    unifac("C1OCCC1",gcstring"CY-CH2:2;THF:1")
    unifac("C[C@H]1O[C@@H](C)O[C@@H](C)O1",gcstring"CH3:3;TRIOXAN:3") #paraldehyde
    ogunifac("C1COCO1",gcstring"CH2:1;THF:2") #1,3-dioxolane
    ogunifac("CC1CCC(C)O1",gcstring"CH:1;CH2:2;CH3:2;CHO:1") #2,5-dimethyltetrahydrofuran

    #epoxides
    unifac("CC[C@@H]1OC1(C)C",gcstring"CH2:1;CH3:3;COCH:1") #3-ethyl-2,2-dimethyloxirane

    #glycol ethers
    ogunifac("C=COCCO",gcstring"C2H5O2:1;CH2=CH:1") #2-(vinyloxy)ethanol

    #morpholines
    ogunifac("CC(C)N1CCOCC1",gcstring"CH:1;CH2:2;CH2N:1;CH3:2;THF:1") #4-isopropylmorpholine

    #pyridines
    ogunifac("Cc1cccnc1",gcstring"C5H4N:1;CH3:1") #3-picoline

    #thiophenes
    ogunifac("c1ccsc1",gcstring"C4H4S:1") #thiophene

    #amides
    unifac("CCC(=O)N(CC)CC",gcstring"AM(CH2)2:1;CH2:1;CH3:3") #N,N-diethylpropionamide
    unifac("CCN1CCCCCC1=O",gcstring"CH3:1;CY-CH2:5;NEP:1") #N-ethylcaprolactam

    #carbonates
    unifac("CCCCOC(=O)OCCCC",gcstring"(CH2)2CB:1;CH2:4;CH3:2") #dibutyl carbonate

    #chlorofluorocarbons
    ogunifac("CC(C)C(F)Cl",gcstring"CH:1;CH3:2;HCCLF:1") #1-chloro-1-fluoro-2-methylpropane

    #sulfides
    unifac("C1CCSC1",gcstring"CH2S:1;CY-CH2:3") #tetrahydrothiophene
    ogunifac("C1CCSC1",gcstring"CH2:3;CH2S:1") #tetrahydrothiophene
    unifac("CCCSSCCC",gcstring"-S-S-:1;CH2:4;CH3:2") #dipropyl disulfide

    #sulfones
    unifac("CCCS(=O)(=O)CCC",gcstring"(CH2)2SU:1;CH2:2;CH3:2") #dipropyl sulfone
    ogunifac("CCCS(=O)(=O)CCC",gcstring"(CH2)2SU:1;CH2:2;CH3:2") #dipropyl sulfone

    #siloxanes
    ogunifac("C[Si](C)(C)O[Si](C)(C)C",gcstring"CH3:6;SI:1;SIO:1") #hexamethyldisiloxane

    #ionic liquids
    unifac("CCCC[n+]1ccc(C)cc1.F[B-](F)(F)F",gcstring"BF4-:1;C5H4N+:1;CH2:3;CH3:2") #1-butyl-4-methylpyridinium tetrafluoroborate
    unifac("COCCOS(=O)(=O)[O-].Cn1cc[n+](C)c1",gcstring"C3H3N2+:1;CH2:2;CH3:2;CH3O:1;SO4:1") #1,3-dimethylimidazolium 2-methoxyethyl sulfate
    ogunifac("CC[n+]1ccn(C)c1.O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F",gcstring"BTI:1;CH2:1;CH3:2;IMIDAZOL:1") #1-ethyl-3-methylimidazolium bis(trifluoromethylsulfonyl)imide

    #non-unique group assignment
    unifac("c1ccccc1COCCOCC",gcstring"ACH:5;ACCH2:1;CH2O:2;CH2:1;CH3:1")
    unifac("CCCCOCc1ccccc1",gcstring"ACCH2:1;ACH:5;CH2:2;CH2O:1;CH3:1") #benzyl butyl ether
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
end
