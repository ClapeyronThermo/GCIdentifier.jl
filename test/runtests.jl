using GCIdentifier
using Test
using GCIdentifier: @gcstring_str

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
    unifac("CC(O)C",gcstring"CH3:2;CH:1;OH(S):1") #2-propanol
    unifac("CCO",gcstring"CH3:1;CH2:1;OH(P):1") #ethanol

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
    #unifac("CCC=O",gcstring"CH3:1;CH2:1;HCO:1") #propionaldehyde , fails at the matching stage

    #acetate group
    #unifac("CCCCOC(=O)C",gcstring"CH3:1;CH2:3;CH3COO:1") #Butyl acetate #fails at the matching stage
    #unifac("O=C(OC)CC",gcstring"CH3:2;CH2COO:1") #methyl propionate #fails at the matching stage
    
    #formate group
    unifac("O=COCC",gcstring"CH3:1;CH2:1;HCOO:1") #ethyl formate

    #ether
    unifac("COC",gcstring"CH3:1;CH3O:1") #dimethyl ether
    unifac("CCOCC",gcstring"CH3:2;CH2:1;CH2O:1") #diethyl ether
    unifac("O(C(C)C)C(C)C",gcstring"CH3:4;CH:1;CHO:1") #diisopropyl ether
    #unifac("C1CCOC1",gcstring"CH2:3;THF:1") #tetrahydrofuran #check?
    
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
    
    #piridine
    unifac("c1ccncc1",gcstring"C5H5N:1")
    
    
    
    #furfural
    unifac("c1cc(oc1)C=O",gcstring"FURFURAL:1") #furfural

end
