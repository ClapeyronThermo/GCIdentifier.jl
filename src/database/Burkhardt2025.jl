const Burkhardt2025Groups = [
GCPair(raw"[CX4H3;!$([CX4H3]-[OH]);!$([CX4H3]-[NH2]);!$([CX4H3]-[CH3]);!$([CX4H3](-N=C=O))]", "-CH3"), # methyl
GCPair(raw"[CX4H2]", "-CH2-"), # methylene
GCPair(raw"[CX4H]", ">CH-"), # C-branch, single
GCPair(raw"[CX4H0]", ">C<"), # C-branch, double
GCPair(raw"[CX3H2]", "=CH2"), # alkenyl, terminal
GCPair(raw"[CX3H1]", "=CH-"), # alkenyl, non-terminal
GCPair(raw"[$([CX3H0,cX3H0]=*)]", ">C="), # alkenyl, non-terminal, branch
GCPair(raw"[$([CX2H0](=*)=*)]", "=C="), # alkenyl, =C=
GCPair(raw"[$([CX2H1]#*)]", "≡CH"), # alkynyl, terminal
GCPair(raw"[$([CX2H0]#*)]", "≡C-"), # alkynyl, non-terminal
GCPair(raw"[cX3H1]", "-CH-(arom)"), # aromatic CH
GCPair(raw"[cX3H0;!$(c(a)(a)a);!$([cX3H0]=[*])]", ">C-(arom)"), # aromatic >C-
GCPair(raw"[cX3H0;$(c(a)(a)a)]", ">C-(arom)(arom)"), # aromatic to aromatic >C-
GCPair(raw"[F]", "-F"), # fluoro
GCPair(raw"[Cl]", "-Cl"), # chloro
GCPair(raw"[Br]", "-Br"), # bromo
GCPair(raw"[I]", "-I"), # iodo
GCPair(raw"[OX2H]", "-OH"), # alcohol
GCPair(raw"[$([OX1]);$([OX1]=*);!$([OX1]=[#7X3,#7X3+]~[O-]);!$([OX1]=[SX4]=[OX1])]", "=O"), # keto/aldeh
GCPair(raw"[OX2H0]", "-O-"), # ether oxygen
GCPair(raw"[oX2H0]", "-O-(arom)"), # -O-(arom)
GCPair(raw"[NX3H2]", "-NH2"), # amine
GCPair(raw"[#7X2H1]", "=NH"), # imine =NH
GCPair(raw"[$([#7X1H0]#*)]", "≡N"), # cyano ≡N
GCPair(raw"[NX3H1]", "-NH-"), # amine -NH-
GCPair(raw"[$([NX2H0]=*)]", "-N="), # imine =N-
GCPair(raw"[nX2H0]", "-N-(arom)"), # -N-(arom)
GCPair(raw"[nX3H1]", "-NH-(arom)"), # -NH-(arom)
GCPair(raw"[$([#7X3H0]);!$([NX3H0](=[O])~[O-])]", ">N-"), # amine, branch or arom
GCPair(raw"[#7X3,#7X3+](=[O])~[O-]", "-NO2"), # nitro -NO2
GCPair(raw"[SX2H]", "-SH"), # thiol -SH
GCPair(raw"[$([SX1H0]=*)]", "=S"), # thioether
GCPair(raw"[SX2H0]", "-S-"), # thioether
GCPair(raw"[$([SX3](=*)(-*)-*)]", ">S="), # sulfoxide
GCPair(raw"[sX2]", "-S-(arom)"), # -S-(arom)
GCPair(raw"[$([SX4]-[C,c])](=[OX1])=[OX1]", ">SO2"), # sulfonyl >SO2
GCPair(raw"[PX3H2]", "-PH2"), # -PH2
GCPair(raw"[$([PX3H1](-*)-*)]", "-PH-"), # -PH-
GCPair(raw"[$([PX3H0](-*)(-*)-*)]", "-P<"), # -P<
GCPair(raw"[$([PX4H1](=*)(-*)-*)]", "=PH<"), # =PH<
GCPair(raw"[$([BX3H1](-*)-*)]", "-BH-"), # -BH-
GCPair(raw"[$([BX3H0](-*)(-*)-*)]", "-B<"), # -B<
GCPair(raw"[InX3H0]", ">In-"), # In-atom
GCPair(raw"[HgX2H0]", "-Hg-"), # Hg-atom
GCPair(raw"[SeX2H1]", "-SeH"), # (-Se-)-atom
GCPair(raw"[SeX2H0]", "-Se-"), # (-Se-)-atom
GCPair(raw"[SeX4H0]", ">Se<"), # (>Se<)-atom
GCPair(raw"[AlX3H0]", ">Al-"), # Al-atom
GCPair(raw"[AsX3H0]", ">As-"), # Al-atom
GCPair(raw"[GaX3H0]", ">Ga-"), # Ga-atom
GCPair(raw"[GeX4H0]", ">Ge<"), # Ge-atom
GCPair(raw"[SiX4H0]", ">Si<"), # Si-group
GCPair(raw"[SiX4H1]", ">SiH-"), # SiH-group
GCPair(raw"[SiX4H2]", "-SiH2-"), # SiH2-group
GCPair(raw"[SiX4H3]", "-SiH3"), # SiH3-group
GCPair(raw"[$([*](-[*;R])(-[*;R])(-[*;R]))]", "SO:multicyclic-site", group_order = 2), # SO:multicyclic-site
GCPair(raw"[CX4H2;R;!r4;!r3]", "SO:-CH2-(cyc5+)", group_order = 2), # SO:-CH2-(cyc5+)
GCPair(raw"[r3,r4;CX4H2]", "SO:-CH2-(cyc34)", group_order = 2), # SO:-CH2-(cyc34)
GCPair(raw"[$([CX4H3]-a)]", "SO:-CH3(@arom)", group_order = 2), # SO:SO:-CH3(@arom)
GCPair(raw"[$([OX2H]-a)]", "SO:-OH(@arom)", group_order = 2), # SO:-OH(@arom)
GCPair(raw"[F;$([F]-[CX4H2,CX3H1,CX2H0,CX4H3])]", "SO:-F(primary to C)", group_order = 2), # SO:-F(primary to C)
GCPair(raw"[CX3H1;R]", "SO:=CH-(cyc)", group_order = 2), # SO:=CH-(cyc)
GCPair(raw"[$([#6X3H0;R](=[OX1]))]=O", "SO:>C=O(keto,ring)", group_order = 2), # SO:>C=O(keto,ring)
GCPair(raw"[OX2H0;R]", "SO:-O-(cyc)", group_order = 2), # SO:-O-(cyc)
GCPair(raw"[Br;$([Br]-[CX4H2,CX3H1,CX2H0,CX4H3])]", "SO:-Br(primary to C)", group_order = 2), # SO:-Br(primary to C)
GCPair(raw"[Cl;$([Cl]-[CX4H2,CX3H1,CX2H0,CX4H3])]", "SO:-Cl(primary to C)", group_order = 2), # SO:-Cl(primary to C)
GCPair(raw"[$([F]-a)]", "SO:-F(@arom)", group_order = 2), # SO:-F(@arom)
GCPair(raw"[OX2H;$([OX2H]-[CX4H1,CX3H0;!$([CX3H0]=O)])]", "SO:-OH(secondary to C)", group_order = 2), # SO:-OH(secondary to C)
GCPair(raw"[$([C](-[F])(-[F])(-[F,Cl,Br,I])-C)]", "SO:-CF2X", group_order = 2), # SO:-CF2X
GCPair(raw"[Cl;$([Cl]-[Si])]", "SO:Cl-(to Si)", group_order = 2), # SO:Cl-(to Si)
GCPair(raw"[CX4H;R]", "SO:>CH-(cyc)", group_order = 2), # SO:>CH-(cyc)
GCPair(raw"[!R;CX4H2;$([CX4H2]-a)]", "SO:-CH2-(@arom)", group_order = 2), # SO:-CH2-(@arom)
GCPair(raw"[CX4H0;R]", "SO:>C<(cyc)", group_order = 2), # SO:>C<(cyc)
GCPair(raw"[$([*](-[Br,I])(-[Br,I])(-[Br,I]))]", "SO:[*]-[Br,I]3", group_order = 2), # SO:[*]-[Br,I]3
GCPair(raw"[NX3H1;R]", "SO:-NH-(cyc)", group_order = 2), # SO:-NH-(cyc)
GCPair(raw"[CX3H1;!$([CX3H1](~O)(~O))](=O)", "SO:-CH=O(aldeh)", group_order = 2), # SO:-CH=O(aldeh)
GCPair(raw"[OH1]-[OX2H0]", "SO:HO-O-(perhydroxide)", group_order = 2), # SO:HO-O-(perhydroxide)
]

export Burkhardt2025Groups

#=
vv = vcat(BG_SMARTS,BG_SMARTS_SEC)
for i in 1:77
    vi = vv[i]
    if i > 55
        println(io,"GCPair(raw\"",vi[3],"\", \"",vi[2],"\", group_order = 2), # ",vi[1])
    else
        println(io,"GCPair(raw\"",vi[3],"\", \"",vi[2],"\"), # ",vi[1])
    end
end

=#