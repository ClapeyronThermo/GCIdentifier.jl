gcPCSAFTGroups = [
    GCPair(raw"[CX4H3]", "CH3"),
    GCPair(raw"[!R;CX4H2]", "CH2"),
    GCPair(raw"[!R;CX4H]", "CH"),
    GCPair(raw"[!R;CX4H0]", "C"),
    GCPair(raw"[CX3H2]", "CH2="),
    GCPair(raw"[!R;CX3H1;!$([CX3H1](=O))]", "CH="),
    GCPair(raw"[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]", "=C<"),
    GCPair(raw"[CX2;H1]#[CX2;H0]", "C#CH"),
    GCPair(raw"[CH2;R1;$(C1CCCC1)]", "cCH2_pen"),
    GCPair(raw"[CH1;R1;$(C1CCCC1)]", "cCH_pen"),
    GCPair(raw"[CH2;R1;$(C1CCCCC1)]", "cCH2_hex"),
    GCPair(raw"[CH1;R1;$(C1CCCCC1)]", "cCH_hex"),
    GCPair(raw"[cX3;H1]", "aCH"),
    GCPair(raw"[cX3;H0]", "aCH"),
    GCPair(raw"[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]", "OH"),
    GCPair(raw"[NX3H2]", "NH2"),
    GCPair(raw"[O;H2]", "H2O")
    ]

export gcPCSAFTGroups
