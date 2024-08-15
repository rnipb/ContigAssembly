package com.rnipb.pab;

public enum SimpleResidue {
    A("Alanine", "Ala", "C3H5NO", 71.037113787),
    R("Arginine", "Arg", "C6H12N4O", 156.101111026),
    N("Asparagine", "Asn", "C4H6N2O2", 114.042927446),
    D("Aspartic Acid", "Asp", "C4H5NO3", 115.026943031),
    C("Cysteine", "Cys", "C3H5NOS", 103.00918447699999),
    E("Glutamic Acid", "Glu", "C5H7NO3", 129.042593095),
    Q("Glutamine", "Gln", "C5H8N2O2", 128.05857751),
    G("Glycine", "Gly", "C2H3NO", 57.021463723),
    H("Histidine", "His", "C6H7N3O", 137.058911861),
    I("Isoleucine", "Ile", "C6H11NO", 113.08406397899999),
    L("Leucine", "Leu", "C6H11NO", 113.08406397899999),
    K("Lysine", "Lys", "C6H12N2O", 128.094963016),
    M("Methionine", "Met", "C5H9NOS", 131.040484605),
    F("Phenyalanine", "Phe", "C9H9NO", 147.068413915),
    O("Pyrrolysine", "Pyl", "C12H19N3O2", 237.147726867),
    P("Proline", "Pro", "C5H7NO", 97.052763851),
    S("Serine", "Ser", "C3H5NO2", 87.03202840899999),
    T("Threonine", "Thr", "C4H7NO2", 101.04767847299999),
    W("Tryptophan", "Trp", "C11H10N2O", 186.079312952),
    Y("Tyrosine", "Tyr", "C9H9NO2", 163.063328537),
    U("Selenocysteine", "Sec", "C3H5NOSe", 150.953635587),
    V("Valine", "Val", "C5H9NO", 99.068413915);

    private String fullName;
    private String code3;
    private char charCode;
    private byte byteCode;
    private static final SimpleResidue[] byteCodeToResidue = new SimpleResidue[26];
    private double mass;

    private SimpleResidue(String fullName, String code3, String composition, double mass) {
        this.code3 = code3;
        this.fullName = fullName;
        this.charCode = this.name().charAt(0);
        this.byteCode = (byte)(this.charCode - 65);
        this.mass = mass;
    }

    public double getMass() {
        return this.mass;
    }

    public static SimpleResidue charCodeToResidue(char code) {
        int k = code - 65;
        return k >= 0 && k < byteCodeToResidue.length ? byteCodeToResidue[k] : null;
    }

    static {
        SimpleResidue[] var0 = values();
        int var1 = var0.length;

        for(int var2 = 0; var2 < var1; ++var2) {
            SimpleResidue residue = var0[var2];
            byteCodeToResidue[residue.name().charAt(0) - 65] = residue;
        }
    }
}
