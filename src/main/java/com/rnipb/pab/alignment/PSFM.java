package com.rnipb.pab.alignment;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class PSFM {

    String name;
    List<Map<Character, Double>> aaRelativeEntropyMatrix; // relative entropy compared to Naturally probability of AA at each position

    PSFM (List<Map<Character, Double>> psmf, String name) {
        this.name = name;
        this.aaRelativeEntropyMatrix = new ArrayList<>(psmf);
    }

    public double aaRelativeEntropy(int pos, char aa) {
        if (aa == 'I')
            aa = 'L';
        return aaRelativeEntropyMatrix.get(pos).get(aa);
    }

    /**
     * @return length of the PSFM matrix
     */
    public int length() {
        return aaRelativeEntropyMatrix.size();
    }

    public String getName() {
        return name;
    }

}
