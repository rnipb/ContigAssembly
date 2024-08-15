package com.rapidnovor.pab.alignment;

import java.util.List;

public class PSFMAligner {
    // alignment score: sum( log( p(x) / q(x)) )
    // q(x) = 1.0 if x=='.' else q(x) = AA natural frequency
    // p(x) is the probability of AA (or '.') in alignedTagSeq at the matched position in PSFM
    // the score equals to RANDOM_AA_ENTROPY*tagSeq.length() - sum(log(p(x))),
    // physical meaning: difference of encoding length between using random AA prob and using AA prob in aligned position on PSFM to encode the tag sequence

    private static final double NEGATIVE_INFINITY = Double.NEGATIVE_INFINITY;

    /** when refId is not provided, try all psfm in the list, the refId is the index of the best psfm
     * @param proteinSeq AA string sequence for alignment
     * @param psfms a list of psfm object
     * @param tagId fraction Id of the denovo tag
     * @param minAlignmentScore minimum sum of relative entropy of tag
     * @return alignment result object
     */
    public static TagPSFMAlignment searchSeq(String proteinSeq, List<PSFM> psfms, int tagId, double minAlignmentScore) {
        TagPSFMAlignment bestTa = null;
        for (int refId = 0; refId < psfms.size(); refId++) {
            if (psfms.get(refId) != null) {
                TagPSFMAlignment ta = searchSeq(refId, proteinSeq, psfms.get(refId), tagId, minAlignmentScore);
                if (bestTa == null || (ta != null && ta.getScore() > bestTa.getScore())) {
                    bestTa = ta;
                }
            }
        }
        return bestTa;
    }

    /**
     * @param refId 0, 1, 2 corresponding to IGH, IGK, IGL respectively, (index of PSFM object in List<PSFM> psfms)
     * @param proteinSeq AA string sequence for alignment
     * @param psfm a psfm object
     * @param tagId fraction Id of the denovo tag
     * @param minAlignmentScore minimum sum of relative entropy of tag
     * @return alignment result object
     */
    public static TagPSFMAlignment searchSeq(int refId, String proteinSeq, PSFM psfm, int tagId, double minAlignmentScore) {
        if (proteinSeq == null)
            return null;
        TagPSFMAlignment align = align(refId, proteinSeq, psfm, tagId, minAlignmentScore);
        return align;
    }

    /**
     * @param refId 0, 1, 2 corresponding to IGH, IGK, IGL respectively, (index of PSFM object in List<PSFM> psfms)
     * @param tagSeq AA string sequence for alignment
     * @param psfm a psfm object
     * @param tagId  fraction Id of the denovo tag
     * @param minAlignmentScore minimum sum of relative entropy of tag
     * @return a TagPSFMAlignment object if the best alignment has score above the threshold
     */
    public static TagPSFMAlignment align(int refId, String tagSeq, PSFM psfm, int tagId, double minAlignmentScore) {
        tagSeq = tagSeq.replace(".", "");
        int tagLen = tagSeq.length();
        int refLen = psfm.length();
        double[][] dp = new double[refLen+1][tagLen+1];
        int[][] backtrace = new int[refLen+1][tagLen+1];
        // prefix of tag must be aligned
        for (int j = 1; j <=tagLen; j++) {
            dp[0][j] = NEGATIVE_INFINITY;
        }

        for (int i = 1; i <= refLen; i++) {
            for (int j = 1; j <= tagLen; j++) {
                dp[i][j] = dp[i-1][j-1] + psfm.aaRelativeEntropy(i-1, tagSeq.charAt(j-1));
                backtrace[i][j] = 1;
                // try insertion in tag seq, the penalty is based on prob of '.'
                double insertScore = dp[i-1][j] + psfm.aaRelativeEntropy(i-1, '.');
                if (insertScore > dp[i][j]) {
                    dp[i][j] = insertScore;
                    backtrace[i][j] = 2;
                }
            }
        }
        // find the best alignment score
        double bestScore = NEGATIVE_INFINITY;
        int refEnd = refLen;
        int tagEnd = tagLen;

        // check last column of dp[][],
        // suffix of tag is aligned: 1. tag into ref (tInr), 3. suffix of tag vs prefix of ref (stpr)
        for (int i = 1; i <= refLen; i++) {
            if (dp[i][tagLen] > bestScore) {
                refEnd = i;
                bestScore = dp[i][tagLen];
            }
        }

        if (bestScore < minAlignmentScore)
            return null;
        // find the best alignment path
        StringBuilder alignedTagBuilder = new StringBuilder();
        int i = refEnd;
        int j = tagEnd;
        while (j > 0 && i > 0) {
            if (backtrace[i][j] == 1) {
                alignedTagBuilder.insert(0, tagSeq.charAt(j - 1));
                j--;
            }else if (backtrace[i][j] == 2) {
                alignedTagBuilder.insert(0, '.');
            }
            i--;
        }
        int refStart = i;
        int tagStart = j;
        String alignedTagSeq = alignedTagBuilder.toString();
        double[] aaScore = new double[alignedTagSeq.length()];
        for (int k = 0; k < alignedTagSeq.length(); k++) {
            aaScore[k] = psfm.aaRelativeEntropy(k+refStart, alignedTagSeq.charAt(k));
        }
        TagPSFMAlignment tagPSFMAlignment = new TagPSFMAlignment(tagId, tagStart, tagEnd, psfm, refId, refStart, refEnd, alignedTagSeq, bestScore, aaScore);
        return tagPSFMAlignment;
    }
}
