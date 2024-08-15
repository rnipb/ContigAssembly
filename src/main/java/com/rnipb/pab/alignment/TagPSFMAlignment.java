package com.rnipb.pab.alignment;

/**
 * an object stores a PFSM alignment result
 */
public class TagPSFMAlignment {
    int tagId;
    int tagStart;
    int tagEnd;
    PSFM ref;
    int refId;
    int refStart; // inclusive
    int refEnd; // exclusive
    String alignedTagSeq;
    double score;
    double[] aaScore;

    public TagPSFMAlignment(int tagId, int tagStart, int tagEnd, PSFM ref, int refId, int refStart, int refEnd, String alignedTagSeq, double score, double[] aaScore) {
        this.tagId = tagId;
        this.tagStart = tagStart;
        this.tagEnd = tagEnd;
        this.ref = ref;
        this.refId = refId;
        this.refStart = refStart;
        this.refEnd = refEnd;
        this.alignedTagSeq = alignedTagSeq;
        this.score = score;
        this.aaScore = aaScore;
    }

    public int getTagId() {
        return tagId;
    }

    public int getTagStart() {
        return tagStart;
    }

    public int getTagEnd() {
        return tagEnd;
    }

    public PSFM getRef() {
        return ref;
    }

    public int getRefId() {
        return refId;
    }

    public int getRefStart() {
        return refStart;
    }

    public int getRefEnd() {
        return refEnd;
    }

    public String getAlignedTagSeq() {
        return alignedTagSeq;
    }

    public double getScore() {
        return score;
    }

    public double[] getAaScore() {
        return aaScore;
    }
}
