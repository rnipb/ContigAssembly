package com.rnipb.pab;

import com.rnipb.pab.alignment.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ContigAssembly {
    double MASS_TOLARANCE = 0.02;
    int PSM_THRESHOLD = 75;

    List<PSM> denovoPSMs;
    PSFMLoader psfmLoader;

    public static void main(String[] args) throws IOException {
        String denovoRoot = args[0]; // path of a folder that has all de novo result files
        String psfmPath = args[1]; // a parameter file for PSFM, used for sequence alignment
        String denovoFileSuffix = args[2]; // suffix of de novo files, all files with the suffix in denovoRoot will be used
        PSFMLoader psfmLoader = new PSFMLoader(psfmPath);
        ContigAssembly inst = new ContigAssembly(denovoRoot, psfmLoader, denovoFileSuffix);
        List<String[]> sequences = inst.overlap(); // find overlap between De Novo peptides to build longer sequences, output results in contigs.fasta
        writeFasta(sequences, new File(denovoRoot, "contigs.fasta")); //output contigs to contigs.fasta
        List<String[]> filtered = inst.filterSequences(sequences);// filter out unlikely sequences based on similarity to germline
        writeFasta(filtered, new File(denovoRoot, "contigs_fil.fasta")); //output filtered results to contig_fil.fasta
    }

    public ContigAssembly(String denovoRoot, PSFMLoader psfmLoader, String denovoFileSuffix) throws IOException {
        denovoPSMs = new ArrayList<>();
        this.psfmLoader = psfmLoader;
        int csvId = 0;
        // read De Novo result files
        for (String s : Objects.requireNonNull(new File(denovoRoot).list())) {
            if (s.endsWith(denovoFileSuffix)) {
                csvId++;
                List<PSM> psms = readCSV(new File(denovoRoot, s), csvId);
                denovoPSMs.addAll(psms);
                System.out.printf("%d De novo PSMs loaded from %s%n", psms.size(), s);
            }
        }
        // exit if no file is found
        if (denovoPSMs.isEmpty()) {
            System.out.println("No PSM was found.");
            System.exit(1);
        } else {
            System.out.printf("%d De novo PSMs loaded in total%n", denovoPSMs.size());
        }
    }

    private List<PSM> readCSV(File file, int csvId) throws IOException {
        List<PSM> denovoPSMs = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(file));
        for (String line = br.readLine(); line != null; line = br.readLine()) {
            line = line.trim().replaceAll("\\s", "");
            if (line.trim().startsWith("#")) {
                continue;
            }
            String[] split = line.split(",");
            if (split.length < 11) {
                continue;
            }
            denovoPSMs.add(new PSM(csvId, Integer.parseInt(split[1]), Double.parseDouble(split[2]), Double.parseDouble(split[3]), Integer.parseInt(split[4]),
                    removeCharactersInBrackets(split[9]), Double.parseDouble(split[8]), Arrays.stream(split[10].split("-")).mapToInt(Integer::parseInt).toArray()));
        }
        return denovoPSMs;
    }

    public static String removeCharactersInBrackets(String input) {
        Pattern pattern = Pattern.compile("\\([^()]*\\)");
        Matcher matcher = pattern.matcher(input);
        StringBuffer output = new StringBuffer();

        while (matcher.find()) {
            matcher.appendReplacement(output, "");
        }
        matcher.appendTail(output);

        return output.toString();
    }

    private List<String[]> overlap() {
        // only PSMs whose score above the threshold will be used
        List<PSM> goodPSMs = denovoPSMs.stream()
                .filter(psm -> psm.score > PSM_THRESHOLD)
                .collect(Collectors.toList());
        // for each unique peptide sequence, find the best De Novo AA score for each AA
        Map<String, int[]> pepSeqsAAScore = new ConcurrentHashMap<>();
        for (PSM psm : goodPSMs) {
            String seq = psm.pep;
            if (pepSeqsAAScore.containsKey(seq)) {
                int[] scores = pepSeqsAAScore.get(seq);
                int[] max = IntStream.range(0, scores.length).boxed().mapToInt(i -> Math.max(scores[i], psm.aaScore[i])).toArray();
                pepSeqsAAScore.put(seq, max);
            } else {
                pepSeqsAAScore.put(seq, psm.aaScore);
            }
        }
        // do alignment for each De Novo peptide sequence, sequences can't be aligned (alignment score < 0) will be discarded
        List<PSFM> psfms = psfmLoader.getPsfms();
        List<TagPSFMAlignment> aligns = pepSeqsAAScore.keySet().stream()
                .map(s -> PSFMAligner.searchSeq(s, psfms, 0, 0.0))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        // save contigs built from overlapping of De Novo peptides in contigsAll
        List<List<Map<String, List<String>>>> contigsAll = new ArrayList<>();
        for (int icdr = 0; icdr <= 2; icdr++) {
            int finalIcdr = icdr;
            Set<String> pepSeqs = aligns.stream()
                    .filter(ta -> ta.getRefId() == finalIcdr) // use peptide from one type of chain (heavy, kappa, lambda)
                    .filter(ta -> ta.getRefStart() < psfmLoader.getkPos()[finalIcdr]) // use only v-region peptides
                    .map(ta -> ta.getAlignedTagSeq().replace(".", ""))
                    .collect(Collectors.toSet());
            Map<String, int[]> contigAAScore = new ConcurrentHashMap<>(pepSeqsAAScore);
            // build longer sequences (contigs) from pepSeqs
            List<String> contigs = extendContigs(pepSeqs, contigAAScore, 3, 3);
            // group different contigs that have same cdr regions
            List<Map<String, List<String>>> groupedContigs = groupingContigs(contigs, icdr);
            contigsAll.add(groupedContigs);
        }
        // naming the contigs, contigs having the same CDR will be given a same group name
        List<String[]> sequences = namingSequences(contigsAll);
        return sequences;
    }

    private List<String[]> namingSequences(List<List<Map<String, List<String>>>> contigsAll) {
        List<String[]> sequences = new ArrayList<>();
        for (int icdr = 0; icdr < contigsAll.size(); icdr++) {
            String chainType = icdr == 0 ? "hcdr" : (icdr == 1 ? "kcdr" : "lcdr");
            for (int i = 0; i < 4; i++) {
                String cdrType = i == 0 ? "1" : (i == 1 ? "2" : (i == 2 ? "12" : "3"));
                Map<String, List<String>> cdrContigMap = contigsAll.get(icdr).get(i);
                List<Map.Entry<String, List<String>>> entries = new ArrayList<>(cdrContigMap.entrySet());
                for (int j = 0; j < entries.size(); j++) {
                    Map.Entry<String, List<String>> entry = entries.get(j);
                    String cdrSeq = entry.getKey().replace(".", "");
                    List<String> contigs = entry.getValue();
                    for (int k = 0; k < contigs.size(); k++) {
                        String contig = contigs.get(k);
                        String name = ">" + chainType + cdrType + "-c" + (j+1) + "." + (k+1) + " " + cdrSeq;
                        sequences.add(new String[]{name, contig});
                    }
                }
            }
        }
        return sequences;
    }

    private List<Map<String, List<String>>> groupingContigs(List<String> contigs, int icdr) {
        List<PSFM> psfms = psfmLoader.getPsfms();
        List<TagPSFMAlignment> aligns = contigs.stream()
                .map(s -> PSFMAligner.searchSeq(icdr, s, psfms.get(icdr), 0, 0.0))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        Map<String, List<String>> cdr1 = new ConcurrentHashMap<>();
        Map<String, List<String>> cdr2 = new ConcurrentHashMap<>();
        Map<String, List<String>> cdr12 = new ConcurrentHashMap<>();
        Map<String, List<String>> cdr3 = new ConcurrentHashMap<>();
        int[][] cdrRegions = psfmLoader.getCdrRegions()[icdr];
        cdrRegions[0] = new int[]{cdrRegions[0][0] - 2, cdrRegions[0][1] + 2}; // use 2 more AA on each side of cdr1 region
        cdrRegions[1] = new int[]{cdrRegions[1][0] - 2, cdrRegions[1][1] + 2}; // use 2 more AA on each side of cdr2 region
        aligns.stream().parallel().forEach(ta -> {
            int refStart = ta.getRefStart();
            int refEnd = ta.getRefEnd();
            String contigSeq = ta.getAlignedTagSeq().replace(".", "");
            if (refStart <= cdrRegions[0][0] && refEnd >= cdrRegions[0][1]) {
                String cdr = ta.getAlignedTagSeq().substring(cdrRegions[0][0] - refStart, cdrRegions[0][1] - refStart);
                List<String> list = cdr1.computeIfAbsent(cdr, key -> new ArrayList<>());
                synchronized (list) {
                    list.add(contigSeq);
                }
            }
            if (refStart <= cdrRegions[1][0] && refEnd >= cdrRegions[1][1]) {
                String cdr = ta.getAlignedTagSeq().substring(cdrRegions[1][0] - refStart, cdrRegions[1][1] - refStart);
                List<String> list = cdr2.computeIfAbsent(cdr, key -> new ArrayList<>());
                synchronized (list) {
                    list.add(contigSeq);
                }
            }
            if (refStart <= cdrRegions[0][0] && refEnd >= cdrRegions[1][1]) {
                String cdr = ta.getAlignedTagSeq().substring(cdrRegions[0][0] - refStart, cdrRegions[1][1] - refStart);
                List<String> list = cdr12.computeIfAbsent(cdr, key -> new ArrayList<>());
                synchronized (list) {
                    list.add(contigSeq);
                }
            }
            if (refStart <= cdrRegions[2][0] && refEnd >= cdrRegions[2][1] + 4) {
                String cdr = ta.getAlignedTagSeq().substring(cdrRegions[2][0] - refStart, cdrRegions[2][1] + 4 - refStart);
                List<String> list = cdr3.computeIfAbsent(cdr, key -> new ArrayList<>());
                synchronized (list) {
                    list.add(contigSeq);
                }
            }
        });
        List<Map<String, List<String>>> allContigs = Arrays.asList(cdr1, cdr2, cdr12, cdr3);
        return allContigs;
    }

    // extend contigs by iteratively attach peptides to current contigs, initial contigs are pepSeqs
    // minOverlap: minimum number of mass block overlap
    // round: a contig can be made up from at most round peptides
    public List<String> extendContigs(Set<String> pepSeqs, Map<String, int[]> pepSeqsAAScore, int minOverlap, int round) {
        List<Set<String>> contigLists = new ArrayList<>();
        contigLists.add(pepSeqs);
        for (int i = 1; i < round; i++) {
            contigLists.add(expand(contigLists.get(i - 1), pepSeqs, pepSeqsAAScore, minOverlap, MASS_TOLARANCE));
        }
        return contigLists.stream().flatMap(Collection::stream).distinct().collect(Collectors.toList());
    }

    // for every contig in baseContigs and every peptide in readSeqs, check if a peptide can be appended to a contig
    // if there are at least minMassBlocksOverlap between a contig and a peptide, extend the contig by the peptide to form a new longer contig
    protected Set<String> expand(Collection<String> baseContigs, Collection<String> readSeqs, Map<String, int[]> pepSeqAAScores, int minMassBlocksOverlap, double massTolerance) {
        List<String> contigs = new ArrayList<>(baseContigs);
        List<String> reads = new ArrayList<>(readSeqs);
        Map<String, double[]> prefixMassMap = reads.parallelStream().collect(Collectors.toConcurrentMap(seq -> seq, ContigAssembly::getPrefixMass));
        Set<String> expandedContigs = ConcurrentHashMap.newKeySet();
        IntStream.range(0, contigs.size()).parallel().forEach(i -> {
            String seq1 = contigs.get(i);
            double[] pep1SuffixMass = getSuffixMass(seq1);
            int[] pep1AAScores = pepSeqAAScores.get(seq1);
            for (String seq2 : reads) {
                double[] pep2PrefixMass = prefixMassMap.get(seq2);
                int[] pep2AAScores = pepSeqAAScores.get(seq2);
                // conditions that seq2 will not be used to extend
                if (seq1.contains(seq2) || seq2.contains(seq1))
                    continue;
                if (seq1.length() < minMassBlocksOverlap || seq2.length() < minMassBlocksOverlap)
                    continue;
                // find max mass overlap
                int p1I = 0;
                int p2I = seq2.length() - 1;
                while (p1I < seq1.length() && p2I >= 0) {
                    if (pep1SuffixMass[p1I] - pep2PrefixMass[p2I] > massTolerance) {
                        p1I++;
                    } else if (pep2PrefixMass[p2I] - pep1SuffixMass[p1I] > massTolerance) {
                        p2I--;
                    } else {// when seq1.substring(p1I) and seq2.substring(0, p2I+1) have equal mass
                        // check number of mass matched blocks and daerau distance of overlapping subseq
                        // if condition satisfied, create new contigs, otherwise move the indicator positions
                        String overlapS1 = seq1.substring(p1I);
                        String overlapS2 = seq2.substring(0, p2I + 1);
                        if (massMatchBlocks(overlapS1, overlapS2, massTolerance) < minMassBlocksOverlap) {// || damerauDistance(overlapS1, overlapS2) > maxOverlapDist
                            p1I++;
                            p2I--;
                            continue;
                        }
                        // concatenate seq1 and seq2
                        int[] overlapAAS1 = Arrays.copyOfRange(pep1AAScores, p1I, pep1AAScores.length);
                        int[] overlapAAS2 = Arrays.copyOfRange(pep2AAScores, 0, p2I + 1);
                        Pair<String, int[]> overlap = mergeMassMatchBlocks(overlapS1, overlapAAS1, overlapS2, overlapAAS2, massTolerance);
                        // add new contigs to contig set
                        String overlapSeq = overlap.first();
                        int[] overlapAA = overlap.second();
                        String newSeq = seq1.substring(0, p1I) + overlapSeq + seq2.substring(p2I + 1);
                        int[] newAAScore = concatenate(concatenate(Arrays.copyOfRange(pep1AAScores, 0, p1I), overlapAA), Arrays.copyOfRange(pep2AAScores, p2I + 1, pep2AAScores.length));
                        expandedContigs.add(newSeq);
                        if (pepSeqAAScores.containsKey(newSeq)) {
                            int[] scores = pepSeqAAScores.get(newSeq);
                            int[] max = IntStream.range(0, scores.length).boxed().mapToInt(ii -> Math.max(scores[ii], newAAScore[ii])).toArray();
                            pepSeqAAScores.put(newSeq, max);
                        } else {
                            pepSeqAAScores.put(newSeq, newAAScore);
                        }
                        break;
                    }
                }
            }
        });
        return expandedContigs;
    }

    private Pair<String, int[]> mergeMassMatchBlocks(String s1, int[] aa1, String s2, int[] aa2, double massTolerance) {
        StringBuilder s = new StringBuilder();
        List<Integer> aaScore = new ArrayList<>();
        double[] s1PrefixMass = getPrefixMass(s1);
        double[] s2PrefixMass = getPrefixMass(s2);
        int s1I = 0, s2I = 0;
        int aaSum1 = 0, aaSum2 = 0;
        int lastS1I = 0, lastS2I = 0;
        while (s1I < s1.length() && s2I < s2.length()) {
            if (s1PrefixMass[s1I] - s2PrefixMass[s2I] > massTolerance) {
                aaSum2 += aa2[s2I];
                s2I++;
            } else if (s2PrefixMass[s2I] - s1PrefixMass[s1I] > massTolerance) {
                aaSum1 += aa1[s1I];
                s1I++;
            } else {
                // matched block
                aaSum1 += aa1[s1I++];
                aaSum2 += aa2[s2I++];
                s.append(aaSum1 >= aaSum2 ? s1.substring(lastS1I, s1I) : s2.substring(lastS2I, s2I));
                aaScore.addAll(aaSum1 >= aaSum2 ? IntStream.range(lastS1I, s1I).boxed().map(i -> aa1[i]).collect(Collectors.toList())
                        : IntStream.range(lastS2I, s2I).boxed().map(i -> aa2[i]).collect(Collectors.toList()));
                lastS1I = s1I;
                lastS2I = s2I;
                aaSum1 = 0;
                aaSum2 = 0;
            }
        }
        return new Pair<>(s.toString(), aaScore.stream().mapToInt(i -> i).toArray());
    }

    // filter out sequences that is not aligned with ref
    private List<String[]> filterSequences(List<String[]> sequences) {
        List<PSFM> psfms = psfmLoader.getPsfms();
        Map<String, TagPSFMAlignment> seqAlignments = new ConcurrentHashMap<>();
        sequences.stream().parallel().forEach(seq -> {
            int chainType = seq[0].startsWith(">h") ? 0 : (seq[0].startsWith(">k") ? 1 : 2);
            TagPSFMAlignment ta = PSFMAligner.searchSeq(chainType, seq[1], psfms.get(chainType), 0, 20.0);
            if (ta != null) {
                seqAlignments.put(seq[0], ta);
            }
        });
        Set<String> sequenceSet = new ConcurrentSkipListSet<>();
        List<String[]> filteredSeqs = seqAlignments.entrySet().stream().parallel()
                .map(entry -> {
                    String seqName = entry.getKey();
                    TagPSFMAlignment ta = entry.getValue();
                    String alignedTagSeq = ta.getAlignedTagSeq();
                    String seq = alignedTagSeq.replace(".", "");
                    String[] split = seqName.split(" ");
                    if (split.length != 2) {
                        return null;
                    }
                    String cdrSeq = split[1];
                    int[] cdrRange = new int[]{seq.indexOf(cdrSeq), seq.indexOf(cdrSeq) + cdrSeq.length()};
                    int[] cdrRangeExt = new int[]{seq.indexOf(cdrSeq), seq.indexOf(cdrSeq) + cdrSeq.length()};
                    if (seqName.contains("cdr3")) {
                        cdrRange[1] -= 4;
                    } else {
                        cdrRange[0] += 2;
                        cdrRange[1] -= 2;
                    }
                    double[] aaScore = ta.getAaScore();
                    double[] aaScoreUnalign = new double[seq.length()];
                    for (int i = 0, j = 0; i < aaScoreUnalign.length; i++, j++) {
                        while (alignedTagSeq.charAt(j) == '.') {
                            j++;
                        }
                        aaScoreUnalign[i] = aaScore[j];
                    }
                    // check the interval before cdr and after cdr, if any one of this interval is obviously different
                    // from germline (all AA have small negative scores), remove this sequence
                    double[] aaBefore = Arrays.copyOfRange(aaScoreUnalign, 0, cdrRange[0]);
                    double[] aaAfter = Arrays.copyOfRange(aaScoreUnalign, cdrRange[1], aaScoreUnalign.length);
                    if (Arrays.stream(aaBefore).filter(sc -> sc >= 3).count() < 2 || Arrays.stream(aaAfter).filter(sc -> sc >= 3).count() < 2) {
                        return null;
                    }
                    int leftTrim = 0;
                    while (aaBefore[leftTrim] <= -3) {
                        leftTrim++;
                    }
                    int rightTrim = aaAfter.length - 1;
                    while (aaAfter[rightTrim] <= -3) {
                        rightTrim--;
                    }
                    seq = seq.substring(Math.min(leftTrim, cdrRangeExt[0]), Math.max(cdrRangeExt[1], cdrRange[1] + rightTrim + 1));
                    if (sequenceSet.contains(seq)) {
                        return null;
                    }
                    sequenceSet.add(seq);
                    return new String[]{seqName, seq};
                })
                .filter(Objects::nonNull)
                .collect(Collectors.toList());

        return filteredSeqs;
    }

    // number of matched sub blocks of two seqs
    protected static int massMatchBlocks(String s1, String s2, double massTolerance) {
        double[] s1PrefixMass = getPrefixMass(s1);
        double[] s2PrefixMass = getPrefixMass(s2);
        int s1I = 0, s2I = 0;
        int matchedBlock = 0;
        while (s1I < s1.length() && s2I < s2.length()) {
            if (s1PrefixMass[s1I] - s2PrefixMass[s2I] > massTolerance) {
                s2I++;
            } else if (s2PrefixMass[s2I] - s1PrefixMass[s1I] > massTolerance) {
                s1I++;
            } else {
                matchedBlock++;
                s1I++;
                s2I++;
            }
        }
        return matchedBlock;
    }

    public static double[] getPrefixMass(String seq) {
        double[] prefixMass = new double[seq.length()];
        prefixMass[0] = getResidueMass(seq.charAt(0));
        for (int i = 1; i < seq.length(); i++) {
            prefixMass[i] = prefixMass[i - 1] + getResidueMass(seq.charAt(i));
        }
        return prefixMass;
    }

    public static double[] getSuffixMass(String seq) {
        double[] suffixMass = new double[seq.length()];
        suffixMass[seq.length() - 1] = getResidueMass(seq.charAt(seq.length() - 1));
        for (int i = seq.length() - 2; i >= 0; i--) {
            suffixMass[i] = suffixMass[i + 1] + getResidueMass(seq.charAt(i));
        }
        return suffixMass;
    }

    private static double getResidueMass(char aa) {
        final double camOffset = 57.021464;
        return aa == 'C' ? SimpleResidue.charCodeToResidue(aa).getMass() + camOffset : SimpleResidue.charCodeToResidue(aa).getMass();
    }

    public static int[] concatenate(int[] a, int[] b) {
        int[] result = new int[a.length + b.length];
        System.arraycopy(a, 0, result, 0, a.length);
        System.arraycopy(b, 0, result, a.length, b.length);
        return result;
    }

    private static void writeFasta(List<String[]> sequences, File file){
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file))) {
            for (String[] sequence : sequences) {
                String name = sequence[0].startsWith(">") ? sequence[0] : ">" + sequence[0];
                bw.write(name + "\n");
                bw.write(sequence[1] + "\n");
                bw.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    // a class that stores peptide spectrum matching information
    public static class PSM {
        int msNum;
        int scanNum;
        double rt;
        double mz;
        int z;
        String pep;
        double score;
        int[] aaScore;

        public PSM(int msNum, int scanNum, double rt, double mz, int z, String pep, double score, int[] aaScore) {
            this.msNum = msNum;
            this.scanNum = scanNum;
            this.rt = rt;
            this.mz = mz;
            this.z = z;
            this.pep = pep;
            this.score = score;
            this.aaScore = aaScore;
        }
    }
}
