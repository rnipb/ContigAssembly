package com.rapidnovor.pab;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class SortContigs {
    Map<String, String[]> proPepLinks;
    Map<String, List<String[]>> psms;
    Map<String, String> sequences;
    List<String[]> proteins;
    String root;

    // sort the sequences in contigs_fil.fasta based on the database search result from novor.cloud
    public static void main(String[] args) throws IOException {
        // root is a folder path contains contigs_fil.fasta from ContigAssembly.java,
        // and proteins.csv, proPepLinks.txt, psms.csv from novor.cloud searched with contigs_fil.fasta
        String root = args[0];
        String proteinsPath = new File(root, "proteins.csv").getAbsolutePath();
        String proPepLinksPath = new File(root, "proPepLinks.txt").getAbsolutePath();
        String psmsPath = new File(root, "psms.csv").getAbsolutePath();
        String fasta = new File(root, "contigs_fil.fasta").getAbsolutePath();

        List<String[]> proteins = readFile(proteinsPath, true);
        Map<String, String[]> proPepLinks = readProPepLinks(proPepLinksPath);
        List<String[]> psms = readFile(psmsPath, true);
        List<String[]> sequences = loadFasta(new File(fasta));
        SortContigs instance = new SortContigs(proteins, proPepLinks, psms, sequences, root);
        instance.analyse();
    }

    public SortContigs(List<String[]> proteins, Map<String, String[]> proPepLinks, List<String[]> psms, List<String[]> sequences, String root) {
        this.root = root;
        this.sequences = sequences.stream().collect(Collectors.toMap(s -> s[0].split("\\s")[0], s -> s[1]));
        this.proPepLinks = proPepLinks;
        this.proteins = proteins;
        this.psms = new HashMap<>();
        for (String[] psm : psms) {
            String psmId = psm[0].split("-")[0];
            this.psms.computeIfAbsent(psmId, key -> new ArrayList<>()).add(psm);
        }
    }

    private void writeFasta(Map<String, List<List<ProteinCoverage>>> sortedGroupByType) {
        File file = new File(root, "sorted.fasta");
        // can be modified
        int groupNum = 35; //the number of groups to be output
        int nInGroup = 5; // number of contigs in each group to be output
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file))) {
            for (String type : sortedGroupByType.keySet()) {
                if (type.contains("12"))
                    continue;
                List<List<ProteinCoverage>> groups = sortedGroupByType.get(type);
                for (int i = 0; i < Math.min(groupNum, groups.size()); i++) {
                    List<ProteinCoverage> groupI = groups.get(i);
                    for (int j = 0; j < Math.min(nInGroup, groupI.size()); j++) {
                        ProteinCoverage pc = groupI.get(j);
                        bw.write(String.format(">%s %s\n", pc.pName, pc.cdrSeq));
                        bw.write(pc.proteinSeq + "\n");
                    }
                    bw.write("\n");
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    // read fasta file
    private static List<String[]> loadFasta(File file) {
        List<String[]> sequences = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line;
            String currentHeader = null;
            StringBuilder sequenceBuilder = new StringBuilder();
            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (currentHeader != null) {
                        sequences.add(new String[]{currentHeader, sequenceBuilder.toString()});
                        sequenceBuilder = new StringBuilder(); // Reset for the next sequence
                    }
                    currentHeader = line.trim();
                } else {
                    sequenceBuilder.append(line.trim()); // Append sequence lines, trimming whitespace
                }
            }
            if (currentHeader != null) {
                sequences.add(new String[]{currentHeader, sequenceBuilder.toString()}); // Add last sequence
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return sequences;
    }

    private void analyse() {
        // create a ProteinCoverage object for each contig sequence, the peptides for each contig are stored in ProteinCoverage.proteinPeps
        List<ProteinCoverage> pcs = proteins.stream().parallel()
                .map(protein -> {
                    String gId = protein[0];
                    String[] split = protein[7].replace("\"", "").split(" ");
                    if (split.length !=  2) {
                        return null;
                    }
                    String pName = split[0];
                    String cdrSeq = split[1];
                    if (!this.sequences.containsKey(">" + pName)) {
                        return null;
                    }
                    String proteinSeq = this.sequences.get(">" + pName);
                    List<List<String[]>> proteinPeps = Arrays.stream(this.proPepLinks.get(gId)).map(s -> this.psms.get(s)).collect(Collectors.toList());
                    return new ProteinCoverage(gId, pName, cdrSeq, proteinSeq, proteinPeps);
                })
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        // group ProteinCoverage objects based on their CDR regions
        Map<String, List<ProteinCoverage>> groupedPcs = new HashMap<>();
        for (ProteinCoverage pc : pcs) {
            groupedPcs.computeIfAbsent(pc.pName.split("\\.")[0], key -> new ArrayList<>()).add(pc);
        }
        // sort the groups and contigs in each group
        List<List<ProteinCoverage>> sortedGroups = groupedPcs.values().stream()
                .map(l -> l.stream()
                        .sorted(Comparator.comparingInt(pc -> -pc.minCdrDepth))
                        .sorted(Comparator.comparingInt(pc -> -pc.longPeptide))
                        .sorted(Comparator.comparingDouble(pc -> -pc.score))
                        .collect(Collectors.toList()))
                .sorted(Comparator.comparingInt(l -> -l.get(0).minCdrDepth))
                .sorted(Comparator.comparingInt(l -> -l.get(0).longPeptide))
                .sorted(Comparator.comparingDouble(l -> -l.get(0).score))
                .collect(Collectors.toList());
        Map<String, List<List<ProteinCoverage>>> sortedGroupByType = new HashMap<>();
        for (List<ProteinCoverage> group : sortedGroups) {
            String type = group.get(0).pName.split("-")[0];
            sortedGroupByType.computeIfAbsent(type, key -> new ArrayList<>()).add(group);
        }
        writeFasta(sortedGroupByType);
    }

    // read psms.csv and proteins.csv from novor.cloud
    private static List<String[]> readFile(String path, boolean skip1Line) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path));
        if (skip1Line) {
            br.readLine();
        }
        List<String[]> res = new ArrayList<>();
        for (String line = br.readLine(); line != null; line = br.readLine()) {
            line = line.replace("\n", "");
            String[] split = line.split(",");
            res.add(split);
        }
        return res;
    }

    // read proPepLinks.txt from novor.cloud
    private static Map<String, String[]> readProPepLinks(String path) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(path));
        Map<String, String[]> res = new HashMap<>();
        for (String line = br.readLine(); line != null; line = br.readLine()) {
            line = line.replaceAll("\\s", "").replace("\n", "");
            String[] split = line.split(":");
            res.put(split[0], split[1].split(","));
        }
        return res;
    }

    static class ProteinCoverage {
        String gId;
        String pName;
        String cdrSeq;
        String proteinSeq;
        int[] coverageDepth;
        int[] cdrRange;
        List<List<String[]>> proteinPeps;
        int minCdrDepth = Integer.MAX_VALUE;
        int longPeptide = 0;
        double cdr3Score = 0;
        double coverageScore = 0;
        int coverageLength = 0;
        double score;
        public ProteinCoverage(String gId, String pName, String cdrSeq, String proteinSeq, List<List<String[]>> proteinPeps) {
            this.gId = gId;
            this.pName = pName;
            this.proteinSeq = proteinSeq;
            this.cdrSeq = cdrSeq;
            int cdrS = proteinSeq.indexOf(cdrSeq);
            this.cdrRange = new int[]{cdrS, cdrS + cdrSeq.length()};
            this.coverageDepth = new int[proteinSeq.length()];
            this.proteinPeps = proteinPeps;
            initCoverage();
        }
        private void initCoverage() {
            // init coverageDepth, for each AA in the protein, how many different peptides cover it
            for (List<String[]> peps : this.proteinPeps) {
                int size = peps.size();
                String seq = ContigAssembly.removeCharactersInBrackets(peps.get(0)[9]);
                int pepStart = proteinSeq.indexOf(seq);
                for (int i = 0; i < seq.length(); i++) {
                    coverageDepth[i + pepStart] += size;
                }
                if (pepStart <= cdrRange[0] && pepStart + seq.length() >= cdrRange[1]) {
                    longPeptide++;
                }
            }
            // if cdr regions is fully covered
            for (int i = cdrRange[0]; i < cdrRange[1]; i++) {
                minCdrDepth = Math.min(minCdrDepth, coverageDepth[i]);
            }
            // if long peptide covers full cdr or connect both cdrs
            //calculate cdr3 score
            int[] cdrDepth = new int[cdrRange[1] - cdrRange[0]];
            System.arraycopy(coverageDepth, cdrRange[0], cdrDepth, 0, cdrDepth.length);
            Arrays.sort(cdrDepth);
            for (int i = 0; i < Math.min(8, cdrDepth.length); i++) {
                if (cdrDepth[i] == 0) {
                    cdr3Score -= 1;
                } else {
                    cdr3Score += normalizeDepthScore(cdrDepth[i]) / (i + 1);
                }
            }

            for (int i = 0; i < proteinSeq.length(); i++) {
                if (coverageDepth[i] > 0) {
                    coverageScore += normalizeDepthScore(coverageDepth[i]);
                    coverageLength++;
                }
            }
            score = coverageScore / 100 + cdr3Score ;
        }

        private double normalizeDepthScore(int depth) {
            if (depth <= 2) {
                return depth;
            } else {
                return 2 + Math.log(Math.min(8, depth) - 1);
            }
        }

    }


}
