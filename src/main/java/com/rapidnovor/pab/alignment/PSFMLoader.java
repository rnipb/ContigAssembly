package com.rapidnovor.pab.alignment;

import java.io.*;
import java.net.URL;
import java.util.*;

public class PSFMLoader {
    List<PSFM> psfms;
    int[][][] cdrRegions;
    int[] kPos;

    public PSFMLoader() {
        File psfmParam = getModelFile("human_psfm.txt");
        psfms = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(psfmParam))) {
            String aaList = br.readLine().replace("\n", "");
            String chainType = null;
            List<Map<Character, Double>> aaRelativeEntropyMatrix = new ArrayList<>();
            cdrRegions = new int[3][3][2];
            kPos = new int[3];
            for (String line = br.readLine(); line != null; line = br.readLine()) {
                line = line.replace("\n", "");
                if (line.startsWith("IG")) {
                    if (chainType != null) {
                        psfms.add(new PSFM(aaRelativeEntropyMatrix, chainType));
                    }
                    chainType = line;
                    aaRelativeEntropyMatrix = new ArrayList<>();
                    continue;
                } else if (line.contains("=")) {
                    String[] split = line.split("=");
                    int[] array = Arrays.stream(split[1].split(",")).mapToInt(Integer::valueOf).toArray();
                    switch (split[0]) {
                        case "kPos":
                            kPos = array;
                            break;
                        case "hcdr":
                            cdrRegions[0] = to2D(array);
                            break;
                        case "kcdr":
                            cdrRegions[1] = to2D(array);
                            break;
                        case "lcdr":
                            cdrRegions[2] = to2D(array);
                            break;
                    }
                    continue;
                }
                double[] array = Arrays.stream(line.split(",")).mapToDouble(Double::valueOf).toArray();
                Map<Character, Double> aaE = new HashMap<>();
                for (int i = 0; i < array.length; i++) {
                    aaE.put(aaList.charAt(i), array[i]);
                }
                aaRelativeEntropyMatrix.add(aaE);
            }
            psfms.add(new PSFM(aaRelativeEntropyMatrix, chainType));
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static int[][] to2D(int[] array) {
        int[][] arr2d = new int[array.length / 2][2];
        for (int i = 0; i < array.length; i++) {
            arr2d[i / 2][i % 2] = array[i];
        }
        return arr2d;
    }

    public List<PSFM> getPsfms() {
        return psfms;
    }

    public int[][][] getCdrRegions() {
        return cdrRegions;
    }

    public int[] getkPos() {
        return kPos;
    }

    private static File getModelFile(String resourceName) {
        try {
            URL url = PSFMLoader.class.getResource(resourceName);
            if (url == null) {
                System.out.println("No input species in database.");
            }
            String type = url.getProtocol();
            if ("file".equals(type)) {
                return new File(url.toURI());
            }
            if ("jar".equals(type)) {
                File file = File.createTempFile("psfm", null);
                file.deleteOnExit();
                InputStream in = PSFMLoader.class.getResourceAsStream(resourceName);
                OutputStream out = new FileOutputStream(file);
                byte[] bytes = new byte[8192];
                int nBytes;
                while ((nBytes = in.read(bytes)) > 0) {
                    out.write(bytes, 0, nBytes);
                }
                return file;
            }
        } catch (Exception e) {
            System.out.println(e);
        }
        return new File(resourceName);
    }
}
