package r_plotting;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.TreeMap;

public class MappingInfo {
    // -----This information needs to be extracted from the read.mappinginfo-file-----
    // K = #mutations, V = #reads
    TreeMap<Integer, Integer> mutations2reads;
    // K = #fragmentLength, V = #reads
    TreeMap<Integer, Integer> fragLength2reads;
    int allReads;
    int nonSplitReads;
    int nonSplitReadsWithoutMismatch;
    int splitReads;
    int splitReadsWithoutMismatch;
    int splitReadsWithoutMismatch5bp;

    private MappingInfo(TreeMap<Integer, Integer> mutations2reads, TreeMap<Integer, Integer> fragLength2reads, int allReads,
                        int nonSplitReads, int nonSplitReadsWithoutMismatch, int splitReads, int splitReadsWithoutMismatch,
                        int splitReadsWithoutMismatch5bp) {
        this.mutations2reads = mutations2reads;
        this.fragLength2reads = fragLength2reads;
        this.allReads = allReads;
        this.nonSplitReads = nonSplitReads;
        this.nonSplitReadsWithoutMismatch = nonSplitReadsWithoutMismatch;
        this.splitReads = splitReads;
        this.splitReadsWithoutMismatch = splitReadsWithoutMismatch;
        this.splitReadsWithoutMismatch5bp = splitReadsWithoutMismatch5bp;
    }

    static MappingInfo getMappingInfoFromFile(File file) {
        System.out.println("|getMappingInfoFromFile| Started...");
        TreeMap<Integer, Integer> mutations2reads = new TreeMap<>();
        TreeMap<Integer, Integer> fragLength2reads = new TreeMap<>();
        int allReads = 0;
        int nonSplitReads = 0;
        int nonSplitReadsWithoutMismatch = 0;
        int splitReads = 0;
        int splitReadsWithoutMismatch = 0;
        int splitReadsWithoutMismatch5bp = 0;

        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = br.readLine(); //skip Headers
            while ((line = br.readLine()) != null) {
                //readid	chr	gene	transcript	t_fw_regvec	t_rw_regvec	fw_regvec	rw_regvec	fw_mut	rw_mut
                // 9	19	ENSG00000104870	ENST00000221466	424-499	528-603	50015960-50016008|50016644-50016671	50016700-50016731|50017139-50017183	2,12	66
                String[] tabSeperated = line.split("\t", -1);
                String fw_regvec = tabSeperated[6]; // 50015960-50016008|50016644-50016671
                String rw_regvec = tabSeperated[7]; // 50016700-50016731|50017139-50017183
                String fw_mut = tabSeperated[8]; // 2,12
                String rw_mut = tabSeperated[9]; // 66

                // --- Gather Information for BoxPlots ---
                allReads++;
                if (!(fw_regvec.contains("|") || rw_regvec.contains("|"))) {
                    nonSplitReads++;
                    if (fw_mut.isEmpty() && rw_mut.isEmpty())
                        nonSplitReadsWithoutMismatch++;
                } else {
                    splitReads++;
                    if (fw_mut.isEmpty() && rw_mut.isEmpty()) {
                        splitReadsWithoutMismatch++;
                        // check if current read is a splitReadsWithoutMismatch5bp
                        ArrayList<Integer> regionLength = new ArrayList<>();
                        for (String region : fw_regvec.split("\\|")) {
                            int start = Integer.parseInt(region.substring(0, region.indexOf("-")));
                            int end = Integer.parseInt(region.substring(region.indexOf("-") + 1, region.length()));
                            regionLength.add(end - start);
                        }
                        for (String region : rw_regvec.split("\\|")) {
                            int start = Integer.parseInt(region.substring(0, region.indexOf("-")));
                            int end = Integer.parseInt(region.substring(region.indexOf("-") + 1, region.length()));
                            regionLength.add(end - start);
                        }
                        boolean containsRegionSmallerThan5bp = false;
                        for (int length : regionLength) {
                            if (length < 5)
                                containsRegionSmallerThan5bp = true;
                        }
                        if (!containsRegionSmallerThan5bp)
                            splitReadsWithoutMismatch5bp++;
                    }
                }
                // --- END: Gather Information for BoxPlots ---

                // --- Gather Information for Distribution Plots ---
                // count number of mutations
                int numberMutations = 0;
                if (!fw_mut.equals(""))
                    numberMutations += fw_mut.split(",").length;
                if (!rw_mut.equals(""))
                    numberMutations += rw_mut.split(",").length;

                if (mutations2reads.containsKey(numberMutations))
                    mutations2reads.put(numberMutations, mutations2reads.get(numberMutations) + 1);
                else
                    mutations2reads.put(numberMutations, 1);

                // count fragment length
                String t_fw_regvec = tabSeperated[4];
                String t_rw_regvec = tabSeperated[5];
                int start = Integer.parseInt(t_fw_regvec.substring(0, t_fw_regvec.indexOf("-")));
                int end = Integer.parseInt(t_rw_regvec.substring(t_rw_regvec.lastIndexOf("-") + 1, t_rw_regvec.length()));
                int fragmentLength = Math.abs(end - start); // absolute value

                if (fragLength2reads.containsKey(fragmentLength))
                    fragLength2reads.put(fragmentLength, fragLength2reads.get(fragmentLength) + 1);
                else
                    fragLength2reads.put(fragmentLength, 1);
                // --- END: Gather Information for Distribution Plots ---
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("|getMappingInfoFromFile| Finished. Result:");
        System.out.println("\tallReads:" + allReads + "\n\tnonSplitReads:" + nonSplitReads + "\n\tnonSplitReadsWithoutMismatch:"
                + nonSplitReadsWithoutMismatch + "\n\tsplitReads:" + splitReads + "\n\tsplitReadsWithoutMismatch:"
                + splitReadsWithoutMismatch + "\n\tsplitReadsWithoutMismatch5bp:" + splitReadsWithoutMismatch5bp
                + "\n\tmutations2reads size:" + mutations2reads.size() + "\n\tfragLength2reads size:" + fragLength2reads.size());
        return new MappingInfo(mutations2reads, fragLength2reads, allReads, nonSplitReads, nonSplitReadsWithoutMismatch,
                splitReads, splitReadsWithoutMismatch, splitReadsWithoutMismatch5bp);
    }
}
