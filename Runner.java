import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class Runner {
    public static void main(String[] args) {
        int readLength = 0;
        int frLengthMean = 0;
        int SD = 0;
        File readcounts = null;
        double mutationrate = 0;
        File fasta = null;
        File fidx = null;
        File gtf = null;
        File od = null;
        if (args.length != 18) {
            String consoleParameters = "";
            for (int i = 0; i < args.length; i++)
                consoleParameters = consoleParameters.concat(" args[" + i + "]=" + args[i]);
            System.err.println("Required 18 parameters. Found: " + args.length + " parameters: " + consoleParameters);
            System.exit(0);
        } else {
            for (int i = 0; i < args.length; i = i + 2) {
                switch (args[i]) {
                    case "-length":
                        readLength = Integer.parseInt(args[i + 1]);
                        continue;
                    case "-frlength":
                        frLengthMean = Integer.parseInt(args[i + 1]);
                        continue;
                    case "-SD":
                        SD = Integer.parseInt(args[i + 1]);
                        continue;
                    case "-readcounts":
                        readcounts = new File(args[i + 1]);
                        continue;
                    case "-mutationrate":
                        mutationrate = Double.parseDouble(args[i + 1]);
                        continue;
                    case "-fasta":
                        fasta = new File(args[i + 1]);
                        continue;
                    case "-fidx":
                        fidx = new File(args[i + 1]);
                        continue;
                    case "-gtf":
                        gtf = new File(args[i + 1]);
                        continue;
                    case "-od":
                        od = new File(args[i + 1]);
                }
            }
            System.out.println("|UserInput|\nreadLength:\t" + readLength + "\nfrLengthMean:\t" + frLengthMean + "\nSD:\t" + SD + "\nreadcounts:\t" + readcounts.toString() + "\nmutationrate:\t" + mutationrate + "\nfasta:\t" + fasta.toString() + "\nfidx:\t" + fidx.toString() + "\ngtf:\t" + gtf.toString() + "\nod:\t" + od.toString());
        }
        // Object[0] = String gene_id, Object[1] = String transcript_id, Object[2] = Integer count
        ArrayList<Object[]> allReadcounts = parseReadcounts(readcounts);

        GSE gse = new GSE(fasta, fidx, gtf);

        // For every ReadcountLine, there's one outputLine containing readsFw, readsRw, readMappingInfo
        ArrayList<OutputLines> outLines = new ArrayList<>();
        for (Object[] obj : allReadcounts) {
            String gene_id = (String) obj[0];
            String transcript_id = (String) obj[1];
            int counts = (Integer) obj[2];
            outLines.add(getOutputLines(gene_id, transcript_id, counts, gse, frLengthMean, SD, readLength, mutationrate));
        }


        Boolean append = false;
        int counter = 0;
        for (OutputLines line : outLines) {
            // append to files
            counter = writeOutputLines(line, od, append, counter);
            append = true;
        }

        // File gzip = new File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt2\\Homo_sapiens.GRCh37.75.cdna.all.fa.gz");
        // GSE.checkTranscript("ENST00000515242", gzip);
        System.out.println("End of main method.");
    }

    private static OutputLines getOutputLines(String gene_id, String transcript_id, int counts, GSE gse, int frLengthMean, int SD, int readLength, double mutationrate) {
        ArrayList<String> readsFw = new ArrayList<>();
        ArrayList<String> readsRw = new ArrayList<>();
        ArrayList<String> readMappingInfo = new ArrayList<>();

        // get transcript sequence
        TreeSet<Exon> exonsInTranscript = gse.getExonsByTranscriptAndGene(transcript_id, gene_id);
        String transcriptSequence = "";
        Exon referenceExon = null;
        for (Exon e : exonsInTranscript) {
            String s = gse.getSequence(e.chr, e.start, e.end, e.strand);
            transcriptSequence = transcriptSequence.concat(s);
            referenceExon = new Exon(e);
        }

        for (int j = 0; j < counts; j++) {
            // TODO
            // j = counts;
            // Get FL from max(readlength, ND(mean, SD))
            int fragmentLength = GSE.getFragmentLength(frLengthMean, SD, readLength, transcriptSequence.length());
            // random fragmentPosition from 0 to length(t) - FL
            int fragmentStart = new Random().nextInt(transcriptSequence.length() - fragmentLength);
            // get fragment sequence
            String fragmentSeq = transcriptSequence.substring(fragmentStart, fragmentStart + fragmentLength);
            // System.out.println("FragmentLength: " + fragmentLength + " fragmentPosition: " + fragmentStart + " fragmentSeq: " + fragmentSeq);
            // get read sequences of length readlength (second read is reverse complement)
            String readFw = GSE.getReadSequence(fragmentSeq, readLength, false);
            String readRw = GSE.getReadSequence(fragmentSeq, readLength, true);
            // simulate mutations with required rate
            String[] mutatedReadFw = GSE.simulateMutations(readFw, mutationrate);
            String[] mutatedReadRw = GSE.simulateMutations(readRw, mutationrate);

            readsFw.add(mutatedReadFw[0]);
            readsRw.add(mutatedReadRw[0]);

            // readid chr gene transcript t_fw_regvec t_rw_regvec fw_regvec rw_regvec fw_mut rw_mut 0 19 ENSG00000104870 ENST00000221466 810-885 854-929 50017390-50017391|50017468-50017542 50017511-50017586 55
            int tFwRegVecStart = fragmentStart - 1;
            int tFwRegVecEnd = fragmentStart + readLength - 1;
            int tRwRegVecStart = fragmentStart + fragmentLength - readLength - 1;
            int tRwRegVecEnd = fragmentStart + fragmentLength - 1;

            int fwRegVecStart = gse.getGenomicPosition(transcript_id, gene_id, fragmentStart);
            // System.out.println("getGenomicPosition for fwRegVecStart: " + (fragmentStart));
            int fwRegVecEnd = gse.getGenomicPosition(transcript_id, gene_id, (fragmentStart + readLength));
            // System.out.println("getGenomicPosition for fwRegVecEnd: " + (fragmentStart + readLength));
            int rwRegVecStart = gse.getGenomicPosition(transcript_id, gene_id, (fragmentStart + fragmentLength - readLength));
            // System.out.println("getGenomicPosition for rwRegVecStart: " + (fragmentStart + fragmentLength - readLength));
            int rwRegVecEnd = gse.getGenomicPosition(transcript_id, gene_id, (fragmentStart + fragmentLength));
            // System.out.println("getGenomicPosition for rwRegVecEnd: " + (fragmentStart + fragmentLength));

            String mappingInfo = referenceExon.chr + "\t" + referenceExon.gene_id + "\t" + referenceExon.transcript_id + "\t" + tFwRegVecStart + "-" + tFwRegVecEnd + "\t"
                    + tRwRegVecStart + "-" + tRwRegVecEnd + "\t" + fwRegVecStart + "-" + fwRegVecEnd + "\t" + rwRegVecStart + "-" + rwRegVecEnd + "\t" + mutatedReadFw[1] + "\t" + mutatedReadRw[1];
            readMappingInfo.add(mappingInfo);
            // System.out.println("readFw: " + readFw + "(" + readFw.length() + ") mutated to " + mutatedReadFw[0] + "(" + mutatedReadFw[0].length() + ")");
            // System.out.println("readRw: " + readRw + "(" + readRw.length() + ") mutated to " + mutatedReadRw[0] + "(" + mutatedReadRw[0].length() + ")");
            // System.out.println("transcriptSeq (" + transcriptSequence.length()/* + ") seq: " + transcriptSequence */);
        }
        return new OutputLines(readsFw, readsRw, readMappingInfo);
    }

    public static int writeOutputLines(OutputLines line, File od, Boolean append, int counter) {
        System.out.println("Begin: writeOutputLines. Append: " + append);
        // System.out.println("\tWriting fw.fastq...");
        try {
            int fwCounter = counter;
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(od + "/fw.fastq", append)));
            for (String string : line.readsFw) {
                out.println("@" + fwCounter);
                out.println(string);
                out.println("+" + fwCounter);
                // write string.length() times 'I'
                char[] chars = new char[string.length()];
                Arrays.fill(chars, 'I');
                out.println(new String(chars));
                fwCounter++;
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        // System.out.println("\tWriting rw.fastq...");
        try {
            int rwCounter = counter;
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(od + "/rw.fastq", append)));
            for (String string : line.readsRw) {
                out.println("@" + rwCounter);
                out.println(string);
                out.println("+" + rwCounter);
                // write string.length() times 'I'
                char[] chars = new char[string.length()];
                Arrays.fill(chars, 'I');
                out.println(new String(chars));
                rwCounter++;
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        // System.out.println("\tWriting read.mappinginfo...");
        int readCounter = counter;
        try {
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(od + "/read.mappinginfo", append)));
            if (!append)
                out.println("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut");
            for (String string : line.readMappingInfo) {
                out.println(readCounter + "\t" + string);
                readCounter++;
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return readCounter;
    }

    private static ArrayList<Object[]> parseReadcounts(File readcounts) {
        System.out.println("|parseReadcounts| Started...");
        ArrayList<Object[]> result = new ArrayList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(readcounts));
            String line = br.readLine(); // skip header
            int linecounter = 0; // for error printing only
            System.out.println("Begin: Parsing readcounts.");
            while ((line = br.readLine()) != null) {
                // gather readcounts line information
                Object[] obj = new Object[3];
                obj[0] = line.split("\\t")[0]; // gene_id
                obj[1] = line.split("\\t")[1]; // transcript_id
                obj[2] = Integer.parseInt(line.split("\\t")[2]); // counts
                System.out.println("Readcount line: " + linecounter + " gene:" + obj[0] + " transcript: " + obj[1] + " counts: " + obj[2]);
                linecounter++;
                result.add(obj.clone());
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("|parseReadcounts| Finished. Found " + result.size() + " entries.");
        return result;
    }


}
