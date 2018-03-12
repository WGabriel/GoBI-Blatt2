import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.TreeSet;

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

        if (args.length == 0) {
            System.err.println("No console parameters found.");
        } else if (args.length != 18) {
            System.err.println("Required 18 parameters. Found: " + args.length);
            String consoleParameters = "";
            for (int i = 0; i < args.length; i++) {
                consoleParameters += " args[" + i + "]=" + args[i];
            }
            System.err.println("Too few or many ConsoleParameters:" + consoleParameters);
            // Console parameters are correct!
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
                        continue;
                }
            }
            System.out.println("|UserInput|\nreadLength:\t" + readLength + "\nfrLengthMean:\t" + frLengthMean
                    + "\nSD:\t" + SD + "\nreadcounts:\t" + readcounts.toString() + "\nmutationrate:\t" + mutationrate
                    + "\nfasta:\t" + fasta.toString() + "\nfidx:\t" + fidx.toString() + "\ngtf:\t" + gtf.toString()
                    + "\nod:\t" + od.toString());
        }

        // Parse file
        GenomeSequenceExtractor gse = new GenomeSequenceExtractor(fasta, fidx, gtf);
        // Go through readcounts and retrieve information

        ArrayList<String> readsFw = new ArrayList<String>();
        ArrayList<String> readsRw = new ArrayList<String>();
        ArrayList<String> readMappingInfo = new ArrayList<String>();
        int readId = 0;

        try {
            BufferedReader br = new BufferedReader(new FileReader(readcounts));
            String line = null;
            int linecounter = 0; // for error printing only
            // int lineLimit = 30;
            System.out.println("Begin: Parsing readcounts.");
            while ((line = br.readLine()) != null) {
                // ignore first line
                if (linecounter != 0 /* && lineLimit > 0 */) {
                    // lineLimit--;
//					System.out.println("#####################");
                    // gather readcounts line information
                    String[] tabSeparated = line.split("\\t");
                    String gene = tabSeparated[0];
                    String transcript = tabSeparated[1];
                    int counts = Integer.parseInt(tabSeparated[2]);
                    System.out.println("Readcount line " + linecounter + ": gene:" + gene + " transcript: " + transcript
                            + " counts: " + counts);

                    // get transcript sequence
                    TreeSet<Exon> exonsInTranscript = gse.getExonsByTranscriptAndGene(transcript, gene);
                    String transcriptSequence = "";

                    Exon referenceExon = null;

                    for (Exon e : exonsInTranscript) {
                        String s = gse.getSequence(e.chr, e.start, e.end, e.strand);
                        transcriptSequence = transcriptSequence + s;
                        referenceExon = new Exon(e);
                    }

                    for (int j = 0; j < counts; j++) {
                        // TODO
                        // j = counts;
                        // Get FL from max(readlength, ND(mean, SD))
                        int fragmentLength = GenomeSequenceExtractor.getFragmentLength(frLengthMean, SD, readLength,
                                transcriptSequence.length());
                        // random fragmentPosition from 0 to length(t) - FL
                        int fragmentStart = new Random().nextInt(transcriptSequence.length() - fragmentLength);
                        // get fragment sequence
                        String fragmentSeq = transcriptSequence.substring(fragmentStart,
                                fragmentStart + fragmentLength);
//						System.out.println("FragmentLength: " + fragmentLength + " fragmentPosition: " + fragmentStart
//								+ " fragmentSeq: " + fragmentSeq);
                        // get read sequences of length readlength (second read is reverse complement)
                        String readFw = GenomeSequenceExtractor.getReadSequence(fragmentSeq, readLength, false);
                        String readRw = GenomeSequenceExtractor.getReadSequence(fragmentSeq, readLength, true);
                        // simulate mutations with required rate
                        String[] mutatedReadFw = GenomeSequenceExtractor.simulateMutations(readFw, mutationrate);
                        String[] mutatedReadRw = GenomeSequenceExtractor.simulateMutations(readRw, mutationrate);

                        readsFw.add(mutatedReadFw[0]);
                        readsRw.add(mutatedReadRw[0]);

                        // readid chr gene transcript t_fw_regvec t_rw_regvec fw_regvec rw_regvec fw_mut rw_mut
                        // 0 19 ENSG00000104870 ENST00000221466 810-885 854-929 50017390-50017391|50017468-50017542
                        // 50017511-50017586 55

                        int tFwRegVecStart = fragmentStart - 1;
                        int tFwRegVecEnd = fragmentStart + readLength - 1;
                        int tRwRegVecStart = fragmentStart + fragmentLength - readLength - 1;
                        int tRwRegVecEnd = fragmentStart + fragmentLength - 1;

                        int fwRegVecStart = gse.getGenomicPosition(transcript, gene, (fragmentStart));
                        // System.out.println("getGenomicPosition for fwRegVecStart: " + (fragmentStart));
                        int fwRegVecEnd = gse.getGenomicPosition(transcript, gene, (fragmentStart + readLength));
                        // System.out.println("getGenomicPosition for fwRegVecEnd: " + (fragmentStart + readLength));
                        int rwRegVecStart = gse.getGenomicPosition(transcript, gene,
                                (fragmentStart + fragmentLength - readLength));
                        // System.out.println("getGenomicPosition for rwRegVecStart: "
                        // + (fragmentStart + fragmentLength - readLength));
                        int rwRegVecEnd = gse.getGenomicPosition(transcript, gene, (fragmentStart + fragmentLength));
                        // System.out.println("getGenomicPosition for rwRegVecEnd: " + (fragmentStart +
                        // fragmentLength));

                        String mappingInfo = readId + "\t" + referenceExon.chr + "\t" + referenceExon.gene_id + "\t"
                                + referenceExon.transcript_id + "\t" + tFwRegVecStart + "-" + tFwRegVecEnd + "\t"
                                + tRwRegVecStart + "-" + tRwRegVecEnd + "\t" + fwRegVecStart + "-" + fwRegVecEnd + "\t"
                                + rwRegVecStart + "-" + rwRegVecEnd + "\t" + mutatedReadFw[1] + "\t" + mutatedReadRw[1];
                        readMappingInfo.add(mappingInfo);
                        readId++;

//						System.out.println("readFw: " + readFw + "(" + readFw.length() + ") mutated to "
//								+ mutatedReadFw[0] + "(" + mutatedReadFw[0].length() + ")");
//						System.out.println("readRw: " + readRw + "(" + readRw.length() + ") mutated to "
//								+ mutatedReadRw[0] + "(" + mutatedReadRw[0].length() + ")");
//						System.out.println(
//								"transcriptSeq (" + transcriptSequence.length()/* + ") seq: " + transcriptSequence */);
                    }
                }
                linecounter++;
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("Begin: fw.fastq file.");
        try {
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(od + "/fw.fastq", false)));
            int counter = 0;
            for (String string : readsFw) {
                out.println("@" + counter);
                out.println(string);
                out.println("+" + counter);
                // write string.length() times 'I'
                char[] chars = new char[string.length()];
                Arrays.fill(chars, 'I');
                out.println(new String(chars));
                counter++;
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.println("Begin: rw.fastq file.");
        try {
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(od + "/rw.fastq", false)));
            int counter = 0;
            for (String string : readsRw) {
                out.println("@" + counter);
                out.println(string);
                out.println("+" + counter);
                // write string.length() times 'I'
                char[] chars = new char[string.length()];
                Arrays.fill(chars, 'I');
                out.println(new String(chars));
                counter++;
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.println("Begin: read.mappinginfo file.");
        try {
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(od + "/read.mappinginfo", false)));
            out.println(
                    "readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut");
            for (String string : readMappingInfo) {
                out.println(string);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        // File gzip = new
        // File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt2\\Homo_sapiens.GRCh37.75.cdna.all.fa.gz");
        // GenomeSequenceExtractor.checkTranscript("ENST00000515242", gzip);

        System.out.println("End of main method.");
    }

}
