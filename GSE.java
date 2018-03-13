import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Random;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

public class GSE {
    // Genome Sequence Extractor
    File fasta;
    File idx;
    File gtf;
    HashMap<String, Gene> genes = new HashMap<>();

    // Constructor
    public GSE(File fasta, File idx, File gtf) {
        super();
        this.fasta = fasta;
        this.idx = idx;
        this.gtf = gtf;
        this.genes = parseGtf(gtf);
    }


    public String getSequence(String chr, int seqStart, int seqEnd, String strand) {
        String result = "";
        // Start & Length of chromosome
        long chrStart = 0;
        long chrLength = 0;
        int lineSize = 0;
        // Go through idx file & fetch start and length for a chromosome
        try {
            BufferedReader br = new BufferedReader(new FileReader(idx));
            String line;
            while ((line = br.readLine()) != null) {
                String[] tabSeparated = line.split("\\t");
                if (tabSeparated[0].equals(chr)) {
                    // Chromosome was found, chrStart and chrLength are filled
                    chrLength = Long.parseLong(tabSeparated[1]);
                    chrStart = Long.parseLong(tabSeparated[2]);
                    lineSize = Integer.parseInt(tabSeparated[3]);
                }
            }
            br.close();
            if (chrStart == 0 || chrLength == 0 || lineSize == 0) {
                // Chromosome was not found.
                System.err.println("Required information for chr " + chr + " was not found in idx file.");
            } else if (seqEnd - seqStart > chrLength) {
                System.err.println("Required seq length is longer than chr itself.");
            } else {
                // ("r"=readOnly)
                RandomAccessFile raf = new RandomAccessFile(fasta, "r");

                // read chr from begin to end
                long rafStart = chrStart + seqStart + (seqStart / lineSize);
                int rafLength = seqEnd - seqStart;

                // Jump to rafStart location
                raf.seek(rafStart);
                // reads seqLength bytes from raf and save it to b
                // Note: occurring newline-characters will be added to bytes
                byte[] bytes = new byte[rafLength];
                raf.readFully(bytes);
                result = new String(bytes, StandardCharsets.UTF_8);
                // Delete newlines /n from temp result-string
                result = result.replaceAll("(\\r|\\n)", "");
                // re-run the reading algorithm with updated byte-length
                raf.seek(rafStart);
                int amountOfNewLineCharacters = rafLength - result.length();
                bytes = new byte[rafLength + amountOfNewLineCharacters];
                raf.readFully(bytes);
                result = new String(bytes, StandardCharsets.UTF_8);
                // Delete newlines /n from final result-string
                result = result.replaceAll("(\\r|\\n)", "");
                // System.out.println(result);
                raf.close();

                if (strand.equals("-")) {
                    // Inverse result string
                    result = new StringBuilder(result).reverse().toString();
                    result = invertSequence(result);
                }
            }
        } catch (IOException e) {
            System.out.println("IOException:");
            e.printStackTrace();
        }
        if (result.isEmpty()) {
            System.err.println("No such sequence on chr:" + chr + "(" + seqStart + "," + seqEnd + ")");
        }
//		System.out.println("|getSequence| Input chr" + chr + "(" + seqStart + "-" + seqEnd + ") strand: " + strand
//				+ ". Result string length: " + result.length() + "/" + (seqEnd - seqStart));
        return result;
    }

    public TreeSet<Exon> getExonsByTranscriptAndGene(String transcriptId, String geneId) {
        // Returns sorted TreeSet of Exons belonging to transcript-/geneId
        // System.out.println("Started getRegionVectorByTranscriptId.");
        TreeSet<Exon> exons = this.genes.get(geneId).transcripts.get(transcriptId).exons;
        if (exons.size() == 0) {
            System.err.println("No exon, which belongs to transcriptID" + transcriptId + "was found!");
            // System.out.println("Finished getRegionVectorByTranscriptId.");
            return null;
        } else {
            // System.out.println("Finished getRegionVectorByTranscriptId.");
            return exons;
        }
    }

    public static int getFragmentLength(int mean, int standardDeviation, int readLength, int transcriptLength) {
        // returns a value from a normal distribution with a given mean and standardDeviation
        Random r = new Random();
        double number = r.nextGaussian() * standardDeviation + mean;
        // re-draws number, if smaller then readLength or larger then transcriptLength
        while (number < readLength || number > transcriptLength)
            number = r.nextGaussian() * standardDeviation + mean;
        return (int) number;
    }

    public static String getReadSequence(String fragment, int readlength, boolean reverseComplement) {
        // returns ReadSequence of a certain length for fragment
        String result;
        if (!reverseComplement)
            result = fragment.substring(0, readlength);
        else {
            // second read is reverse complement
            result = fragment.substring(fragment.length() - readlength, fragment.length());
            result = new StringBuilder(result).reverse().toString();
            result = invertSequence(result);
        }
        return result;
    }

    public static String[] simulateMutations(String sequence, double mutationratePercent) {
        // get mutated sequence and positions of mutations as string
        // result[0] = mutatedSeq, result[1] = cs positions of mutated chars
        String[] result = {"", ""};

        // go through all chars in sequence. select a random number [0,100]
        // if this number <= mutationrate: mutate the current char
        StringBuilder mutatedSeq = new StringBuilder(sequence);
        for (int i = 0; i < sequence.length(); i++) {
            int randomNumber = new Random().nextInt(100);
            if (randomNumber <= mutationratePercent) {
                result[1] = result[1].concat(i + ",");
                if (mutatedSeq.charAt(i) == 'A') {
                    mutatedSeq.setCharAt(i, 'G');
                } else if (mutatedSeq.charAt(i) == 'T') {
                    mutatedSeq.setCharAt(i, 'C');
                } else if (mutatedSeq.charAt(i) == 'G') {
                    mutatedSeq.setCharAt(i, 'T');
                } else if (mutatedSeq.charAt(i) == 'C') {
                    mutatedSeq.setCharAt(i, 'A');
                } else {
                    System.err.println("|simulateMutations| Found inappropriate char: " + mutatedSeq.charAt(i) + " in "
                            + sequence + ".");
                }
            }
        }
        if (!result[1].isEmpty()) {
            // trimm comma
            result[1] = result[1].substring(0, result[1].length() - 1);
        }
        result[0] = new String(mutatedSeq);
        return result;
    }

    public int getGenomicPosition(String transcriptId, String geneId, int positionInTranscript) {
        // return the position on the gene of a position on a transcript
        // go through all exons and reduce positionInTranscript every time we pass an exon
        int result = 0;
        TreeSet<Exon> exons = this.genes.get(geneId).transcripts.get(transcriptId).exons;
        for (Exon exon : exons) {
            int exonLength = exon.end - exon.start;
            if (exonLength < positionInTranscript)
                positionInTranscript -= exonLength;
            else
                result = exon.start + positionInTranscript;
        }
        if (result <= 0)
            System.err.println("|getGenomicPosition| positionInTranscript " + positionInTranscript + " not found. geneId: " + geneId + " transcriptId: " + transcriptId);
        return result;
    }

    public static String invertSequence(String transcriptSeq) {
        // invert given Sequence: A<->T, G<->C
        char[] result = new char[transcriptSeq.length()];
        for (int i = 0; i < transcriptSeq.length(); i++) {
            if (transcriptSeq.charAt(i) == (char) 'A')
                result[i] = 'T';
            else if (transcriptSeq.charAt(i) == (char) 'T')
                result[i] = 'A';
            else if (transcriptSeq.charAt(i) == (char) 'G')
                result[i] = 'C';
            else if (transcriptSeq.charAt(i) == (char) 'C')
                result[i] = 'G';
            else
                System.err.println("Invalid character found in transcriptSeq (" + (transcriptSeq.charAt(i) + ")"));
        }
        return new String(result);
    }

    public static void checkTranscriptOccurenceInGzip(String transcriptId, File gzip) {
        try {
            // looking for lines starting with >transcriptId
            InputStream gzipStream = new GZIPInputStream(new FileInputStream(gzip));
            BufferedReader br = new BufferedReader(new InputStreamReader(gzipStream, "UTF-8"));
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    // System.out.println(line);
                    if (line.startsWith(">" + transcriptId))
                        System.out.println("Line: " + line);
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static HashMap<String, Gene> parseGtf(File gtfInput) {
        // Key is the always the respective ID
        HashMap<String, Gene> allGenes = new HashMap<>();
        int exonCounter = 0;
        int transcriptCounter = 0;
        int geneCounter = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(gtfInput));
            String line;
            int linecounter = 0; // for error printing only
            System.out.println("Begin: Parsing gtf file.");
            while ((line = br.readLine()) != null) {
                linecounter++;
                // ignore comments (beginning with "#")
                if (!line.startsWith("#")) {
                    // For every line in input
                    String[] tabSeparated = line.split("\\t");
                    String seqname = tabSeparated[0];
                    // String source = tabSeparated[1];
                    String feature = tabSeparated[2];
                    String start = tabSeparated[3];
                    String end = tabSeparated[4];
                    // String score = tabSeparated[5];
                    String strand = tabSeparated[6];
                    // String frame = tabSeparated[7];
                    String attribute = tabSeparated[8];
                    // -------For lines, which are exons:-------
                    if (feature.equalsIgnoreCase("exon")) {
                        exonCounter++;
                        // parameters needed to construct new Exon()
                        String exon_id = "";
                        String exon_number = "";
                        String transcript_id = "";
                        String transcript_name = "";
                        String gene_id = "";
                        String gene_name = "";
                        // gather parameters from String "attribute"
                        String[] attributeSeparated = attribute.split(";");
                        for (int i = 0; i < attributeSeparated.length; i++) {
                            if (attributeSeparated[i].contains("exon_id")) {
                                // get only value between quotation marks
                                exon_id = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            } else if (attributeSeparated[i].contains("transcript_id"))
                                transcript_id = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            else if (attributeSeparated[i].contains("transcript_name"))
                                transcript_name = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            else if (attributeSeparated[i].contains("exon_number"))
                                exon_number = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            else if (attributeSeparated[i].contains("gene_id"))
                                gene_id = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                            else if (attributeSeparated[i].contains("gene_name"))
                                gene_name = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1, attributeSeparated[i].lastIndexOf("\""));
                        }
                        Exon e = new Exon(seqname, Integer.parseInt(start) - 1, Integer.parseInt(end), strand, exon_id, Integer.parseInt(exon_number), transcript_id, transcript_name, gene_id, gene_name);
                        // Check, if all exon values are actually filled. Note: e.gene_name can be empty, is not checked
                        if (e.start == 0 || e.end == 0 || e.chr.isEmpty() || e.strand.isEmpty() || e.exon_id.isEmpty() || e.exon_number == 0 || e.transcript_id.isEmpty() || e.gene_id.isEmpty())
                            System.err.println("Exon in line " + linecounter + " has an empty value!");
                        // Create data structure allGenes, contains transcripts, contains exons
                        if (!allGenes.containsKey(e.gene_id)) {
                            geneCounter++;
                            allGenes.put(e.gene_id, new Gene());
                        }
                        if (!allGenes.get(e.gene_id).transcripts.containsKey(e.transcript_id)) {
                            transcriptCounter++;
                            allGenes.get(e.gene_id).transcripts.put(e.transcript_id, new Transcript());
                        }
                        exonCounter++;
                        allGenes.get(e.gene_id).transcripts.get(e.transcript_id).exons.add(new Exon(e));
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Finished: Parsing file.");
        System.out.println("\texonCounter: " + exonCounter);
        System.out.println("\ttranscriptCounter:" + transcriptCounter);
        System.out.println("\tgeneCounter:" + geneCounter);
        return allGenes;
    }

}
