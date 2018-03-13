import java.util.TreeSet;

public class Transcript {
    public String transcript_id;
    // String = exon Ids
    public TreeSet<Exon> exons = new TreeSet<>(new Exon());

    // Default Constructor
    public Transcript() {
        super();
    }

    // Constructor
    public Transcript(String transcript_id, TreeSet<Exon> exons) {
        super();
        this.exons = exons;
        this.transcript_id = transcript_id;
    }

    // Constructor for cloning a Transcript
    public Transcript(Transcript another) {
        //Deep copy of Exon list required
        for (Exon e : another.exons) {
            this.exons.add(new Exon(e));
        }
        this.transcript_id = another.transcript_id;
    }

    public static void printTranscriptInfo(Transcript t) {
        System.out.print("|printTranscriptInfo| transcript_id: " + t.transcript_id + " contains exons: ");
        for (Exon e : t.exons) {
            System.out.print(e.exon_id + " | ");
        }
        System.out.println();
    }

}
