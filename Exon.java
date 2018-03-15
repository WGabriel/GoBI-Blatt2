import java.util.Comparator;

public class Exon implements Comparator<Exon> {
    public String chr;
    public int start;
    public int end;
    public String strand;
    public String exon_id;
    public int exon_number;
    public String transcript_id;
    public String transcript_name;
    public String gene_id;
    public String gene_name;

    // Constructor
    public Exon(String chr, int start, int end, String strand, String exon_id, int exon_number, String transcript_id,
                String transcript_name, String gene_id, String gene_name) {
        super();
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.exon_id = exon_id;
        this.exon_number = exon_number;
        this.transcript_id = transcript_id;
        this.transcript_name = transcript_name;
        this.gene_id = gene_id;
        this.gene_name = gene_name;
    }

    // Constructor for compare
    public Exon() {
        super();
    }

    // This constructor is used to create a deep copy of another exon
    public Exon(Exon another) {
        this.chr = another.chr;
        this.start = another.start;
        this.end = another.end;
        this.strand = another.strand;
        this.exon_id = another.exon_id;
        this.exon_number = another.exon_number;
        this.transcript_id = another.transcript_id;
        this.transcript_name = another.transcript_name;
        this.gene_id = another.gene_id;
        this.gene_name = another.gene_name;
    }

    public static void printExonInfo(Exon e) {
        System.out.println("|printExonInfo| exon_id: " + e.exon_id + " | transcript_id: " + e.transcript_id
                + " | transcript_name: " + e.transcript_name + " | gene_id: " + e.gene_id + " | gene_name: "
                + e.gene_name + " | chr: " + e.chr + " | start: " + e.start
                + " | end: " + e.end + " | strand: " + e.strand + " |");
    }

    @Override
    public int compare(Exon e1, Exon e2) {
        // Sorts Exons ascending to their start-points
        if (e1.exon_number > e2.exon_number) {
            return 1;
        } else if (e1.exon_number < e2.exon_number) {
            return -1;
        } else {
            if (!e1.exon_id.equals(e2.exon_id))
                System.err.println("Found two equally enumerated exons: " + e1.exon_id + " and " + e2.exon_id);
            return 0;
        }
    }

}
